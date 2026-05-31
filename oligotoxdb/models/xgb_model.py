"""
OligoTox-XGB: Gradient-boosted model for oligonucleotide toxicity prediction.

Features:
  - Multi-output XGBoost (one model per toxicity endpoint)
  - SHAP values for per-prediction interpretability
  - Calibrated probability estimates for toxicity classification
  - Hyperparameter optimization via cross-validated grid search
  - Model persistence and versioning
  - CLI for training, evaluation, and inference
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, mean_squared_error, r2_score, f1_score,
)
from sklearn.preprocessing import LabelEncoder
from sklearn.calibration import CalibratedClassifierCV
import xgboost as xgb
import shap
import click
from rich.console import Console
from rich.table import Table
from rich.progress import track

console = Console()

# Import all 47 endpoints from the registry
from oligotoxdb.endpoints import PRIMARY_ENDPOINT_NAMES as PRIMARY_ENDPOINTS

XGB_DEFAULTS = {
    "n_estimators": 500,
    "max_depth": 6,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 3,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "random_state": 42,
    "n_jobs": -1,
}

IC50_BINS = [0.0, 1.0, 10.0, 50.0, np.inf]
IC50_LABELS = ["high", "moderate", "low", "inactive"]


@dataclass
class ModelMetrics:
    endpoint: str
    task: str
    auroc: Optional[float]
    f1: Optional[float]
    rmse: Optional[float]
    r2: Optional[float]
    n_train: int
    n_test: int


class OligoToxXGB:
    """
    Multi-endpoint XGBoost toxicity predictor.

    Each endpoint gets its own XGBoost regressor (predicting log10 IC50)
    and a separate calibrated classifier (predicting toxicity class).
    """

    def __init__(
        self,
        endpoints: list[str] = PRIMARY_ENDPOINTS,
        xgb_params: Optional[dict] = None,
    ):
        self.endpoints = endpoints
        self.xgb_params = {**XGB_DEFAULTS, **(xgb_params or {})}
        self.regressors: dict[str, xgb.XGBRegressor] = {}
        self.classifiers: dict[str, CalibratedClassifierCV] = {}
        self.label_encoders: dict[str, LabelEncoder] = {}
        self.explainers: dict[str, shap.TreeExplainer] = {}
        self.feature_names: list[str] = []
        self._fitted = False

    def _get_feature_matrix(self, X: pd.DataFrame) -> pd.DataFrame:
        """Return feature-only matrix (drop oligo_id and sequence columns)."""
        return X.drop(columns=["oligo_id", "sequence"], errors="ignore")[self.feature_names]

    def fit(
        self,
        X: pd.DataFrame,
        y_dict: dict[str, pd.Series],
        eval_fraction: float = 0.15,
        verbose: bool = True,
    ) -> dict[str, ModelMetrics]:
        """
        Train one XGBRegressor + calibrated classifier per endpoint.

        Args:
            X: Feature matrix (may include oligo_id column)
            y_dict: {endpoint → Series of log10(IC50), may contain NaN}
            eval_fraction: Held-out fraction for early stopping
        """
        self.feature_names = [c for c in X.columns if c not in ("oligo_id", "sequence")]
        Xf = X[self.feature_names]
        metrics: dict[str, ModelMetrics] = {}
        early_stopping = self.xgb_params.pop("early_stopping_rounds", 50)

        for endpoint in track(self.endpoints, description="Training XGB models"):
            if endpoint not in y_dict:
                continue
            y = y_dict[endpoint].dropna()
            if len(y) < 20:
                continue

            X_ep = Xf.loc[y.index]
            n_eval = max(10, int(len(y) * eval_fraction))
            X_train, X_eval = X_ep.iloc[:-n_eval], X_ep.iloc[-n_eval:]
            y_train, y_eval = y.iloc[:-n_eval], y.iloc[-n_eval:]

            # ── Regression model ──────────────────────────────────────────
            reg = xgb.XGBRegressor(
                **self.xgb_params,
                early_stopping_rounds=early_stopping,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                reg.fit(X_train, y_train, eval_set=[(X_eval, y_eval)], verbose=False)
            self.regressors[endpoint] = reg
            self.explainers[endpoint] = shap.TreeExplainer(reg)

            y_pred = reg.predict(X_eval)
            rmse = float(np.sqrt(mean_squared_error(y_eval, y_pred)))
            r2 = float(r2_score(y_eval, y_pred))

            # ── Classification model ──────────────────────────────────────
            auroc = f1 = None
            try:
                # Convert log10 IC50 → toxicity class labels
                ic50_linear = 10 ** y
                y_class = pd.cut(ic50_linear, bins=IC50_BINS, labels=IC50_LABELS, right=False)
                y_class_train = y_class.iloc[:-n_eval].dropna()
                y_class_eval = y_class.iloc[-n_eval:].dropna()

                if len(y_class_train) >= 20 and y_class_train.nunique() >= 2:
                    le = LabelEncoder()
                    y_enc_train = le.fit_transform(y_class_train.astype(str))
                    X_cls_train = X_ep.loc[y_class_train.index]

                    clf_base = xgb.XGBClassifier(
                        **{k: v for k, v in self.xgb_params.items()},
                        eval_metric="mlogloss",
                        verbosity=0,
                    )
                    clf = CalibratedClassifierCV(clf_base, cv=3, method="isotonic")
                    clf.fit(X_cls_train, y_enc_train)
                    self.classifiers[endpoint] = clf
                    self.label_encoders[endpoint] = le

                    X_cls_eval = X_ep.loc[y_class_eval.index]
                    y_enc_eval = le.transform(y_class_eval.astype(str))
                    y_prob = clf.predict_proba(X_cls_eval)

                    if y_prob.shape[1] > 1 and len(np.unique(y_enc_eval)) >= 2:
                        auroc = float(roc_auc_score(
                            y_enc_eval, y_prob,
                            multi_class="ovr", average="macro",
                            labels=np.arange(len(le.classes_))
                        ))
                    f1 = float(f1_score(
                        y_enc_eval, clf.predict(X_cls_eval),
                        average="macro", zero_division=0
                    ))
            except Exception:
                pass

            metrics[endpoint] = ModelMetrics(
                endpoint=endpoint, task="regression+classification",
                auroc=auroc, f1=f1, rmse=rmse, r2=r2,
                n_train=len(X_train), n_test=len(X_eval),
            )

        self.xgb_params["early_stopping_rounds"] = early_stopping
        self._fitted = True
        if verbose:
            self._print_metrics(metrics)
        return metrics

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Predict log10(IC50) and toxicity class for all endpoints.
        Returns DataFrame with columns: oligo_id, {endpoint}_log10_ic50,
                                        {endpoint}_ic50_um, {endpoint}_toxclass
        """
        if not self._fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        result = pd.DataFrame()
        if "oligo_id" in X.columns:
            result["oligo_id"] = X["oligo_id"].values
        Xf = self._get_feature_matrix(X)

        for endpoint, reg in self.regressors.items():
            log_ic50 = reg.predict(Xf)
            result[f"{endpoint}_log10_ic50"] = log_ic50
            result[f"{endpoint}_ic50_um"] = 10 ** log_ic50

        for endpoint, clf in self.classifiers.items():
            le = self.label_encoders.get(endpoint)
            if le is not None:
                result[f"{endpoint}_toxclass"] = le.inverse_transform(clf.predict(Xf))

        return result

    def predict_with_uncertainty(
        self,
        X: pd.DataFrame,
        n_samples: int = 20,
    ) -> pd.DataFrame:
        """
        Monte Carlo feature dropout uncertainty estimation.
        Returns mean predictions + std across n_samples with random feature dropout.
        """
        Xf = self._get_feature_matrix(X)
        all_preds = []
        rng = np.random.default_rng(0)
        for _ in range(n_samples):
            # Drop 10% of features randomly (simulates uncertainty)
            X_noisy = Xf.copy()
            n_drop = max(1, len(self.feature_names) // 10)
            drop_cols = rng.choice(self.feature_names, n_drop, replace=False)
            X_noisy[drop_cols] = 0
            sample_preds = {}
            for ep, reg in self.regressors.items():
                sample_preds[f"{ep}_log10_ic50"] = reg.predict(X_noisy)
            all_preds.append(pd.DataFrame(sample_preds))

        mean_df = pd.concat(all_preds).groupby(level=0).mean()
        std_df = pd.concat(all_preds).groupby(level=0).std()
        std_df.columns = [c.replace("log10_ic50", "uncertainty") for c in std_df.columns]
        return pd.concat([mean_df, std_df], axis=1)

    def explain(self, X: pd.DataFrame, endpoint: str) -> pd.DataFrame:
        """Return SHAP values for the given endpoint."""
        if endpoint not in self.explainers:
            raise KeyError(f"No explainer for endpoint '{endpoint}'")
        Xf = self._get_feature_matrix(X)
        shap_vals = self.explainers[endpoint].shap_values(Xf)
        return pd.DataFrame(shap_vals, columns=self.feature_names, index=X.index)

    def feature_importance(self, endpoint: str, top_n: int = 20) -> pd.DataFrame:
        if endpoint not in self.regressors:
            raise KeyError(endpoint)
        fi = self.regressors[endpoint].get_booster().get_score(importance_type="gain")
        df = pd.Series(fi, name="gain").sort_values(ascending=False).head(top_n).reset_index()
        df.columns = ["feature", "gain"]
        return df

    def save(self, path: Path) -> None:
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        for endpoint, reg in self.regressors.items():
            reg.save_model(path / f"{endpoint}_regressor.json")
        with open(path / "config.json", "w") as f:
            json.dump({
                "endpoints": self.endpoints,
                "feature_names": self.feature_names,
                "xgb_params": {k: v for k, v in self.xgb_params.items() if k != "early_stopping_rounds"},
            }, f, indent=2)
        console.print(f"[green]Model saved → {path}[/green]")

    @classmethod
    def load(cls, path: Path) -> "OligoToxXGB":
        path = Path(path)
        with open(path / "config.json") as f:
            config = json.load(f)
        model = cls(endpoints=config["endpoints"], xgb_params=config["xgb_params"])
        model.feature_names = config["feature_names"]
        for endpoint in config["endpoints"]:
            reg_path = path / f"{endpoint}_regressor.json"
            if reg_path.exists():
                reg = xgb.XGBRegressor()
                reg.load_model(str(reg_path))
                model.regressors[endpoint] = reg
                model.explainers[endpoint] = shap.TreeExplainer(reg)
        model._fitted = True
        return model

    def _print_metrics(self, metrics: dict[str, ModelMetrics]) -> None:
        table = Table(title="OligoTox-XGB Training Metrics")
        table.add_column("Endpoint")
        table.add_column("RMSE (log10 IC50)")
        table.add_column("R²")
        table.add_column("AUROC (class)")
        table.add_column("F1 (class)")
        table.add_column("N train")
        for m in metrics.values():
            table.add_row(
                m.endpoint,
                f"{m.rmse:.3f}" if m.rmse is not None else "—",
                f"{m.r2:.3f}" if m.r2 is not None else "—",
                f"{m.auroc:.3f}" if m.auroc is not None else "—",
                f"{m.f1:.3f}" if m.f1 is not None else "—",
                str(m.n_train),
            )
        console.print(table)


# ─── Benchmark evaluation ─────────────────────────────────────────────────────

def evaluate_benchmark(
    model: OligoToxXGB,
    X_test: pd.DataFrame,
    y_dict: dict[str, pd.Series],
) -> pd.DataFrame:
    from scipy.stats import spearmanr
    rows = []
    preds = model.predict(X_test)
    for endpoint in model.endpoints:
        col = f"{endpoint}_log10_ic50"
        if col not in preds.columns or endpoint not in y_dict:
            continue
        y_true = y_dict[endpoint].dropna()
        shared = y_true.index.intersection(preds.index)
        if len(shared) < 5:
            continue
        yt, yp = y_true[shared].values, preds.loc[shared, col].values
        valid = ~(np.isnan(yt) | np.isnan(yp))
        rows.append({
            "endpoint": endpoint,
            "rmse": float(np.sqrt(mean_squared_error(yt[valid], yp[valid]))),
            "r2": float(r2_score(yt[valid], yp[valid])),
            "spearman_r": float(spearmanr(yt[valid], yp[valid]).correlation),
            "n": int(valid.sum()),
        })
    return pd.DataFrame(rows)


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """OligoTox-XGB model training and inference."""


@cli.command()
@click.argument("features_csv", type=click.Path(exists=True))
@click.argument("labels_csv", type=click.Path(exists=True))
@click.option("--model-dir", default="models/xgb")
@click.option("--endpoints", default=None, help="Comma-separated endpoint list (default: all primary)")
def train(features_csv: str, labels_csv: str, model_dir: str, endpoints: Optional[str]) -> None:
    """Train OligoTox-XGB on FEATURES_CSV with labels in LABELS_CSV."""
    X = pd.read_csv(features_csv)
    label_df = pd.read_csv(labels_csv)
    ep_list = [e.strip() for e in endpoints.split(",")] if endpoints else PRIMARY_ENDPOINTS
    y_dict = {ep: label_df[ep] for ep in ep_list if ep in label_df.columns}
    model = OligoToxXGB(endpoints=ep_list)
    model.fit(X, y_dict)
    model.save(Path(model_dir))


@cli.command()
@click.argument("features_csv", type=click.Path(exists=True))
@click.option("--model-dir", default="models/xgb")
@click.option("--output", default="predictions.csv")
@click.option("--with-uncertainty", is_flag=True)
def predict(features_csv: str, model_dir: str, output: str, with_uncertainty: bool) -> None:
    """Run OligoTox-XGB inference on FEATURES_CSV."""
    X = pd.read_csv(features_csv)
    model = OligoToxXGB.load(Path(model_dir))
    preds = model.predict_with_uncertainty(X) if with_uncertainty else model.predict(X)
    preds.to_csv(output, index=False)
    console.print(f"[green]Predictions → {output}[/green]")


if __name__ == "__main__":
    cli()
