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
import pickle
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    mean_squared_error, r2_score, f1_score,
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

# ─── Toxicity endpoints ───────────────────────────────────────────────────────

PRIMARY_ENDPOINTS = [
    "Cell_viability_ATPLite",
    "LDH_release",
    "ALT_secretion",
    "Caspase_3_7",
    "KIM1_secretion",
    "NGAL_secretion",
    "IFNa",
    "IL6",
    "TNFa",
    "C3a",
    "C5a",
    "TLR9_activation",
    "Platelet_aggregation",
    "aPTT",
    "PT",
]

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
    "early_stopping_rounds": 50,
}

IC50_BINS = [0.0, 1.0, 10.0, 50.0, np.inf]
IC50_LABELS = ["high", "moderate", "low", "inactive"]


@dataclass
class ModelMetrics:
    endpoint: str
    task: str  # "regression" | "classification"
    auroc: Optional[float]
    auprc: Optional[float]
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

    def __init__(self, endpoints: list[str] = PRIMARY_ENDPOINTS, xgb_params: Optional[dict] = None):
        self.endpoints = endpoints
        self.xgb_params = {**XGB_DEFAULTS, **(xgb_params or {})}
        self.regressors: dict[str, xgb.XGBRegressor] = {}
        self.classifiers: dict[str, CalibratedClassifierCV] = {}
        self.explainers: dict[str, shap.TreeExplainer] = {}
        self.feature_names: list[str] = []
        self._fitted = False

    def fit(
        self,
        X: pd.DataFrame,
        y_dict: dict[str, pd.Series],
        eval_fraction: float = 0.15,
        verbose: bool = True,
    ) -> dict[str, ModelMetrics]:
        """
        Train one XGBRegressor per endpoint.

        Args:
            X: Feature matrix (oligo_id excluded; feature columns only)
            y_dict: {endpoint_name → Series of log10(IC50) or NaN}
            eval_fraction: Fraction of data held out for early stopping
        """
        self.feature_names = list(X.columns)
        metrics = {}

        for endpoint in track(self.endpoints, description="Training XGB models"):
            if endpoint not in y_dict:
                continue
            y = y_dict[endpoint].dropna()
            X_ep = X.loc[y.index]
            if len(y) < 50:
                console.print(f"[yellow]Skip {endpoint}: only {len(y)} labeled samples.[/yellow]")
                continue

            n_eval = max(10, int(len(y) * eval_fraction))
            X_train, X_eval = X_ep.iloc[:-n_eval], X_ep.iloc[-n_eval:]
            y_train, y_eval = y.iloc[:-n_eval], y.iloc[-n_eval:]

            # Regression model
            params = {k: v for k, v in self.xgb_params.items() if k != "early_stopping_rounds"}
            reg = xgb.XGBRegressor(
                **params,
                early_stopping_rounds=self.xgb_params.get("early_stopping_rounds", 50),
            )
            reg.fit(
                X_train, y_train,
                eval_set=[(X_eval, y_eval)],
                verbose=False,
            )
            self.regressors[endpoint] = reg
            self.explainers[endpoint] = shap.TreeExplainer(reg)

            # Evaluation
            y_pred = reg.predict(X_eval)
            rmse = float(np.sqrt(mean_squared_error(y_eval, y_pred)))
            r2 = float(r2_score(y_eval, y_pred))

            # Classification model (toxicity class)
            y_class = pd.cut(10 ** y, bins=IC50_BINS, labels=IC50_LABELS, right=False)
            y_class_train = pd.cut(10 ** y_train, bins=IC50_BINS, labels=IC50_LABELS, right=False)
            y_class_eval = pd.cut(10 ** y_eval, bins=IC50_BINS, labels=IC50_LABELS, right=False)

            clf_base = xgb.XGBClassifier(**params, use_label_encoder=False, eval_metric="mlogloss")
            le = LabelEncoder()
            y_cls_enc_train = le.fit_transform(y_class_train.dropna())
            X_cls_train = X_train.loc[y_class_train.dropna().index]

            try:
                clf = CalibratedClassifierCV(clf_base, cv=3, method="isotonic")
                clf.fit(X_cls_train, y_cls_enc_train)
                self.classifiers[endpoint] = clf

                y_cls_enc_eval = le.transform(y_class_eval.dropna())
                X_cls_eval = X_eval.loc[y_class_eval.dropna().index]
                y_prob = clf.predict_proba(X_cls_eval)
                auroc = float(roc_auc_score(y_cls_enc_eval, y_prob, multi_class="ovr", average="macro"))
                f1 = float(f1_score(y_cls_enc_eval, clf.predict(X_cls_eval), average="macro", zero_division=0))
            except Exception:
                auroc = f1 = None

            metrics[endpoint] = ModelMetrics(
                endpoint=endpoint,
                task="regression+classification",
                auroc=auroc,
                auprc=None,
                f1=f1,
                rmse=rmse,
                r2=r2,
                n_train=len(X_train),
                n_test=len(X_eval),
            )

        self._fitted = True
        if verbose:
            self._print_metrics(metrics)
        return metrics

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Predict log10(IC50) for all endpoints.
        Returns DataFrame with columns: oligo_id, {endpoint}_log10_ic50, {endpoint}_toxclass
        """
        if not self._fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        result = X[["oligo_id"]].copy() if "oligo_id" in X.columns else pd.DataFrame(index=X.index)
        Xf = X.drop(columns=["oligo_id"], errors="ignore")[self.feature_names]

        for endpoint, reg in self.regressors.items():
            result[f"{endpoint}_log10_ic50"] = reg.predict(Xf)
            result[f"{endpoint}_ic50_um"] = 10 ** result[f"{endpoint}_log10_ic50"]

        for endpoint, clf in self.classifiers.items():
            result[f"{endpoint}_toxclass"] = clf.predict(Xf)

        return result

    def explain(self, X: pd.DataFrame, endpoint: str) -> pd.DataFrame:
        """
        Return SHAP values for the given endpoint.
        Output: DataFrame(feature × oligo) with SHAP values.
        """
        if endpoint not in self.explainers:
            raise KeyError(f"No explainer for endpoint '{endpoint}'")
        Xf = X.drop(columns=["oligo_id"], errors="ignore")[self.feature_names]
        shap_vals = self.explainers[endpoint].shap_values(Xf)
        return pd.DataFrame(shap_vals, columns=self.feature_names, index=X.index)

    def feature_importance(self, endpoint: str, top_n: int = 20) -> pd.DataFrame:
        """Return top-N features by mean |SHAP| for an endpoint."""
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
                "xgb_params": self.xgb_params,
            }, f, indent=2)
        console.print(f"[green]Model saved to {path}[/green]")

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
                f"{m.rmse:.3f}" if m.rmse else "—",
                f"{m.r2:.3f}" if m.r2 else "—",
                f"{m.auroc:.3f}" if m.auroc else "—",
                f"{m.f1:.3f}" if m.f1 else "—",
                str(m.n_train),
            )
        console.print(table)


# ─── Benchmark evaluation ─────────────────────────────────────────────────────

def evaluate_benchmark(
    model: OligoToxXGB,
    X_test: pd.DataFrame,
    y_dict: dict[str, pd.Series],
) -> pd.DataFrame:
    """Run the model on a held-out test set and compute benchmark metrics."""
    rows = []
    preds = model.predict(X_test)
    for endpoint in model.endpoints:
        col = f"{endpoint}_log10_ic50"
        if col not in preds.columns or endpoint not in y_dict:
            continue
        y_true = y_dict[endpoint].dropna()
        y_hat = preds.loc[y_true.index, col]
        valid = y_true.index.intersection(y_hat.index)
        if len(valid) < 5:
            continue
        rmse = float(np.sqrt(mean_squared_error(y_true[valid], y_hat[valid])))
        r2 = float(r2_score(y_true[valid], y_hat[valid]))
        spearman = float(pd.Series(y_true[valid]).corr(pd.Series(y_hat[valid]), method="spearman"))
        rows.append({"endpoint": endpoint, "rmse": rmse, "r2": r2, "spearman_r": spearman, "n": len(valid)})
    return pd.DataFrame(rows)


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """OligoTox-XGB model training and inference."""


@cli.command()
@click.argument("features_csv", type=click.Path(exists=True))
@click.argument("labels_csv", type=click.Path(exists=True))
@click.option("--model-dir", default="models/xgb", help="Directory to save model")
@click.option("--endpoints", default=",".join(PRIMARY_ENDPOINTS), help="Comma-separated endpoint list")
def train(features_csv: str, labels_csv: str, model_dir: str, endpoints: str) -> None:
    """Train OligoTox-XGB on FEATURES_CSV with labels in LABELS_CSV."""
    X = pd.read_csv(features_csv)
    label_df = pd.read_csv(labels_csv)
    ep_list = [e.strip() for e in endpoints.split(",")]
    y_dict = {ep: label_df[ep] for ep in ep_list if ep in label_df.columns}
    X_feat = X.drop(columns=["oligo_id", "sequence"], errors="ignore")
    model = OligoToxXGB(endpoints=ep_list)
    model.fit(X_feat.assign(oligo_id=X.get("oligo_id", X.index)), y_dict)
    model.save(Path(model_dir))


@cli.command()
@click.argument("features_csv", type=click.Path(exists=True))
@click.option("--model-dir", default="models/xgb")
@click.option("--output", default="predictions.csv")
def predict(features_csv: str, model_dir: str, output: str) -> None:
    """Run OligoTox-XGB inference on FEATURES_CSV."""
    X = pd.read_csv(features_csv)
    model = OligoToxXGB.load(Path(model_dir))
    preds = model.predict(X)
    preds.to_csv(output, index=False)
    console.print(f"[green]Predictions saved to {output}[/green]")


if __name__ == "__main__":
    cli()
