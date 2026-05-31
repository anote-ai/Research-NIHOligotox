"""
Benchmark evaluation for OligoTox-Predict models.

Computes standardized metrics on each benchmark split:
  - Regression: RMSE, R², Spearman ρ (on log10 IC50)
  - Classification: AUROC, AUPRC, F1-macro (on toxicity class)
  - Per-endpoint breakdown + aggregate summary

Output: metrics CSV + markdown leaderboard table.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    mean_squared_error, r2_score, f1_score,
    classification_report,
)
import click
from rich.console import Console
from rich.table import Table

console = Console()

TOX_CLASSES = ["inactive", "low", "moderate", "high"]


def evaluate_regression(
    y_true: np.ndarray,
    y_pred: np.ndarray,
) -> dict:
    valid = ~(np.isnan(y_true) | np.isnan(y_pred))
    if valid.sum() < 5:
        return {"rmse": np.nan, "r2": np.nan, "spearman_r": np.nan, "n": int(valid.sum())}
    yt, yp = y_true[valid], y_pred[valid]
    return {
        "rmse": float(np.sqrt(mean_squared_error(yt, yp))),
        "r2": float(r2_score(yt, yp)),
        "spearman_r": float(spearmanr(yt, yp).correlation),
        "n": int(valid.sum()),
    }


def evaluate_classification(
    y_true_class: np.ndarray,
    y_pred_class: np.ndarray,
    y_prob: Optional[np.ndarray] = None,
) -> dict:
    valid = ~pd.isnull(y_true_class)
    yt, yp = y_true_class[valid], y_pred_class[valid]
    if len(yt) < 5:
        return {"auroc": np.nan, "auprc": np.nan, "f1_macro": np.nan, "n": int(valid.sum())}

    f1 = float(f1_score(yt, yp, average="macro", zero_division=0))
    result = {"f1_macro": f1, "n": int(valid.sum()), "auroc": np.nan, "auprc": np.nan}

    if y_prob is not None:
        yp_prob = y_prob[valid]
        try:
            result["auroc"] = float(roc_auc_score(yt, yp_prob, multi_class="ovr", average="macro"))
            result["auprc"] = float(average_precision_score(
                pd.get_dummies(yt).values, yp_prob, average="macro"
            ))
        except Exception:
            pass

    return result


def evaluate_model_on_split(
    split_df: pd.DataFrame,
    predictions_df: pd.DataFrame,
    endpoints: list[str],
    split_name: str = "test",
) -> pd.DataFrame:
    """
    Evaluate predictions against ground truth on the named split.

    Args:
        split_df: Has columns: oligo_id, split, {endpoint}_log10_ic50, {endpoint}_toxclass
        predictions_df: Has columns: oligo_id, {endpoint}_log10_ic50_pred, {endpoint}_toxclass_pred
        endpoints: List of toxicity endpoints to evaluate
        split_name: Which split to filter to ("test", "val", etc.)
    """
    test_df = split_df[split_df["split"] == split_name].copy()
    merged = test_df.merge(predictions_df, on="oligo_id", how="inner")
    if len(merged) == 0:
        console.print(f"[red]No predictions matched for split '{split_name}'[/red]")
        return pd.DataFrame()

    rows = []
    for endpoint in endpoints:
        true_col = f"{endpoint}_log10_ic50"
        pred_col = f"{endpoint}_log10_ic50_pred"
        true_cls_col = f"{endpoint}_toxclass"
        pred_cls_col = f"{endpoint}_toxclass_pred"

        reg_metrics = {}
        if true_col in merged.columns and pred_col in merged.columns:
            reg_metrics = evaluate_regression(
                merged[true_col].values, merged[pred_col].values
            )

        cls_metrics = {}
        if true_cls_col in merged.columns and pred_cls_col in merged.columns:
            cls_metrics = evaluate_classification(
                merged[true_cls_col].values, merged[pred_cls_col].values
            )

        rows.append({
            "endpoint": endpoint,
            "split": split_name,
            "n": reg_metrics.get("n", cls_metrics.get("n", 0)),
            "rmse": reg_metrics.get("rmse"),
            "r2": reg_metrics.get("r2"),
            "spearman_r": reg_metrics.get("spearman_r"),
            "auroc": cls_metrics.get("auroc"),
            "auprc": cls_metrics.get("auprc"),
            "f1_macro": cls_metrics.get("f1_macro"),
        })

    return pd.DataFrame(rows)


def aggregate_metrics(metrics_df: pd.DataFrame) -> dict:
    """Compute macro-average metrics across all endpoints."""
    return {
        "mean_rmse": float(metrics_df["rmse"].dropna().mean()),
        "mean_r2": float(metrics_df["r2"].dropna().mean()),
        "mean_spearman": float(metrics_df["spearman_r"].dropna().mean()),
        "mean_auroc": float(metrics_df["auroc"].dropna().mean()),
        "mean_f1": float(metrics_df["f1_macro"].dropna().mean()),
        "n_endpoints": len(metrics_df),
    }


def print_leaderboard(metrics_df: pd.DataFrame, model_name: str) -> None:
    table = Table(title=f"Evaluation Results — {model_name}")
    table.add_column("Endpoint", style="bold")
    table.add_column("N")
    table.add_column("RMSE")
    table.add_column("R²")
    table.add_column("Spearman ρ")
    table.add_column("AUROC")
    table.add_column("F1 (macro)")

    for _, row in metrics_df.iterrows():
        table.add_row(
            row["endpoint"],
            str(int(row["n"])) if not pd.isna(row["n"]) else "—",
            f"{row['rmse']:.3f}" if not pd.isna(row.get("rmse")) else "—",
            f"{row['r2']:.3f}" if not pd.isna(row.get("r2")) else "—",
            f"{row['spearman_r']:.3f}" if not pd.isna(row.get("spearman_r")) else "—",
            f"{row['auroc']:.3f}" if not pd.isna(row.get("auroc")) else "—",
            f"{row['f1_macro']:.3f}" if not pd.isna(row.get("f1_macro")) else "—",
        )
    console.print(table)

    agg = aggregate_metrics(metrics_df)
    console.print(f"\n[bold]Aggregate (macro-avg across endpoints):[/bold]")
    for k, v in agg.items():
        console.print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")


def generate_leaderboard_markdown(results: dict[str, pd.DataFrame]) -> str:
    """Generate a markdown leaderboard table from multiple model results."""
    lines = ["# OligoToxDB Benchmark Leaderboard\n"]
    lines.append("| Model | Mean RMSE | Mean R² | Mean Spearman ρ | Mean AUROC | Mean F1 |")
    lines.append("|---|---|---|---|---|---|")
    for model_name, metrics_df in results.items():
        agg = aggregate_metrics(metrics_df)
        lines.append(
            f"| {model_name} | {agg['mean_rmse']:.3f} | {agg['mean_r2']:.3f} | "
            f"{agg['mean_spearman']:.3f} | {agg['mean_auroc']:.3f} | {agg['mean_f1']:.3f} |"
        )
    return "\n".join(lines)


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.command()
@click.argument("split_csv", type=click.Path(exists=True))
@click.argument("predictions_csv", type=click.Path(exists=True))
@click.option("--endpoints-json", required=True)
@click.option("--split-name", default="test")
@click.option("--model-name", default="OligoTox-XGB")
@click.option("--output", default="evaluation_metrics.csv")
def cli(split_csv, predictions_csv, endpoints_json, split_name, model_name, output):
    """Evaluate model PREDICTIONS_CSV against ground truth SPLIT_CSV."""
    split_df = pd.read_csv(split_csv)
    pred_df = pd.read_csv(predictions_csv)
    with open(endpoints_json) as f:
        endpoints = json.load(f)

    metrics_df = evaluate_model_on_split(split_df, pred_df, endpoints, split_name)
    metrics_df.to_csv(output, index=False)
    print_leaderboard(metrics_df, model_name)
    console.print(f"\n[green]Metrics saved to {output}[/green]")


if __name__ == "__main__":
    cli()
