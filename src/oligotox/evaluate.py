"""Evaluation utilities for OligoTox pipeline."""
from __future__ import annotations

from typing import List, Dict

from .core import ToxicityRecord


def aucroc(scores: List[float], labels: List[int]) -> float:
    """Compute AUC-ROC via trapezoidal rule (manual implementation)."""
    paired = sorted(zip(scores, labels), key=lambda x: -x[0])
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5

    tpr_list = [0.0]
    fpr_list = [0.0]
    tp = 0
    fp = 0
    prev_score = None
    for score, label in paired:
        if prev_score is not None and score != prev_score:
            tpr_list.append(tp / n_pos)
            fpr_list.append(fp / n_neg)
        if label == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score
    tpr_list.append(tp / n_pos)
    fpr_list.append(fp / n_neg)

    auc = 0.0
    for i in range(1, len(fpr_list)):
        auc += (fpr_list[i] - fpr_list[i - 1]) * (tpr_list[i] + tpr_list[i - 1]) / 2.0
    return float(auc)


def calibration_error(
    predicted_probs: List[float],
    true_labels: List[int],
    n_bins: int = 10,
) -> float:
    """Expected Calibration Error (ECE)."""
    bins = [[] for _ in range(n_bins)]
    bin_labels: List[List[int]] = [[] for _ in range(n_bins)]
    for p, l in zip(predicted_probs, true_labels):
        idx = min(int(p * n_bins), n_bins - 1)
        bins[idx].append(p)
        bin_labels[idx].append(l)

    n = len(predicted_probs)
    ece = 0.0
    for b_probs, b_labels in zip(bins, bin_labels):
        if not b_probs:
            continue
        mean_prob = sum(b_probs) / len(b_probs)
        frac_pos = sum(b_labels) / len(b_labels)
        ece += (len(b_probs) / n) * abs(mean_prob - frac_pos)
    return float(ece)


def model_comparison(results: List[Dict]) -> Dict:
    """Return best model by AUC-ROC and by ECE."""
    best_auc_model = max(results, key=lambda r: r["aucroc"])["model"]
    best_ece_model = min(results, key=lambda r: r["ece"])["model"]
    return {"best_aucroc": best_auc_model, "best_ece": best_ece_model}


def endpoint_breakdown(records: List[ToxicityRecord]) -> Dict:
    """Aggregate statistics per toxicity endpoint."""
    result: Dict[str, Dict] = {}
    for rec in records:
        key = rec.endpoint.value
        if key not in result:
            result[key] = {"n": 0, "sum_value": 0.0, "sum_confidence": 0.0}
        result[key]["n"] += 1
        result[key]["sum_value"] += rec.value
        result[key]["sum_confidence"] += rec.confidence
    return {
        k: {
            "n": v["n"],
            "mean_value": v["sum_value"] / v["n"],
            "mean_confidence": v["sum_confidence"] / v["n"],
        }
        for k, v in result.items()
    }


def feature_importance_summary(
    feature_names: List[str], importances: List[float]
) -> List[Dict]:
    """Return features sorted by importance descending."""
    paired = sorted(
        zip(feature_names, importances), key=lambda x: -x[1]
    )
    return [{"feature": name, "importance": imp} for name, imp in paired]
