"""Evaluation metrics for oligonucleotide toxicity prediction models."""

from __future__ import annotations

import math
from typing import Optional


def calibration_error(
    predicted_probs: list[float],
    true_labels: list[int],
    n_bins: int = 10,
) -> float:
    """Compute Expected Calibration Error (ECE).

    Bins predictions by confidence and measures mean absolute difference
    between average confidence and accuracy within each bin.

    Args:
        predicted_probs: Predicted probabilities in [0, 1].
        true_labels: Binary ground truth labels (0 or 1).
        n_bins: Number of equal-width bins.

    Returns:
        ECE scalar value in [0, 1].
    """
    if len(predicted_probs) != len(true_labels):
        raise ValueError("predicted_probs and true_labels must have the same length.")
    if not predicted_probs:
        return 0.0

    n = len(predicted_probs)
    bin_width = 1.0 / n_bins
    ece = 0.0

    for b in range(n_bins):
        low = b * bin_width
        high = low + bin_width
        indices = [
            i for i, p in enumerate(predicted_probs)
            if (low <= p < high) or (b == n_bins - 1 and p == 1.0)
        ]
        if not indices:
            continue
        avg_conf = sum(predicted_probs[i] for i in indices) / len(indices)
        avg_acc = sum(true_labels[i] for i in indices) / len(indices)
        ece += (len(indices) / n) * abs(avg_conf - avg_acc)

    return ece


def aucroc(scores: list[float], labels: list[int]) -> float:
    """Compute AUC-ROC using the trapezoidal rule via stdlib sort.

    Args:
        scores: Predicted scores/probabilities.
        labels: Binary ground truth labels (0 or 1).

    Returns:
        AUC-ROC value in [0, 1].

    Raises:
        ValueError: If labels contain only one class.
    """
    if len(scores) != len(labels):
        raise ValueError("scores and labels must have the same length.")
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        raise ValueError("labels must contain both positive and negative examples.")

    sorted_pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    tps = 0
    fps = 0
    auc = 0.0
    prev_fps = 0
    prev_tps = 0

    for _score, label in sorted_pairs:
        if label == 1:
            tps += 1
        else:
            fps += 1
        # Trapezoidal rule contribution
        auc += (fps - prev_fps) * (tps + prev_tps) / 2.0
        prev_fps = fps
        prev_tps = tps

    return auc / (n_pos * n_neg)


def model_comparison(results: list[dict]) -> dict:
    """Compare models and return the best by AUCROC and ECE.

    Args:
        results: List of dicts with keys: 'model_name', 'aucroc', 'ece'.

    Returns:
        Dict with 'best_by_aucroc' (model name) and 'best_by_ece' (model name).

    Raises:
        ValueError: If results list is empty.
    """
    if not results:
        raise ValueError("results list cannot be empty.")

    required_keys = {"model_name", "aucroc", "ece"}
    for r in results:
        missing = required_keys - set(r.keys())
        if missing:
            raise ValueError(f"Result dict missing keys: {missing}")

    best_aucroc = max(results, key=lambda r: r["aucroc"])
    best_ece = min(results, key=lambda r: r["ece"])

    return {
        "best_by_aucroc": best_aucroc["model_name"],
        "best_by_ece": best_ece["model_name"],
        "aucroc_value": best_aucroc["aucroc"],
        "ece_value": best_ece["ece"],
    }
