from __future__ import annotations
import math
from typing import Sequence

from .core import BackboneClass, Oligonucleotide, ToxicityRecord


def aucroc_score(labels: list[int], scores: list[float]) -> float:
    """Area under the ROC curve via the trapezoidal rule."""
    if len(set(labels)) < 2:
        return 0.5
    pairs = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    tp = fp = 0
    prev_fp = prev_tp = 0
    auc = 0.0
    prev_score = None
    for score, label in pairs:
        if score != prev_score and prev_score is not None:
            auc += (fp - prev_fp) * (tp + prev_tp) / 2
            prev_fp, prev_tp = fp, tp
        if label == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score
    auc += (fp - prev_fp) * (tp + prev_tp) / 2
    return auc / (n_pos * n_neg)


def calibration_error(
    labels: list[int],
    scores: list[float],
    n_bins: int = 10,
) -> float:
    """Expected Calibration Error (ECE) over equal-width probability bins."""
    bin_size = 1.0 / n_bins
    ece = 0.0
    n = len(labels)
    for b in range(n_bins):
        low, high = b * bin_size, (b + 1) * bin_size
        indices = [i for i, s in enumerate(scores) if low <= s < high]
        if not indices:
            continue
        avg_conf = sum(scores[i] for i in indices) / len(indices)
        avg_acc = sum(labels[i] for i in indices) / len(indices)
        ece += (len(indices) / n) * abs(avg_conf - avg_acc)
    return ece


def model_comparison(
    labels: list[int],
    scores_a: list[float],
    scores_b: list[float],
) -> dict:
    """Compare two models by AUROC and ECE."""
    return {
        "model_a": {"auroc": aucroc_score(labels, scores_a), "ece": calibration_error(labels, scores_a)},
        "model_b": {"auroc": aucroc_score(labels, scores_b), "ece": calibration_error(labels, scores_b)},
    }


def endpoint_breakdown(
    records: list[ToxicityRecord],
) -> dict[str, dict[str, float]]:
    """Per-endpoint descriptive stats (mean, min, max, count)."""
    from collections import defaultdict

    buckets: dict[str, list[float]] = defaultdict(list)
    for rec in records:
        buckets[rec.endpoint.value].append(rec.value)
    result = {}
    for ep, vals in buckets.items():
        result[ep] = {
            "mean": sum(vals) / len(vals),
            "min": min(vals),
            "max": max(vals),
            "count": float(len(vals)),
        }
    return result


def feature_importance_summary(
    feature_names: list[str],
    importances: list[float],
) -> list[dict]:
    """Sort features by importance descending."""
    pairs = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
    return [{"feature": f, "importance": v} for f, v in pairs]


def sequence_complexity(sequence: str) -> float:
    """Linguistic complexity of a nucleotide sequence via Shannon entropy.

    Returns entropy in bits (max ~2 bits for perfectly uniform ATGC distribution).
    Higher values indicate more diverse/complex sequences.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    counts = {base: seq.count(base) for base in "ATGC"}
    n = len(seq)
    entropy = 0.0
    for c in counts.values():
        if c > 0:
            p = c / n
            entropy -= p * math.log2(p)
    return entropy


def backbone_risk_tier(backbone: BackboneClass) -> str:
    """Classify backbone chemistry by hepatotoxicity risk tier.

    Based on published case data:
      high  - PS (phosphorothioate), F2_ANA (2'-F arabino)
      medium - LNA, PNA
      low   - PO, PMO, MORPHOLINO
    """
    high = {BackboneClass.PS, BackboneClass.F2_ANA}
    medium = {BackboneClass.LNA, BackboneClass.PNA}
    if backbone in high:
        return "high"
    if backbone in medium:
        return "medium"
    return "low"


def population_toxicity_summary(
    oligos: Sequence[Oligonucleotide],
) -> dict[str, float]:
    """Dataset-level summary: GC stats, mean complexity, risk-tier distribution."""
    if not oligos:
        return {}
    gc_vals = [o.gc_content for o in oligos]
    complexities = [sequence_complexity(o.sequence) for o in oligos]
    tier_counts: dict[str, int] = {"high": 0, "medium": 0, "low": 0}
    for o in oligos:
        tier_counts[backbone_risk_tier(o.backbone)] += 1
    n = len(oligos)
    return {
        "mean_gc": sum(gc_vals) / n,
        "min_gc": min(gc_vals),
        "max_gc": max(gc_vals),
        "mean_complexity": sum(complexities) / n,
        "frac_high_risk": tier_counts["high"] / n,
        "frac_medium_risk": tier_counts["medium"] / n,
        "frac_low_risk": tier_counts["low"] / n,
    }
