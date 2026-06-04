"""Tests for oligotox.evaluate."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oligotox.evaluate import (
    aucroc,
    calibration_error,
    model_comparison,
    endpoint_breakdown,
    feature_importance_summary,
    sequence_complexity,
    backbone_risk_tier,
    population_toxicity_summary,
)
from oligotox.core import BackboneClass, Oligonucleotide, ToxicityRecord, ToxicityEndpoint


def test_aucroc_perfect():
    scores = [0.9, 0.8, 0.2, 0.1]
    labels = [1, 1, 0, 0]
    assert abs(aucroc(scores, labels) - 1.0) < 1e-9


def test_aucroc_random_approx():
    import random
    rng = random.Random(0)
    scores = [rng.random() for _ in range(200)]
    labels = [rng.randint(0, 1) for _ in range(200)]
    auc = aucroc(scores, labels)
    assert 0.3 < auc < 0.7


def test_aucroc_all_same_score():
    scores = [0.5, 0.5, 0.5, 0.5]
    labels = [1, 0, 1, 0]
    assert aucroc(scores, labels) == 0.5


def test_calibration_error_perfect():
    probs2 = [0.0] * 10
    labels2 = [0] * 10
    ece = calibration_error(probs2, labels2, n_bins=10)
    assert ece == 0.0


def test_model_comparison_returns_both_keys():
    results = [
        {"model": "A", "aucroc": 0.9, "ece": 0.1},
        {"model": "B", "aucroc": 0.8, "ece": 0.05},
    ]
    out = model_comparison(results)
    assert "best_aucroc" in out and "best_ece" in out
    assert out["best_aucroc"] == "A"
    assert out["best_ece"] == "B"


def test_endpoint_breakdown_structure():
    records = [
        ToxicityRecord("o1", ToxicityEndpoint.HEPATOTOXICITY, 0.3, "PHH", 0.9),
        ToxicityRecord("o2", ToxicityEndpoint.HEPATOTOXICITY, 0.7, "PHH", 0.8),
        ToxicityRecord("o3", ToxicityEndpoint.NEPHROTOXICITY, 0.5, "KO", 0.7),
    ]
    out = endpoint_breakdown(records)
    assert "HEPATOTOXICITY" in out
    assert out["HEPATOTOXICITY"]["n"] == 2
    assert abs(out["HEPATOTOXICITY"]["mean_value"] - 0.5) < 1e-9


def test_feature_importance_summary_sorted():
    names = ["gc", "len", "cpg"]
    imps = [0.3, 0.5, 0.2]
    out = feature_importance_summary(names, imps)
    assert out[0]["feature"] == "len"
    assert out[-1]["feature"] == "cpg"


def test_sequence_complexity_uniform():
    score = sequence_complexity("ATGCATGC")
    assert score == pytest.approx(2.0, abs=0.01)


def test_sequence_complexity_single_base():
    score = sequence_complexity("AAAAAAA")
    assert score == pytest.approx(0.0)


def test_sequence_complexity_empty():
    assert sequence_complexity("") == 0.0


def test_backbone_risk_tier_ps_high():
    assert backbone_risk_tier(BackboneClass.PS) == "high"


def test_backbone_risk_tier_lna_medium():
    assert backbone_risk_tier(BackboneClass.LNA) == "medium"


def test_backbone_risk_tier_pmo_low():
    assert backbone_risk_tier(BackboneClass.PMO) == "low"


def test_population_toxicity_summary_keys():
    oligos = [
        Oligonucleotide(oligo_id="o1", sequence="ATGCATGCAT", backbone=BackboneClass.PS),
        Oligonucleotide(oligo_id="o2", sequence="GCTAGCTAGC", backbone=BackboneClass.PMO),
    ]
    summary = population_toxicity_summary(oligos)
    for key in ("mean_gc", "mean_complexity", "frac_high_risk", "frac_low_risk"):
        assert key in summary


def test_population_toxicity_summary_empty():
    assert population_toxicity_summary([]) == {}


import pytest
