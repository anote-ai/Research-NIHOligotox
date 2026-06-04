"""Tests for oligotox.evaluate."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pytest
from oligotox.core import BackboneClass, Oligonucleotide
from oligotox.evaluate import (
    aucroc_score,
    calibration_error,
    model_comparison,
    endpoint_breakdown,
    feature_importance_summary,
    sequence_complexity,
    backbone_risk_tier,
    population_toxicity_summary,
)


def test_aucroc_perfect():
    labels = [1, 1, 0, 0]
    scores = [0.9, 0.8, 0.2, 0.1]
    assert aucroc_score(labels, scores) == pytest.approx(1.0)


def test_aucroc_random_approx():
    labels = [1, 0, 1, 0]
    scores = [0.5, 0.5, 0.5, 0.5]
    assert aucroc_score(labels, scores) == pytest.approx(0.5)


def test_aucroc_all_same_score():
    labels = [1, 0]
    scores = [0.5, 0.5]
    assert 0.0 <= aucroc_score(labels, scores) <= 1.0


def test_calibration_error_perfect():
    labels = [1, 1, 0, 0]
    scores = [1.0, 1.0, 0.0, 0.0]
    assert calibration_error(labels, scores) == pytest.approx(0.0)


def test_model_comparison_returns_both_keys():
    labels = [1, 0, 1, 0]
    result = model_comparison(labels, [0.8, 0.2, 0.7, 0.3], [0.6, 0.4, 0.9, 0.1])
    assert "model_a" in result
    assert "model_b" in result


def test_endpoint_breakdown_structure():
    from oligotox.core import ToxicityEndpoint, ToxicityRecord
    records = [
        ToxicityRecord(oligo_id="o1", endpoint=ToxicityEndpoint.HEPATOTOXICITY, value=0.7),
        ToxicityRecord(oligo_id="o2", endpoint=ToxicityEndpoint.HEPATOTOXICITY, value=0.3),
    ]
    bd = endpoint_breakdown(records)
    assert "HEPATOTOXICITY" in bd
    assert bd["HEPATOTOXICITY"]["count"] == 2


def test_feature_importance_summary_sorted():
    names = ["gc_content", "length", "cpg_ratio"]
    imps = [0.3, 0.5, 0.2]
    result = feature_importance_summary(names, imps)
    assert result[0]["feature"] == "length"


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
