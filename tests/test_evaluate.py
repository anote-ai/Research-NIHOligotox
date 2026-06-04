"""Tests for oligotox.evaluate module."""

import pytest
from oligotox.evaluate import aucroc, calibration_error, model_comparison


def test_aucroc_perfect_predictor() -> None:
    scores = [0.9, 0.8, 0.2, 0.1]
    labels = [1, 1, 0, 0]
    assert aucroc(scores, labels) == pytest.approx(1.0)


def test_aucroc_random_predictor() -> None:
    # Random predictor should give ~0.5
    scores = [0.5, 0.5, 0.5, 0.5]
    labels = [1, 0, 1, 0]
    result = aucroc(scores, labels)
    assert 0.0 <= result <= 1.0


def test_aucroc_only_one_class_raises() -> None:
    with pytest.raises(ValueError):
        aucroc([0.9, 0.8], [1, 1])


def test_calibration_error_perfect_calibration() -> None:
    # When predictions match true labels exactly, ECE should be very small
    probs = [1.0] * 5 + [0.0] * 5
    labels = [1] * 5 + [0] * 5
    ece = calibration_error(probs, labels, n_bins=10)
    assert ece == pytest.approx(0.0, abs=1e-9)


def test_calibration_error_mismatched_lengths_raises() -> None:
    with pytest.raises(ValueError):
        calibration_error([0.5, 0.6], [1], n_bins=10)


def test_calibration_error_range() -> None:
    probs = [0.3, 0.7, 0.6, 0.4]
    labels = [0, 1, 0, 1]
    ece = calibration_error(probs, labels)
    assert 0.0 <= ece <= 1.0


def test_model_comparison_basic() -> None:
    results = [
        {"model_name": "XGB", "aucroc": 0.85, "ece": 0.12},
        {"model_name": "LR", "aucroc": 0.78, "ece": 0.05},
        {"model_name": "RF", "aucroc": 0.82, "ece": 0.08},
    ]
    comparison = model_comparison(results)
    assert comparison["best_by_aucroc"] == "XGB"
    assert comparison["best_by_ece"] == "LR"


def test_model_comparison_empty_raises() -> None:
    with pytest.raises(ValueError):
        model_comparison([])


def test_model_comparison_single_model() -> None:
    results = [{"model_name": "MLP", "aucroc": 0.75, "ece": 0.10}]
    comparison = model_comparison(results)
    assert comparison["best_by_aucroc"] == "MLP"
    assert comparison["best_by_ece"] == "MLP"
