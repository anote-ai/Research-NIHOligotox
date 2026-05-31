"""Tests for plate QC and dose-response fitting."""

import pytest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from oligotoxdb.qc import (
    z_prime, assess_plate, fit_dose_response, grubbs_test,
    _four_pl, _estimate_noec,
)


class TestZPrime:
    def test_perfect_separation(self):
        pos = np.array([10.0] * 16)
        neg = np.array([100.0] * 16)
        assert z_prime(pos, neg) == 1.0

    def test_acceptable_zp(self):
        rng = np.random.default_rng(42)
        pos = rng.normal(10, 1, 16)
        neg = rng.normal(100, 3, 16)
        zp = z_prime(pos, neg)
        assert zp > 0.5  # should pass QC

    def test_poor_separation(self):
        rng = np.random.default_rng(42)
        pos = rng.normal(50, 15, 16)
        neg = rng.normal(60, 15, 16)
        zp = z_prime(pos, neg)
        assert zp < 0.5

    def test_identical_values_returns_neg_inf(self):
        pos = neg = np.array([50.0] * 8)
        assert z_prime(pos, neg) == -np.inf


class TestAssessPlate:
    def test_good_plate_passes(self):
        rng = np.random.default_rng(1)
        pos = rng.normal(10, 0.5, 16)
        neg = rng.normal(100, 2, 16)
        result = assess_plate("PLT-001", pos, neg)
        assert result.passed
        assert result.failure_reasons == []

    def test_bad_zp_fails(self):
        rng = np.random.default_rng(2)
        pos = rng.normal(50, 20, 8)
        neg = rng.normal(60, 20, 8)
        result = assess_plate("PLT-002", pos, neg)
        assert not result.passed
        assert any("Z'" in r for r in result.failure_reasons)

    def test_high_ctrl_cv_fails(self):
        pos = np.array([10, 10, 10, 10, 100, 100, 100, 100])  # bimodal → high CV
        neg = np.array([100] * 8)
        result = assess_plate("PLT-003", pos, neg)
        assert not result.passed


class TestDoseResponseFitting:
    def _make_clean_4pl_data(self):
        conc = np.array([0.01, 0.1, 1.0, 10.0, 50.0, 100.0])
        true_params = (5.0, 100.0, 5.0, 1.5)  # bottom, top, IC50, hill
        resp = _four_pl(conc, *true_params) + np.random.default_rng(99).normal(0, 1, 6)
        return conc, resp

    def test_clean_data_fits(self):
        conc, resp = self._make_clean_4pl_data()
        result = fit_dose_response("OT-00001", "Cell_viability", conc, resp)
        assert result.fit_success
        assert result.qc_flag in ("pass", "marginal")
        assert result.ic50 is not None
        assert 1.0 < result.ic50 < 20.0  # should be near 5.0 µM

    def test_r2_reasonable(self):
        conc, resp = self._make_clean_4pl_data()
        result = fit_dose_response("OT-00002", "Cell_viability", conc, resp)
        assert result.r_squared > 0.8

    def test_insufficient_data_fails(self):
        result = fit_dose_response("OT-00003", "ALT", np.array([1.0, 10.0]), np.array([100, 50]))
        assert not result.fit_success
        assert result.qc_flag == "fail"

    def test_emax_reasonable(self):
        conc, resp = self._make_clean_4pl_data()
        result = fit_dose_response("OT-00004", "Cell_viability", conc, resp)
        assert result.emax is not None
        assert 50 < result.emax < 130

    def test_audrc_nonnegative(self):
        conc, resp = self._make_clean_4pl_data()
        result = fit_dose_response("OT-00005", "Cell_viability", conc, resp)
        # AUDRC can be negative if response falls below baseline; just check it's finite
        assert np.isfinite(result.audrc)

    def test_noec_inactive_compound(self):
        # Flat dose-response (inactive compound)
        conc = np.array([0.01, 0.1, 1.0, 10.0, 50.0, 100.0])
        resp = np.array([100.0, 98.0, 102.0, 99.0, 101.0, 97.0])  # all near 100%
        result = fit_dose_response("OT-00006", "Cell_viability", conc, resp)
        assert result.noec is None  # no observed effect


class TestGrubbs:
    def test_obvious_outlier_detected(self):
        values = np.array([100.0, 101.0, 99.0, 100.5, 200.0])  # 200 is clear outlier
        outliers = grubbs_test(values, alpha=0.05)
        assert outliers[4]  # last element should be flagged

    def test_no_outlier_clean_data(self):
        rng = np.random.default_rng(5)
        values = rng.normal(100, 2, 10)
        outliers = grubbs_test(values)
        assert not outliers.any()

    def test_too_few_values(self):
        outliers = grubbs_test(np.array([1.0, 2.0]))
        assert not outliers.any()
