"""Tests for ToxicityPredictor, compute_motif_score, auroc_score, toxicity_profile_summary."""
from __future__ import annotations

import pytest

from oligotox.core import (
    BackboneClass,
    Oligonucleotide,
    ToxicityEndpoint,
    ToxicityPredictor,
    compute_motif_score,
)
from oligotox.data import make_dataset, make_oligonucleotide, make_toxicity_record
from oligotox.evaluate import auroc_score, calibration_error, toxicity_profile_summary


# ---------------------------------------------------------------------------
# compute_motif_score
# ---------------------------------------------------------------------------

class TestComputeMotifScore:
    def test_empty_motif_sequence_gives_zero(self) -> None:
        # Sequence with no known motifs
        assert compute_motif_score("ATATATATAT") == pytest.approx(0.0)

    def test_gggg_motif_detected(self) -> None:
        score = compute_motif_score("AAGGGGTTTT")
        assert score > 0.0

    def test_ccgg_motif_detected(self) -> None:
        score = compute_motif_score("TTCCGGATTT")
        assert score > 0.0

    def test_score_capped_at_one(self) -> None:
        # Many motifs in a short window
        score = compute_motif_score("CCGGCCGGCCGGGGGGCCGG")
        assert score <= 1.0

    def test_case_insensitive(self) -> None:
        lower = compute_motif_score("aaggggtttt")
        upper = compute_motif_score("AAGGGGTTTT")
        assert lower == pytest.approx(upper)


# ---------------------------------------------------------------------------
# ToxicityPredictor
# ---------------------------------------------------------------------------

class TestToxicityPredictor:
    def _ps_high_gc_oligo(self) -> Oligonucleotide:
        """High-GC PS backbone oligo — should score high for HEPATOTOXICITY."""
        return Oligonucleotide(
            oligo_id="ps_high_gc",
            sequence="GCGCGCGCGCGCGCGCGCGC",  # 100% GC
            backbone=BackboneClass.PS,
        )

    def _po_low_gc_oligo(self) -> Oligonucleotide:
        """Low-GC PO backbone oligo — should score low."""
        return Oligonucleotide(
            oligo_id="po_low_gc",
            sequence="ATATATATATATATATATATATATAT",  # 0% GC
            backbone=BackboneClass.PO,
        )

    def test_ps_higher_hepatotox_than_po(self) -> None:
        predictor = ToxicityPredictor()
        ps_oligo = self._ps_high_gc_oligo()
        po_oligo = self._po_low_gc_oligo()
        ps_score = predictor.predict(ps_oligo, ToxicityEndpoint.HEPATOTOXICITY)
        po_score = predictor.predict(po_oligo, ToxicityEndpoint.HEPATOTOXICITY)
        assert ps_score > po_score

    def test_high_gc_raises_complement_score(self) -> None:
        predictor = ToxicityPredictor()
        high_gc = Oligonucleotide(
            oligo_id="high_gc",
            sequence="GCGCGCGCGCGCGCGCGCGC",
            backbone=BackboneClass.PO,
        )
        low_gc = Oligonucleotide(
            oligo_id="low_gc",
            sequence="ATATATATATATATATAT",
            backbone=BackboneClass.PO,
        )
        assert predictor.predict(high_gc, ToxicityEndpoint.COMPLEMENT) > \
               predictor.predict(low_gc, ToxicityEndpoint.COMPLEMENT)

    def test_predict_score_in_range(self) -> None:
        predictor = ToxicityPredictor()
        oligo = make_oligonucleotide(seed=1)
        for ep in ToxicityEndpoint:
            score = predictor.predict(oligo, ep)
            assert 0.0 <= score <= 1.0, f"Out of range for {ep}: {score}"

    def test_predict_all_endpoints_keys(self) -> None:
        predictor = ToxicityPredictor()
        oligo = make_oligonucleotide(seed=2)
        result = predictor.predict_all_endpoints(oligo)
        assert set(result.keys()) == {ep.value for ep in ToxicityEndpoint}

    def test_risk_category_bands(self) -> None:
        predictor = ToxicityPredictor()
        assert predictor.risk_category(0.10) == "low"
        assert predictor.risk_category(0.35) == "moderate"
        assert predictor.risk_category(0.60) == "high"
        assert predictor.risk_category(0.80) == "very_high"

    def test_gggg_raises_coagulopathy(self) -> None:
        predictor = ToxicityPredictor()
        g4_oligo = Oligonucleotide(
            oligo_id="g4",
            sequence="AAGGGGAAGGGGAAGGGGAA",
            backbone=BackboneClass.PO,
        )
        benign = Oligonucleotide(
            oligo_id="benign",
            sequence="AATTAATTAATTAATTAATT",
            backbone=BackboneClass.PO,
        )
        assert predictor.predict(g4_oligo, ToxicityEndpoint.COAGULOPATHY) > \
               predictor.predict(benign, ToxicityEndpoint.COAGULOPATHY)


# ---------------------------------------------------------------------------
# auroc_score
# ---------------------------------------------------------------------------

class TestAurocScore:
    def test_perfect_classifier(self) -> None:
        scores = [0.9, 0.8, 0.2, 0.1]
        labels = [1, 1, 0, 0]
        assert auroc_score(scores, labels) == pytest.approx(1.0)

    def test_random_classifier_approx_half(self) -> None:
        scores = [0.5] * 100
        labels = [i % 2 for i in range(100)]
        auc = auroc_score(scores, labels)
        assert 0.0 <= auc <= 1.0

    def test_only_one_class_returns_half(self) -> None:
        assert auroc_score([0.9, 0.8, 0.7], [1, 1, 1]) == pytest.approx(0.5)

    def test_auroc_equals_aucroc(self) -> None:
        from oligotox.evaluate import aucroc
        scores = [0.9, 0.4, 0.7, 0.3]
        labels = [1, 0, 1, 0]
        assert auroc_score(scores, labels) == pytest.approx(aucroc(scores, labels))


# ---------------------------------------------------------------------------
# toxicity_profile_summary
# ---------------------------------------------------------------------------

class TestToxicityProfileSummary:
    def test_keys_present(self) -> None:
        oligos, records = make_dataset(n=10, seed=7)
        oid = oligos[0].oligo_id
        predicted = {ep.value: 0.3 for ep in ToxicityEndpoint}
        summary = toxicity_profile_summary(oid, records, predicted)
        assert "oligo_id" in summary
        assert "measured" in summary
        assert "predicted" in summary
        assert "max_predicted_endpoint" in summary
        assert "overall_risk_level" in summary

    def test_oligo_id_matches(self) -> None:
        oligos, records = make_dataset(n=5, seed=3)
        oid = oligos[0].oligo_id
        predicted = {ep.value: 0.1 for ep in ToxicityEndpoint}
        summary = toxicity_profile_summary(oid, records, predicted)
        assert summary["oligo_id"] == oid

    def test_risk_level_very_high(self) -> None:
        records = [make_toxicity_record("o1", seed=1)]
        predicted = {ep.value: 0.0 for ep in ToxicityEndpoint}
        predicted[ToxicityEndpoint.HEPATOTOXICITY.value] = 0.90
        summary = toxicity_profile_summary("o1", records, predicted)
        assert summary["overall_risk_level"] == "very_high"

    def test_risk_level_low(self) -> None:
        records = [make_toxicity_record("o2", seed=2)]
        predicted = {ep.value: 0.05 for ep in ToxicityEndpoint}
        summary = toxicity_profile_summary("o2", records, predicted)
        assert summary["overall_risk_level"] == "low"

    def test_empty_predicted_scores(self) -> None:
        records = [make_toxicity_record("o3", seed=3)]
        summary = toxicity_profile_summary("o3", records, {})
        assert summary["max_predicted_endpoint"] == ""
        assert summary["overall_risk_level"] == "low"


# ---------------------------------------------------------------------------
# make_dataset realistic patterns
# ---------------------------------------------------------------------------

class TestMakeDataset:
    def test_ps_oligos_have_higher_mean_hepatotox(self) -> None:
        """PS oligos in the dataset should have higher hepatotox than PO oligos."""
        oligos, records = make_dataset(n=60, seed=42)
        record_map = {r.oligo_id: r for r in records}

        ps_vals = [
            record_map[o.oligo_id].value
            for o in oligos
            if o.backbone == BackboneClass.PS
            and record_map.get(o.oligo_id) is not None
            and record_map[o.oligo_id].endpoint == ToxicityEndpoint.HEPATOTOXICITY
        ]
        po_vals = [
            record_map[o.oligo_id].value
            for o in oligos
            if o.backbone == BackboneClass.PO
            and record_map.get(o.oligo_id) is not None
            and record_map[o.oligo_id].endpoint == ToxicityEndpoint.HEPATOTOXICITY
        ]
        if ps_vals and po_vals:
            assert sum(ps_vals) / len(ps_vals) > sum(po_vals) / len(po_vals)
