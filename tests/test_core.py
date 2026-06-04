"""Tests for oligotox.core module."""

import pytest
from oligotox.core import (
    BackboneClass,
    ToxicityEndpoint,
    Oligonucleotide,
    ToxicityRecord,
    FeatureVector,
    compute_gc_content,
    compute_cpg_ratio,
    extract_features,
)


def test_backbone_class_enum_values() -> None:
    assert BackboneClass.PS == "PS"
    assert BackboneClass.PO == "PO"
    assert BackboneClass.LNA == "LNA"
    assert len(BackboneClass) == 7


def test_toxicity_endpoint_enum_values() -> None:
    assert ToxicityEndpoint.HEPATOTOXICITY == "HEPATOTOXICITY"
    assert ToxicityEndpoint.COMPLEMENT == "COMPLEMENT"
    assert len(ToxicityEndpoint) == 5


def test_oligonucleotide_construction() -> None:
    oligo = Oligonucleotide(
        oligo_id="OLG001",
        sequence="ATGCGC",
        backbone=BackboneClass.PS,
        modifications=["2'-OMe"],
    )
    assert oligo.oligo_id == "OLG001"
    assert oligo.backbone == BackboneClass.PS
    assert oligo.modifications == ["2'-OMe"]
    # gc_content auto-computed
    assert abs(oligo.gc_content - 0.5) < 1e-9


def test_toxicity_record_construction() -> None:
    record = ToxicityRecord(
        oligo_id="OLG001",
        endpoint=ToxicityEndpoint.HEPATOTOXICITY,
        value=0.72,
        cell_system="HepG2",
        confidence=0.9,
    )
    assert record.endpoint == ToxicityEndpoint.HEPATOTOXICITY
    assert record.value == 0.72


def test_toxicity_record_invalid_confidence() -> None:
    with pytest.raises(ValueError):
        ToxicityRecord(
            oligo_id="OLG002",
            endpoint=ToxicityEndpoint.NEPHROTOXICITY,
            value=0.5,
            cell_system="HEK293",
            confidence=1.5,
        )


def test_compute_gc_content_known_sequence() -> None:
    assert compute_gc_content("ATGCGC") == pytest.approx(0.5, abs=1e-9)


def test_compute_gc_content_all_at() -> None:
    assert compute_gc_content("AAATTT") == pytest.approx(0.0)


def test_compute_gc_content_empty() -> None:
    assert compute_gc_content("") == 0.0


def test_compute_cpg_ratio_basic() -> None:
    # "ACGT" has one CpG out of 3 dinucleotides
    ratio = compute_cpg_ratio("ACGT")
    assert ratio == pytest.approx(1 / 3, abs=1e-9)


def test_compute_cpg_ratio_no_cpg() -> None:
    assert compute_cpg_ratio("AAAA") == pytest.approx(0.0)


def test_compute_cpg_ratio_short_sequence() -> None:
    assert compute_cpg_ratio("A") == 0.0


def test_extract_features_returns_feature_vector() -> None:
    oligo = Oligonucleotide(
        oligo_id="OLG003",
        sequence="ATCGCGAT",
        backbone=BackboneClass.PS,
    )
    fv = extract_features(oligo)
    assert isinstance(fv, FeatureVector)
    assert fv.oligo_id == "OLG003"
    assert "gc_content" in fv.features
    assert "cpg_ratio" in fv.features
    assert "sequence_length" in fv.features
    assert "has_ps" in fv.features


def test_extract_features_has_ps_flag() -> None:
    oligo_ps = Oligonucleotide(oligo_id="PS1", sequence="ATCG", backbone=BackboneClass.PS)
    oligo_po = Oligonucleotide(oligo_id="PO1", sequence="ATCG", backbone=BackboneClass.PO)
    assert extract_features(oligo_ps).features["has_ps"] == 1.0
    assert extract_features(oligo_po).features["has_ps"] == 0.0
