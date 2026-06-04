"""Tests for oligotox.core."""
import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oligotox.core import (
    BackboneClass, ToxicityEndpoint, Oligonucleotide, ToxicityRecord,
    compute_gc_content, compute_cpg_ratio, extract_features, FeatureVector,
)


def test_backbone_class_values():
    values = {b.value for b in BackboneClass}
    assert values == {"PS", "PO", "PMO", "LNA", "PNA", "F2_ANA", "MORPHOLINO"}


def test_backbone_class_count():
    assert len(BackboneClass) == 7


def test_toxicity_endpoint_values():
    values = {e.value for e in ToxicityEndpoint}
    assert values == {"HEPATOTOXICITY", "NEPHROTOXICITY", "IMMUNOTOXICITY", "COMPLEMENT", "COAGULOPATHY"}


def test_toxicity_endpoint_count():
    assert len(ToxicityEndpoint) == 5


def test_oligonucleotide_auto_gc():
    oligo = Oligonucleotide(oligo_id="o1", sequence="ATGCGC", backbone=BackboneClass.PS)
    assert abs(oligo.gc_content - 4/6) < 1e-9


def test_compute_gc_content_known():
    assert abs(compute_gc_content("ATGCGC") - 4/6) < 1e-9


def test_compute_gc_content_empty_raises():
    with pytest.raises(ValueError):
        compute_gc_content("")


def test_compute_cpg_ratio_no_cpg():
    assert compute_cpg_ratio("ATGC") == 0.0


def test_compute_cpg_ratio_with_cpg():
    assert compute_cpg_ratio("CGCG") > 0


def test_extract_features_keys():
    oligo = Oligonucleotide(oligo_id="o1", sequence="ATGCATGCAT", backbone=BackboneClass.LNA)
    fv = extract_features(oligo)
    assert isinstance(fv, FeatureVector)
    assert len(fv.features) == 6
    assert fv.oligo_id == "o1"


def test_toxicity_record_confidence_in_range():
    rec = ToxicityRecord(oligo_id="o1", endpoint=ToxicityEndpoint.HEPATOTOXICITY,
                         value=0.5, cell_system="PHH", confidence=0.8)
    assert 0.0 <= rec.confidence <= 1.0
