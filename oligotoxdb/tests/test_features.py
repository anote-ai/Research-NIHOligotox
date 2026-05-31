"""Tests for oligonucleotide feature computation."""

import pytest
import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from oligotoxdb.features import (
    compute_features, compute_features_batch,
    _max_run, _reverse_complement, _nearest_neighbor_tm, _estimate_mfe,
)


class TestComputeFeatures:
    def test_basic_sequence(self):
        feat = compute_features("test-01", "ATCGATCGATCG")
        assert feat.length == 12
        assert abs(feat.gc_content - 0.5) < 0.01
        assert feat.oligo_id == "test-01"

    def test_gc_content_pure_gc(self):
        feat = compute_features("gc-01", "GCGCGCGCGC")
        assert abs(feat.gc_content - 1.0) < 0.01

    def test_gc_content_pure_at(self):
        feat = compute_features("at-01", "ATATATATATAT")
        assert feat.gc_content == 0.0

    def test_cpg_count(self):
        feat = compute_features("cpg-01", "ACGACGACG")
        assert feat.cpg_count == 3  # three CG dinucleotides

    def test_cpg_count_none(self):
        feat = compute_features("nocpg-01", "AAAAAAAAAAAA")
        assert feat.cpg_count == 0

    def test_g4_score_high(self):
        # 4 G-tracts of length 4 = classic G4 motif
        feat = compute_features("g4-01", "GGGGTTTTGGGGTTTTGGGGTTTTGGGG")
        assert feat.g4_score >= 4

    def test_g4_score_zero(self):
        feat = compute_features("nog4-01", "ATCGATCGATCGATCG")
        assert feat.g4_score == 0

    def test_poly_g_run(self):
        feat = compute_features("polyg-01", "ATGGGGGAT")
        assert feat.max_poly_g == 5

    def test_backbone_one_hot_ps(self):
        feat = compute_features("ps-01", "ATCGATCG", backbone="PS")
        assert feat.backbone_ps == 1
        assert feat.backbone_po == 0

    def test_backbone_one_hot_pmo(self):
        feat = compute_features("pmo-01", "ATCGATCG", backbone="PMO")
        assert feat.backbone_pmo == 1
        assert feat.backbone_ps == 0

    def test_conjugate_galnac(self):
        feat = compute_features("gn-01", "ATCGATCG", conjugate="GalNAc")
        assert feat.conjugate_galnac == 1
        assert feat.conjugate_none == 0

    def test_rna_u_to_t_conversion(self):
        feat = compute_features("rna-01", "AUCGAUCGAUCG")
        assert "U" not in feat.sequence
        assert feat.length == 12

    def test_gapmer_encoding(self):
        feat = compute_features("gap-01", "ATCGATCGATCGATCGATCG",
                                gapmer_gap=10, gapmer_wing_5p=5, gapmer_wing_3p=5)
        assert feat.is_gapmer == 1
        assert feat.gap_length == 10
        assert feat.wing_length_5p == 5

    def test_feature_vector_length_consistent(self):
        feat1 = compute_features("t1", "ATCGATCG", backbone="PS")
        feat2 = compute_features("t2", "GCGCGCGC", backbone="PMO")
        assert len(feat1.to_vector()) == len(feat2.to_vector())

    def test_molecular_weight_positive(self):
        feat = compute_features("mw-01", "ATCGATCGATCG")
        assert feat.molecular_weight > 0

    def test_tm_reasonable_range(self):
        # Pure GC 20-mer will have high Tm; AT-rich will have low Tm
        feat_gc = compute_features("tm-gc", "GCGCGCGCGCGCGCGCGCGC")
        feat_at = compute_features("tm-at", "ATATATATATATATATATATAT")
        assert feat_gc.tm_estimate > 20
        assert feat_at.tm_estimate < feat_gc.tm_estimate

    def test_batch_computation(self):
        df = pd.DataFrame([
            {"oligo_id": "b-01", "sequence": "ATCGATCG", "backbone": "PS",
             "sugar_mod": "DNA", "conjugate": "none",
             "gapmer_gap": 0, "gapmer_wing_5p": 0, "gapmer_wing_3p": 0},
            {"oligo_id": "b-02", "sequence": "GCGCGCGC", "backbone": "PMO",
             "sugar_mod": "2OMe", "conjugate": "GalNAc",
             "gapmer_gap": 0, "gapmer_wing_5p": 0, "gapmer_wing_3p": 0},
        ])
        result = compute_features_batch(df)
        assert len(result) == 2
        assert "gc_content" in result.columns
        assert result.iloc[1]["backbone_pmo"] == 1


class TestHelpers:
    def test_max_run_basic(self):
        assert _max_run("ATGGGCAT", "G") == 3
        assert _max_run("AAAAAAA", "A") == 7
        assert _max_run("GCGCGC", "A") == 0

    def test_reverse_complement(self):
        assert _reverse_complement("ATCG") == "CGAT"
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("GCGC") == "GCGC"

    def test_mfe_negative(self):
        mfe = _estimate_mfe("GCGCGCGCGCGCGCGC")
        assert mfe < 0

    def test_tm_gc_dependence(self):
        tm_low_gc = _nearest_neighbor_tm("ATATATATAT")
        tm_high_gc = _nearest_neighbor_tm("GCGCGCGCGC")
        assert tm_high_gc > tm_low_gc
