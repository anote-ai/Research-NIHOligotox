"""Tests for the OTMRS validator."""

import sys
import io
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from oligotoxdb.otmrs_validator import (
    validate_compounds, validate_results, validate_submission, OTMRSReport,
)


def make_valid_compound_df(n: int = 3) -> pd.DataFrame:
    return pd.DataFrame([{
        "oligo_id": f"OT-{i+1:05d}",
        "sequence": "ATCGATCGATCG",
        "backbone_class": "PS",
        "sugar_mod": "DNA",
        "conjugate": "none",
        "synthesis_vendor": "IDT",
        "purity_percent": 98.0,
        "mass_confirmed": True,
        "endotoxin_eu_ml": 0.1,
    } for i in range(n)])


def make_valid_results_df(oligo_ids: list[str]) -> pd.DataFrame:
    rows = []
    for oid in oligo_ids:
        for conc in [0.1, 1.0, 10.0, 50.0]:
            for rep in [1, 2]:
                rows.append({
                    "oligo_id": oid,
                    "assay_system": "PHH",
                    "endpoint": "Cell_viability_ATPLite",
                    "concentration_um": conc,
                    "timepoint_h": 24,
                    "replicate_bio": rep,
                    "replicate_tech": 1,
                    "batch_id": "BATCH-A",
                    "plate_id": "PLT-001",
                    "cell_donor_id": "DONOR-01",
                    "value_raw": 95.0,
                    "unit": "% vehicle",
                })
    return pd.DataFrame(rows)


class TestCompoundValidation:
    def test_valid_compounds_pass(self):
        df = make_valid_compound_df()
        errors = validate_compounds(df)
        assert len([e for e in errors if e.severity == "error"]) == 0

    def test_invalid_sequence_flagged(self):
        df = make_valid_compound_df(1)
        df.loc[0, "sequence"] = "ATCGXYZ"  # X, Y, Z are invalid
        errors = validate_compounds(df)
        assert any(e.field == "sequence" and e.severity == "error" for e in errors)

    def test_invalid_backbone_flagged(self):
        df = make_valid_compound_df(1)
        df.loc[0, "backbone_class"] = "FAKE_BACKBONE"
        errors = validate_compounds(df)
        assert any(e.field == "backbone_class" and e.severity == "error" for e in errors)

    def test_low_purity_flagged_as_warning(self):
        df = make_valid_compound_df(1)
        df.loc[0, "purity_percent"] = 82.0  # below MIN_PURITY=85 but above 80 → warning
        errors = validate_compounds(df)
        purity_errs = [e for e in errors if e.field == "purity_percent"]
        assert len(purity_errs) > 0
        assert purity_errs[0].severity == "warning"

    def test_very_low_purity_flagged_as_error(self):
        df = make_valid_compound_df(1)
        df.loc[0, "purity_percent"] = 75.0  # below 80%
        errors = validate_compounds(df)
        purity_errs = [e for e in errors if e.field == "purity_percent" and e.severity == "error"]
        assert len(purity_errs) > 0

    def test_duplicate_oligo_id_flagged(self):
        df = make_valid_compound_df(2)
        df.loc[1, "oligo_id"] = df.loc[0, "oligo_id"]  # duplicate
        errors = validate_compounds(df)
        assert any(e.field == "oligo_id" and e.severity == "error" for e in errors)

    def test_missing_required_column(self):
        df = make_valid_compound_df().drop(columns=["sequence"])
        errors = validate_compounds(df)
        assert any(e.field == "sequence" and e.severity == "error" for e in errors)


class TestResultValidation:
    def test_valid_results_pass(self):
        compounds = make_valid_compound_df(2)
        results = make_valid_results_df(["OT-00001", "OT-00002"])
        errors = validate_results(results, compound_ids=set(compounds["oligo_id"]))
        assert len([e for e in errors if e.severity == "error"]) == 0

    def test_unknown_oligo_id_flagged(self):
        results = make_valid_results_df(["OT-00001"])
        errors = validate_results(results, compound_ids={"OT-00002"})  # OT-00001 not in registry
        assert any(e.field == "oligo_id" and e.severity == "error" for e in errors)

    def test_insufficient_replicates_warned(self):
        results = make_valid_results_df(["OT-00001"])
        # Remove all replicate 2 rows
        results = results[results["replicate_bio"] == 1]
        errors = validate_results(results)
        assert any("replicate" in e.field.lower() for e in errors)

    def test_invalid_assay_system_flagged(self):
        results = make_valid_results_df(["OT-00001"])
        results.loc[0, "assay_system"] = "FAKE_SYSTEM"
        errors = validate_results(results)
        assert any(e.field == "assay_system" and e.severity == "error" for e in errors)

    def test_negative_concentration_flagged(self):
        results = make_valid_results_df(["OT-00001"])
        results.loc[0, "concentration_um"] = -1.0
        errors = validate_results(results)
        assert any(e.field == "concentration_um" and e.severity == "error" for e in errors)
