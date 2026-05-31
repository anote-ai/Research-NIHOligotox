"""
ETL pipeline: ingest raw experimental data from CRO deliverables into OligoToxDB.

CROs deliver data in CSV format following the OligoTox Minimum Reporting Standard
(OTMRS). This module validates, normalizes, and loads those files into DuckDB.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import numpy as np
from pydantic import BaseModel, field_validator, ValidationError
import click
from rich.console import Console
from rich.table import Table

from .database import OligoToxDB
from .features import compute_features_batch
from .qc import run_qc_pipeline, normalize_batch_drift

console = Console()

# ─── OTMRS validation schema ─────────────────────────────────────────────────

VALID_BACKBONES = {"PS", "PO", "PMO", "LNA_mix", "PNA", "2F_ANA", "morpholino"}
VALID_SUGAR_MODS = {"DNA", "2OMe", "2F", "LNA", "mixed"}
VALID_CONJUGATES = {"none", "GalNAc", "cholesterol", "lipid", "antibody"}
VALID_ASSAY_SYSTEMS = {"PHH", "KidneyOrganoid", "PBMC", "Platelet", "MPS"}
VALID_QC_FLAGS = {"pass", "marginal", "fail"}

NUCLEOTIDES = set("ACGTU")


class CompoundRecord(BaseModel):
    oligo_id: str
    sequence: str
    backbone_class: str
    sugar_mod: str
    conjugate: str
    synthesis_vendor: str
    purity_percent: float
    mass_confirmed: bool
    endotoxin_eu_ml: float

    @field_validator("sequence")
    @classmethod
    def valid_sequence(cls, v: str) -> str:
        v = v.upper().strip()
        invalid = set(v) - NUCLEOTIDES
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        return v

    @field_validator("backbone_class")
    @classmethod
    def valid_backbone(cls, v: str) -> str:
        if v not in VALID_BACKBONES:
            raise ValueError(f"backbone_class must be one of {VALID_BACKBONES}")
        return v

    @field_validator("sugar_mod")
    @classmethod
    def valid_sugar(cls, v: str) -> str:
        if v not in VALID_SUGAR_MODS:
            raise ValueError(f"sugar_mod must be one of {VALID_SUGAR_MODS}")
        return v

    @field_validator("conjugate")
    @classmethod
    def valid_conjugate(cls, v: str) -> str:
        if v not in VALID_CONJUGATES:
            raise ValueError(f"conjugate must be one of {VALID_CONJUGATES}")
        return v

    @field_validator("purity_percent")
    @classmethod
    def valid_purity(cls, v: float) -> float:
        if not (0 < v <= 100):
            raise ValueError("purity_percent must be in (0, 100]")
        return v

    @field_validator("endotoxin_eu_ml")
    @classmethod
    def valid_endotoxin(cls, v: float) -> float:
        if v < 0:
            raise ValueError("endotoxin_eu_ml must be >= 0")
        return v


class ResultRecord(BaseModel):
    oligo_id: str
    assay_system: str
    endpoint: str
    concentration_um: float
    timepoint_h: int
    replicate_bio: int
    replicate_tech: int
    batch_id: str
    plate_id: str
    cell_donor_id: str
    value_raw: float
    unit: str

    @field_validator("assay_system")
    @classmethod
    def valid_assay(cls, v: str) -> str:
        if v not in VALID_ASSAY_SYSTEMS:
            raise ValueError(f"assay_system must be one of {VALID_ASSAY_SYSTEMS}")
        return v

    @field_validator("concentration_um")
    @classmethod
    def positive_conc(cls, v: float) -> float:
        if v < 0:
            raise ValueError("concentration_um must be >= 0")
        return v


# ─── Ingestion pipeline ───────────────────────────────────────────────────────

@dataclass
class IngestionReport:
    compound_rows_loaded: int
    result_rows_loaded: int
    validation_errors: list[dict]
    plate_qc_pass_rate: float
    dr_fit_success_rate: float


def ingest_compound_file(
    csv_path: Path,
    db: OligoToxDB,
    skip_invalid: bool = True,
) -> tuple[pd.DataFrame, list[dict]]:
    """Validate and load a compound registry CSV into OligoToxDB."""
    df = pd.read_csv(csv_path)
    errors = []
    valid_rows = []

    for idx, row in df.iterrows():
        try:
            rec = CompoundRecord(**row.to_dict())
            valid_rows.append(rec.model_dump())
        except ValidationError as e:
            errors.append({"row": int(idx), "oligo_id": row.get("oligo_id", "?"), "errors": str(e)})
            if not skip_invalid:
                raise

    valid_df = pd.DataFrame(valid_rows)
    if len(valid_df) > 0:
        # Add derived columns
        valid_df["length"] = valid_df["sequence"].str.len()
        valid_df["gapmer_design"] = valid_df.get("gapmer_design", "")
        valid_df["batch_group"] = valid_df.get("batch_group", "")
        db.insert_compounds(valid_df)

        # Compute and store features
        feat_df = compute_features_batch(valid_df)
        feat_df["feature_json"] = feat_df.apply(
            lambda r: json.dumps(r.drop(["oligo_id", "sequence"]).to_dict()), axis=1
        )
        db.insert_features(feat_df[["oligo_id", "gc_content", "cpg_count",
                                    "mfe_estimate", "tm_estimate", "g4_score",
                                    "self_comp_index", "molecular_weight",
                                    "net_charge_ph7", "hydrophobicity_index", "feature_json"]])

    return valid_df, errors


def ingest_results_file(
    results_csv: Path,
    controls_csv: Path,
    db: OligoToxDB,
    skip_invalid: bool = True,
) -> IngestionReport:
    """Validate raw results, run QC pipeline, and load everything into OligoToxDB."""
    raw_df = pd.read_csv(results_csv)
    ctrl_df = pd.read_csv(controls_csv)

    # Validate each result row
    errors = []
    valid_rows = []
    for idx, row in raw_df.iterrows():
        try:
            rec = ResultRecord(**row.to_dict())
            valid_rows.append(rec.model_dump())
        except ValidationError as e:
            errors.append({"row": int(idx), "errors": str(e)})
            if not skip_invalid:
                raise

    valid_df = pd.DataFrame(valid_rows) if valid_rows else pd.DataFrame()
    if valid_df.empty:
        return IngestionReport(0, 0, errors, 0.0, 0.0)

    # Normalize to % vehicle control
    vehicle_medians = (
        valid_df[valid_df["concentration_um"] == 0.0]
        .groupby(["plate_id", "endpoint"])["value_raw"]
        .median()
        .rename("vehicle_median")
    )
    valid_df = valid_df.merge(vehicle_medians, on=["plate_id", "endpoint"], how="left")
    valid_df["value_normalized"] = np.where(
        valid_df["vehicle_median"] > 0,
        (valid_df["value_raw"] / valid_df["vehicle_median"]) * 100,
        np.nan,
    )
    valid_df.drop(columns=["vehicle_median"], inplace=True)

    # Assign result IDs
    valid_df["result_id"] = range(1, len(valid_df) + 1)

    # QC pipeline
    plate_qc_df, dr_df = run_qc_pipeline(valid_df, ctrl_df)

    # Batch normalization using anchor compounds
    if "is_anchor" in valid_df.columns:
        valid_df = normalize_batch_drift(valid_df)

    # Set QC flags on result rows
    valid_df["qc_flag"] = valid_df["plate_qc_pass"].map({True: "pass", False: "fail"})

    # Load into DB
    db.insert_results(valid_df)
    db.insert_plate_qc(plate_qc_df)

    # Prepare dose_response table
    dr_df["dr_id"] = range(1, len(dr_df) + 1)
    dr_df = dr_df.rename(columns={"ic50": "ic50_um"})
    dr_df["assay_system"] = dr_df["endpoint"].apply(_infer_assay_system)
    dr_df["n_bio_reps"] = (
        valid_df.groupby(["oligo_id", "endpoint"])["replicate_bio"].nunique().reindex(
            dr_df.set_index(["oligo_id", "endpoint"]).index
        ).values
    )
    db.insert_dose_response(dr_df)

    pass_rate = float(plate_qc_df["passed"].mean()) if len(plate_qc_df) > 0 else 0.0
    dr_pass_rate = float(dr_df["fit_success"].mean()) if len(dr_df) > 0 else 0.0

    return IngestionReport(
        compound_rows_loaded=0,  # handled by ingest_compound_file
        result_rows_loaded=len(valid_df),
        validation_errors=errors,
        plate_qc_pass_rate=pass_rate,
        dr_fit_success_rate=dr_pass_rate,
    )


# ─── Helpers ──────────────────────────────────────────────────────────────────

_ENDPOINT_ASSAY_MAP = {
    "Cell_viability": "PHH",
    "LDH": "PHH",
    "ALT": "PHH",
    "AST": "PHH",
    "KIM1": "KidneyOrganoid",
    "NGAL": "KidneyOrganoid",
    "IFN": "PBMC",
    "IL": "PBMC",
    "TNF": "PBMC",
    "C3a": "PBMC",
    "C5a": "PBMC",
    "Platelet": "Platelet",
    "aPTT": "Platelet",
    "PT": "Platelet",
    "Chip": "MPS",
}

def _infer_assay_system(endpoint: str) -> str:
    for prefix, system in _ENDPOINT_ASSAY_MAP.items():
        if endpoint.startswith(prefix):
            return system
    return "Unknown"


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """OligoToxDB data ingestion pipeline."""


@cli.command()
@click.argument("compounds_csv", type=click.Path(exists=True))
@click.option("--db", default="oligotoxdb.duckdb", help="Path to DuckDB database")
def load_compounds(compounds_csv: str, db: str) -> None:
    """Load and validate a compound registry CSV into the database."""
    with OligoToxDB(db) as odb:
        valid_df, errors = ingest_compound_file(Path(compounds_csv), odb)
        console.print(f"[green]Loaded {len(valid_df)} compounds.[/green]")
        if errors:
            console.print(f"[yellow]{len(errors)} validation errors (skipped).[/yellow]")


@cli.command()
@click.argument("results_csv", type=click.Path(exists=True))
@click.argument("controls_csv", type=click.Path(exists=True))
@click.option("--db", default="oligotoxdb.duckdb")
def load_results(results_csv: str, controls_csv: str, db: str) -> None:
    """Run QC pipeline and load experimental results into the database."""
    with OligoToxDB(db) as odb:
        report = ingest_results_file(Path(results_csv), Path(controls_csv), odb)
        table = Table(title="Ingestion Report")
        table.add_column("Metric")
        table.add_column("Value")
        table.add_row("Results loaded", str(report.result_rows_loaded))
        table.add_row("Validation errors", str(len(report.validation_errors)))
        table.add_row("Plate QC pass rate", f"{report.plate_qc_pass_rate:.1%}")
        table.add_row("DR fit success rate", f"{report.dr_fit_success_rate:.1%}")
        console.print(table)


if __name__ == "__main__":
    cli()
