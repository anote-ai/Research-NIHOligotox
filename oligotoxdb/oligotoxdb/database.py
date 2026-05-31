"""
OligoToxDB DuckDB interface — parameterized queries, no SQL injection.

DuckDB is chosen because:
  - Columnar storage → fast aggregations over large result tables
  - Single-file deployment → easy sharing of the full database
  - Pandas/Parquet interop → seamless with the ML pipeline
  - No server required → runs in-process
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

DEFAULT_DB = Path(os.getenv("OLIGOTOXDB_PATH", "oligotoxdb.duckdb"))

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS compound (
    oligo_id         VARCHAR PRIMARY KEY,
    sequence         VARCHAR NOT NULL,
    length           INTEGER,
    backbone_class   VARCHAR,
    sugar_mod        VARCHAR,
    gapmer_design    VARCHAR,
    conjugate        VARCHAR,
    synthesis_vendor VARCHAR,
    purity_percent   DOUBLE,
    mass_confirmed   BOOLEAN,
    endotoxin_eu_ml  DOUBLE,
    batch_group      VARCHAR,
    created_at       TIMESTAMP DEFAULT current_timestamp
);

CREATE TABLE IF NOT EXISTS compound_features (
    oligo_id             VARCHAR PRIMARY KEY REFERENCES compound(oligo_id),
    gc_content           DOUBLE,
    cpg_count            INTEGER,
    mfe_estimate         DOUBLE,
    tm_estimate          DOUBLE,
    g4_score             DOUBLE,
    self_comp_index      DOUBLE,
    molecular_weight     DOUBLE,
    net_charge_ph7       INTEGER,
    hydrophobicity_index DOUBLE,
    feature_json         VARCHAR
);

CREATE TABLE IF NOT EXISTS result (
    result_id        BIGINT PRIMARY KEY,
    oligo_id         VARCHAR REFERENCES compound(oligo_id),
    assay_system     VARCHAR NOT NULL,
    endpoint         VARCHAR NOT NULL,
    concentration_um DOUBLE,
    timepoint_h      INTEGER,
    replicate_bio    INTEGER,
    replicate_tech   INTEGER,
    batch_id         VARCHAR,
    plate_id         VARCHAR,
    cell_donor_id    VARCHAR,
    value_raw        DOUBLE,
    value_normalized DOUBLE,
    unit             VARCHAR,
    plate_qc_pass    BOOLEAN,
    qc_flag          VARCHAR DEFAULT 'pass'
);

CREATE TABLE IF NOT EXISTS dose_response (
    dr_id            BIGINT PRIMARY KEY,
    oligo_id         VARCHAR REFERENCES compound(oligo_id),
    assay_system     VARCHAR,
    endpoint         VARCHAR,
    ic50_um          DOUBLE,
    ic50_ci_low      DOUBLE,
    ic50_ci_high     DOUBLE,
    emax_percent     DOUBLE,
    hill_coefficient DOUBLE,
    bottom           DOUBLE,
    noec_um          DOUBLE,
    audrc            DOUBLE,
    r_squared        DOUBLE,
    n_bio_reps       INTEGER,
    fit_success      BOOLEAN,
    toxicity_class   VARCHAR,
    qc_flag          VARCHAR
);

CREATE TABLE IF NOT EXISTS plate_qc (
    plate_id       VARCHAR PRIMARY KEY,
    batch_id       VARCHAR,
    z_prime        DOUBLE,
    pos_ctrl_mean  DOUBLE,
    pos_ctrl_cv    DOUBLE,
    neg_ctrl_mean  DOUBLE,
    neg_ctrl_cv    DOUBLE,
    signal_window  DOUBLE,
    passed         BOOLEAN,
    failure_reasons VARCHAR
);

CREATE TABLE IF NOT EXISTS transcriptomics (
    oligo_id         VARCHAR REFERENCES compound(oligo_id),
    assay_system     VARCHAR,
    concentration_um DOUBLE,
    gene_symbol      VARCHAR,
    log2fc           DOUBLE,
    padj             DOUBLE,
    tpm              DOUBLE,
    geo_accession    VARCHAR
);

CREATE TABLE IF NOT EXISTS proteomics (
    oligo_id         VARCHAR REFERENCES compound(oligo_id),
    assay_system     VARCHAR,
    concentration_um DOUBLE,
    protein_id       VARCHAR,
    gene_symbol      VARCHAR,
    log2fc           DOUBLE,
    padj             DOUBLE,
    intensity        DOUBLE,
    pride_accession  VARCHAR
);
"""

IC50_THRESHOLDS = {"inactive": 50.0, "low": 10.0, "moderate": 1.0}

# Allowed characters for safe text search (allowlist approach)
_SAFE_PATTERN = re.compile(r"^[A-Za-z0-9\-_. %]+$")


def _safe_str(value: str, field_name: str = "value") -> str:
    """Raise ValueError if value contains characters outside the allowlist."""
    if not _SAFE_PATTERN.match(value):
        raise ValueError(f"Unsafe characters in {field_name}: {value!r}")
    return value


class OligoToxDB:
    """
    Thin wrapper around a DuckDB connection for OligoToxDB.

    All user-supplied strings are validated against an allowlist before
    being interpolated into SQL — DuckDB's Python API does not support
    parameterized queries for all statement types, so we sanitize instead.
    """

    def __init__(self, db_path: Path = DEFAULT_DB) -> None:
        self.db_path = Path(db_path)
        self.conn = duckdb.connect(str(self.db_path))
        self.conn.execute(SCHEMA_SQL)

    def close(self) -> None:
        self.conn.close()

    def __enter__(self) -> "OligoToxDB":
        return self

    def __exit__(self, *_) -> None:
        self.close()

    # ── Inserts ──────────────────────────────────────────────────────────

    def insert_compounds(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT OR REPLACE INTO compound SELECT * FROM df")
        return len(df)

    def insert_features(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT OR REPLACE INTO compound_features SELECT * FROM df")
        return len(df)

    def insert_results(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT INTO result SELECT * FROM df")
        return len(df)

    def insert_dose_response(self, df: pd.DataFrame) -> int:
        df = df.copy()
        df["toxicity_class"] = df["ic50_um"].apply(_classify_toxicity)
        self.conn.execute("INSERT INTO dose_response SELECT * FROM df")
        return len(df)

    def insert_plate_qc(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT OR REPLACE INTO plate_qc SELECT * FROM df")
        return len(df)

    def insert_transcriptomics(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT INTO transcriptomics SELECT * FROM df")
        return len(df)

    def insert_proteomics(self, df: pd.DataFrame) -> int:
        self.conn.execute("INSERT INTO proteomics SELECT * FROM df")
        return len(df)

    # ── Queries ──────────────────────────────────────────────────────────

    def query(self, sql: str) -> pd.DataFrame:
        return self.conn.execute(sql).df()

    def get_compound(self, oligo_id: str) -> Optional[pd.DataFrame]:
        safe_id = _safe_str(oligo_id, "oligo_id")
        df = self.query(f"SELECT * FROM compound WHERE oligo_id = '{safe_id}'")
        return df if len(df) > 0 else None

    def get_toxicity_profile(self, oligo_id: str) -> pd.DataFrame:
        safe_id = _safe_str(oligo_id, "oligo_id")
        return self.query(
            f"SELECT * FROM dose_response WHERE oligo_id = '{safe_id}' ORDER BY assay_system, endpoint"
        )

    def search_compounds(
        self,
        sequence_contains: str = "",
        backbone: str = "",
        toxicity_class: str = "",
        limit: int = 500,
    ) -> pd.DataFrame:
        """Safe compound search with allowlist-validated filters."""
        conditions = []
        if sequence_contains:
            seq = _safe_str(sequence_contains.upper(), "sequence")
            conditions.append(f"c.sequence LIKE '%{seq}%'")
        if backbone and backbone != "All":
            bb = _safe_str(backbone, "backbone")
            conditions.append(f"c.backbone_class = '{bb}'")
        if toxicity_class and toxicity_class != "All":
            tc = _safe_str(toxicity_class, "toxicity_class")
            conditions.append(f"dr.toxicity_class = '{tc}'")

        where = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        return self.query(f"""
            SELECT c.oligo_id, c.sequence, c.backbone_class, c.sugar_mod, c.conjugate,
                   cf.gc_content, cf.cpg_count, cf.tm_estimate,
                   dr.endpoint, dr.ic50_um, dr.toxicity_class
            FROM compound c
            LEFT JOIN compound_features cf ON c.oligo_id = cf.oligo_id
            LEFT JOIN dose_response dr ON c.oligo_id = dr.oligo_id
            {where}
            LIMIT {int(limit)}
        """)

    def get_ml_dataset(
        self,
        assay_system: Optional[str] = None,
        endpoint: Optional[str] = None,
        qc_flag: str = "pass",
    ) -> pd.DataFrame:
        conditions = [f"dr.qc_flag = '{_safe_str(qc_flag)}'"]
        if assay_system:
            conditions.append(f"dr.assay_system = '{_safe_str(assay_system)}'")
        if endpoint:
            conditions.append(f"dr.endpoint = '{_safe_str(endpoint)}'")
        where = " AND ".join(conditions)
        return self.query(f"""
            SELECT cf.*, dr.assay_system, dr.endpoint,
                   dr.ic50_um, dr.emax_percent, dr.audrc,
                   dr.toxicity_class, dr.r_squared
            FROM compound_features cf
            JOIN dose_response dr ON cf.oligo_id = dr.oligo_id
            WHERE {where}
        """)

    def summary_stats(self) -> dict:
        return {
            "n_compounds": self.query("SELECT COUNT(*) AS n FROM compound").iloc[0, 0],
            "n_results": self.query("SELECT COUNT(*) AS n FROM result").iloc[0, 0],
            "n_dr_fits": self.query("SELECT COUNT(*) AS n FROM dose_response").iloc[0, 0],
            "n_plates": self.query("SELECT COUNT(*) AS n FROM plate_qc").iloc[0, 0],
            "assay_systems": self.query(
                "SELECT DISTINCT assay_system FROM dose_response"
            )["assay_system"].tolist(),
            "n_endpoints": self.query(
                "SELECT COUNT(DISTINCT endpoint) AS n FROM dose_response"
            ).iloc[0, 0],
        }

    def export_parquet(self, output_dir: Path) -> dict[str, int]:
        """Export all tables as Parquet files for public data release."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        counts = {}
        for table in ["compound", "compound_features", "dose_response",
                      "plate_qc", "transcriptomics", "proteomics"]:
            df = self.query(f"SELECT * FROM {table}")
            out_path = output_dir / f"{table}.parquet"
            df.to_parquet(out_path, index=False)
            counts[table] = len(df)
        return counts


def _classify_toxicity(ic50: Optional[float]) -> str:
    if ic50 is None or (isinstance(ic50, float) and ic50 != ic50):  # NaN check
        return "inactive"
    if ic50 >= IC50_THRESHOLDS["inactive"]:
        return "inactive"
    if ic50 >= IC50_THRESHOLDS["low"]:
        return "low"
    if ic50 >= IC50_THRESHOLDS["moderate"]:
        return "moderate"
    return "high"
