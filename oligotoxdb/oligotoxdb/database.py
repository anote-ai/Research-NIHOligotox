"""
OligoToxDB DuckDB interface.

DuckDB is chosen because:
  - Columnar storage → fast aggregations over large result tables
  - Single-file deployment → easy sharing of the full database
  - Pandas/Parquet interop → seamless with the ML pipeline
  - No server required → runs in-process
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

DEFAULT_DB = Path(os.getenv("OLIGOTOXDB_PATH", "oligotoxdb.duckdb"))

SCHEMA_SQL = """
-- Compound registry
CREATE TABLE IF NOT EXISTS compound (
    oligo_id        VARCHAR PRIMARY KEY,
    sequence        VARCHAR NOT NULL,
    length          INTEGER,
    backbone_class  VARCHAR,
    sugar_mod       VARCHAR,
    gapmer_design   VARCHAR,
    conjugate       VARCHAR,
    synthesis_vendor VARCHAR,
    purity_percent  DOUBLE,
    mass_confirmed  BOOLEAN,
    endotoxin_eu_ml DOUBLE,
    batch_group     VARCHAR,
    created_at      TIMESTAMP DEFAULT current_timestamp
);

-- Pre-computed features
CREATE TABLE IF NOT EXISTS compound_features (
    oligo_id        VARCHAR PRIMARY KEY REFERENCES compound(oligo_id),
    gc_content      DOUBLE,
    cpg_count       INTEGER,
    mfe_estimate    DOUBLE,
    tm_estimate     DOUBLE,
    g4_score        DOUBLE,
    self_comp_index DOUBLE,
    molecular_weight DOUBLE,
    net_charge_ph7  INTEGER,
    hydrophobicity_index DOUBLE,
    feature_json    VARCHAR  -- full JSON blob of all features for ML ingestion
);

-- Raw experimental results
CREATE TABLE IF NOT EXISTS result (
    result_id       BIGINT PRIMARY KEY,
    oligo_id        VARCHAR REFERENCES compound(oligo_id),
    assay_system    VARCHAR NOT NULL,  -- PHH | KidneyOrganoid | PBMC | Platelet | MPS
    endpoint        VARCHAR NOT NULL,
    concentration_um DOUBLE,
    timepoint_h     INTEGER,
    replicate_bio   INTEGER,
    replicate_tech  INTEGER,
    batch_id        VARCHAR,
    plate_id        VARCHAR,
    cell_donor_id   VARCHAR,
    value_raw       DOUBLE,
    value_normalized DOUBLE,   -- % of vehicle control
    unit            VARCHAR,
    plate_qc_pass   BOOLEAN,
    qc_flag         VARCHAR DEFAULT 'pass'
);

-- Dose-response summary
CREATE TABLE IF NOT EXISTS dose_response (
    dr_id           BIGINT PRIMARY KEY,
    oligo_id        VARCHAR REFERENCES compound(oligo_id),
    assay_system    VARCHAR,
    endpoint        VARCHAR,
    ic50_um         DOUBLE,
    ic50_ci_low     DOUBLE,
    ic50_ci_high    DOUBLE,
    emax_percent    DOUBLE,
    hill_coefficient DOUBLE,
    bottom          DOUBLE,
    noec_um         DOUBLE,
    audrc           DOUBLE,
    r_squared       DOUBLE,
    n_bio_reps      INTEGER,
    fit_success     BOOLEAN,
    toxicity_class  VARCHAR,   -- inactive | low | moderate | high
    qc_flag         VARCHAR
);

-- Plate QC log
CREATE TABLE IF NOT EXISTS plate_qc (
    plate_id        VARCHAR PRIMARY KEY,
    batch_id        VARCHAR,
    z_prime         DOUBLE,
    pos_ctrl_mean   DOUBLE,
    pos_ctrl_cv     DOUBLE,
    neg_ctrl_mean   DOUBLE,
    neg_ctrl_cv     DOUBLE,
    signal_window   DOUBLE,
    passed          BOOLEAN,
    failure_reasons VARCHAR
);

-- Transcriptomics (RNA-seq gene-level summary; full data in GEO)
CREATE TABLE IF NOT EXISTS transcriptomics (
    oligo_id        VARCHAR REFERENCES compound(oligo_id),
    assay_system    VARCHAR,
    concentration_um DOUBLE,
    gene_symbol     VARCHAR,
    log2fc          DOUBLE,
    padj            DOUBLE,
    tpm             DOUBLE,
    geo_accession   VARCHAR
);
"""

IC50_THRESHOLDS = {
    "inactive": 50.0,
    "low": 10.0,
    "moderate": 1.0,
    # < 1 µM → high
}


class OligoToxDB:
    """
    Thin wrapper around a DuckDB connection for OligoToxDB.
    Usage:
        db = OligoToxDB()
        db.insert_compounds(compound_df)
        ic50s = db.query("SELECT * FROM dose_response WHERE assay_system='PHH'")
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
        """Insert rows from a compound DataFrame. Returns rows inserted."""
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

    # ── Queries ──────────────────────────────────────────────────────────

    def query(self, sql: str) -> pd.DataFrame:
        return self.conn.execute(sql).df()

    def get_compound(self, oligo_id: str) -> Optional[pd.DataFrame]:
        df = self.query(f"SELECT * FROM compound WHERE oligo_id = '{oligo_id}'")
        return df if len(df) > 0 else None

    def get_toxicity_profile(self, oligo_id: str) -> pd.DataFrame:
        """All dose-response summaries for one oligo."""
        return self.query(
            f"SELECT * FROM dose_response WHERE oligo_id = '{oligo_id}' ORDER BY assay_system, endpoint"
        )

    def get_ml_dataset(
        self,
        assay_system: Optional[str] = None,
        endpoint: Optional[str] = None,
        qc_flag: str = "pass",
    ) -> pd.DataFrame:
        """
        Join compound_features with dose_response for ML model training.
        Returns one row per (oligo, assay_system, endpoint) with features + label.
        """
        where_clauses = [f"dr.qc_flag = '{qc_flag}'"]
        if assay_system:
            where_clauses.append(f"dr.assay_system = '{assay_system}'")
        if endpoint:
            where_clauses.append(f"dr.endpoint = '{endpoint}'")
        where = " AND ".join(where_clauses)

        return self.query(f"""
            SELECT cf.*, dr.assay_system, dr.endpoint,
                   dr.ic50_um, dr.emax_percent, dr.audrc,
                   dr.toxicity_class, dr.r_squared
            FROM compound_features cf
            JOIN dose_response dr ON cf.oligo_id = dr.oligo_id
            WHERE {where}
        """)

    def summary_stats(self) -> dict:
        """Quick database summary."""
        return {
            "n_compounds": self.query("SELECT COUNT(*) AS n FROM compound").iloc[0, 0],
            "n_results": self.query("SELECT COUNT(*) AS n FROM result").iloc[0, 0],
            "n_dr_fits": self.query("SELECT COUNT(*) AS n FROM dose_response").iloc[0, 0],
            "n_plates": self.query("SELECT COUNT(*) AS n FROM plate_qc").iloc[0, 0],
            "assay_systems": self.query("SELECT DISTINCT assay_system FROM dose_response")["assay_system"].tolist(),
            "endpoints": self.query("SELECT COUNT(DISTINCT endpoint) AS n FROM dose_response").iloc[0, 0],
        }

    def export_parquet(self, output_dir: Path) -> None:
        """Export all tables as Parquet files for public data release."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for table in ["compound", "compound_features", "dose_response", "plate_qc", "transcriptomics"]:
            df = self.query(f"SELECT * FROM {table}")
            out_path = output_dir / f"{table}.parquet"
            df.to_parquet(out_path, index=False)
            print(f"Exported {len(df)} rows → {out_path}")


def _classify_toxicity(ic50: Optional[float]) -> str:
    if ic50 is None:
        return "inactive"
    if ic50 >= IC50_THRESHOLDS["inactive"]:
        return "inactive"
    if ic50 >= IC50_THRESHOLDS["low"]:
        return "low"
    if ic50 >= IC50_THRESHOLDS["moderate"]:
        return "moderate"
    return "high"
