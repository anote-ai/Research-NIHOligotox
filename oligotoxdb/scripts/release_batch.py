"""
OligoToxDB rolling batch release script.

Packages a completed experimental batch into:
  1. FAIR-compliant Parquet + CSV data files
  2. Zenodo metadata JSON (ready for upload via Zenodo REST API)
  3. HuggingFace dataset card (dataset_card.md)
  4. SHA-256 checksums for all files
  5. Release changelog entry

Usage:
    python scripts/release_batch.py \\
        --batch-label A \\
        --db oligotoxdb.duckdb \\
        --output releases/batch_A/ \\
        --zenodo-token $ZENODO_TOKEN

Zenodo API docs: https://developers.zenodo.org/
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from datetime import date
from pathlib import Path
from typing import Optional

import pandas as pd
import click
from rich.console import Console
from rich.table import Table

sys.path.insert(0, str(Path(__file__).parent.parent))
from oligotoxdb.database import OligoToxDB

console = Console()

ZENODO_API = "https://zenodo.org/api"
ZENODO_SANDBOX_API = "https://sandbox.zenodo.org/api"

DATASET_METADATA = {
    "title": "OligoToxDB: Large-scale human in vitro oligonucleotide toxicity dataset",
    "description": (
        "OligoToxDB is the largest open dataset of human in vitro oligonucleotide toxicity profiles. "
        "Generated as part of the NIH NCATS OligoTox Open Data Challenge (Phase 2). "
        "Contains dose-response data for 2,800 oligonucleotides across 47 toxicity endpoints "
        "in 5 human in vitro model systems."
    ),
    "creators": [{"name": "Vidra, Natan", "affiliation": "Anote, Inc.", "orcid": ""}],
    "keywords": [
        "oligonucleotide", "toxicity", "antisense", "siRNA", "drug safety",
        "in vitro", "machine learning", "open data", "NIH NCATS",
    ],
    "license": {"id": "cc-by-4.0"},
    "access_right": "open",
    "upload_type": "dataset",
    "related_identifiers": [
        {"relation": "isPartOf", "identifier": "https://github.com/anote-ai/nih-oligotox"},
        {"relation": "isSupplementTo", "identifier": "https://github.com/anote-ai/oligotox-predict"},
    ],
    "communities": [{"identifier": "oligotoxdb"}],
    "grants": [{"id": "NIH-NCATS-OligoTox-Phase2"}],
}


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def export_batch(
    db: OligoToxDB,
    batch_label: str,
    output_dir: Path,
) -> dict[str, int]:
    """Export all data for a batch to Parquet and CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)

    counts = {}

    # Compounds in this batch
    compounds = db.query(
        f"SELECT * FROM compound WHERE batch_group = '{batch_label}'"
    )
    counts["compounds"] = len(compounds)
    compounds.to_parquet(output_dir / "compounds.parquet", index=False)
    compounds.to_csv(output_dir / "compounds.csv", index=False)

    if len(compounds) == 0:
        console.print(f"[yellow]No compounds found for batch '{batch_label}'[/yellow]")
        return counts

    oligo_ids = tuple(compounds["oligo_id"].tolist())
    id_clause = f"({','.join(repr(x) for x in oligo_ids)})"

    # Features
    features = db.query(f"SELECT * FROM compound_features WHERE oligo_id IN {id_clause}")
    counts["features"] = len(features)
    features.to_parquet(output_dir / "features.parquet", index=False)

    # Raw results
    results = db.query(f"SELECT * FROM result WHERE oligo_id IN {id_clause}")
    counts["results"] = len(results)
    results.to_parquet(output_dir / "results.parquet", index=False)
    results.to_csv(output_dir / "results.csv", index=False)

    # Dose-response summaries
    dr = db.query(f"SELECT * FROM dose_response WHERE oligo_id IN {id_clause}")
    counts["dose_response"] = len(dr)
    dr.to_parquet(output_dir / "dose_response.parquet", index=False)
    dr.to_csv(output_dir / "dose_response.csv", index=False)

    # Plate QC
    plates = db.query(
        f"SELECT * FROM plate_qc WHERE batch_id = 'BATCH-{batch_label}'"
    )
    counts["plate_qc"] = len(plates)
    plates.to_parquet(output_dir / "plate_qc.parquet", index=False)

    # Transcriptomics (if any)
    tx = db.query(f"SELECT * FROM transcriptomics WHERE oligo_id IN {id_clause}")
    if len(tx) > 0:
        counts["transcriptomics"] = len(tx)
        tx.to_parquet(output_dir / "transcriptomics.parquet", index=False)

    return counts


def write_checksums(output_dir: Path) -> Path:
    """Compute SHA-256 checksums for all output files."""
    checksums = {}
    for f in sorted(output_dir.glob("*")):
        if f.is_file() and f.suffix in (".parquet", ".csv", ".json", ".md"):
            checksums[f.name] = compute_sha256(f)

    checksum_path = output_dir / "SHA256SUMS.txt"
    with open(checksum_path, "w") as fh:
        for fname, chk in checksums.items():
            fh.write(f"{chk}  {fname}\n")
    return checksum_path


def write_zenodo_metadata(
    output_dir: Path,
    batch_label: str,
    counts: dict[str, int],
    version: str = "1.0",
) -> Path:
    """Write Zenodo deposition metadata JSON."""
    meta = {
        **DATASET_METADATA,
        "title": f"OligoToxDB Batch {batch_label} — {DATASET_METADATA['title']}",
        "version": version,
        "publication_date": str(date.today()),
        "notes": (
            f"Batch {batch_label} release. "
            f"Contains {counts.get('compounds', 0)} compounds, "
            f"{counts.get('results', 0)} result rows, "
            f"{counts.get('dose_response', 0)} dose-response summaries."
        ),
    }
    meta_path = output_dir / "zenodo_metadata.json"
    with open(meta_path, "w") as f:
        json.dump({"metadata": meta}, f, indent=2)
    return meta_path


def write_dataset_card(
    output_dir: Path,
    batch_label: str,
    counts: dict[str, int],
) -> Path:
    """Write a HuggingFace-compatible dataset card (dataset_card.md)."""
    card = f"""---
language:
  - en
license: cc-by-4.0
tags:
  - oligonucleotide
  - toxicity
  - drug-discovery
  - biology
  - genomics
task_categories:
  - tabular-regression
  - tabular-classification
size_categories:
  - 1K<n<10K
---

# OligoToxDB — Batch {batch_label}

## Dataset Description

**OligoToxDB** is the largest open dataset of human in vitro oligonucleotide toxicity profiles,
part of the NIH NCATS OligoTox Open Data Challenge (Phase 2).

This release (Batch {batch_label}) contains:
- **{counts.get("compounds", 0):,} oligonucleotides** with full chemical metadata
- **{counts.get("results", 0):,} experimental measurements** across 5 human in vitro systems
- **{counts.get("dose_response", 0):,} dose-response summaries** with IC50, Emax, Hill coefficient

### Supported Tasks

- **Regression**: Predict IC50 (µM) for hepatotoxicity, nephrotoxicity, immunotoxicity, complement, coagulopathy
- **Classification**: Classify oligos as inactive / low / moderate / high toxicity
- **Active learning**: Use uncertainty estimates to guide experimental compound selection

### Languages

N/A (scientific data)

## Dataset Structure

```
batch_{batch_label}/
├── compounds.parquet        # Compound registry (sequence, chemistry metadata)
├── features.parquet         # Pre-computed sequence/physicochemical features
├── results.parquet          # Raw dose-response measurements (all replicates)
├── dose_response.parquet    # IC50, Emax, Hill, NOEC, AUDRC per compound/endpoint
├── plate_qc.parquet         # Plate-level QC metrics (Z'-factor, pass/fail)
├── SHA256SUMS.txt           # File integrity checksums
└── zenodo_metadata.json     # Zenodo deposition metadata
```

## Data Fields

### compounds.parquet

| Column | Type | Description |
|---|---|---|
| oligo_id | string | Unique compound ID (OT-XXXXX) |
| sequence | string | Nucleotide sequence (5'→3') |
| backbone_class | string | PS / PO / PMO / LNA_mix / PNA / 2F_ANA / morpholino |
| sugar_mod | string | DNA / 2OMe / 2F / LNA / mixed |
| conjugate | string | none / GalNAc / cholesterol / lipid / antibody |
| purity_percent | float | HPLC purity (%) |
| endotoxin_eu_ml | float | Endotoxin level (EU/mL) |

### dose_response.parquet

| Column | Type | Description |
|---|---|---|
| oligo_id | string | Compound ID |
| assay_system | string | PHH / KidneyOrganoid / PBMC / Platelet / MPS |
| endpoint | string | Toxicity endpoint name |
| ic50_um | float | Half-maximal inhibitory concentration (µM) |
| emax_percent | float | Maximum effect at highest concentration (%) |
| hill_coefficient | float | Hill / slope coefficient |
| audrc | float | Area under dose-response curve |
| toxicity_class | string | inactive / low / moderate / high |
| r_squared | float | 4PL curve fit R² |

## Data Collection

All data generated using human in vitro model systems (no animal studies):
- **PHH**: Primary human hepatocytes (3 independent donors)
- **KidneyOrganoid**: Human kidney proximal tubule organoids (HUB certified)
- **PBMC**: Human peripheral blood mononuclear cells (5 donors)
- **Platelet**: Human platelet-rich plasma (5 donors)
- **MPS**: Liver-Kidney organ-on-chip (Emulate, Inc.)

## Licensing

Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
Code: [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)

## Citation

```bibtex
@dataset{{oligotoxdb_{batch_label.lower()},
  author    = {{Vidra, Natan and Anote, Inc.}},
  title     = {{OligoToxDB Batch {batch_label}}},
  year      = {{2026}},
  publisher = {{Zenodo}},
  license   = {{CC BY 4.0}},
  url       = {{https://zenodo.org/record/oligotoxdb-batch-{batch_label.lower()}}}
}}
```

## Contact

nvidra@anote.ai | [GitHub](https://github.com/anote-ai/nih-oligotox)
"""
    card_path = output_dir / "dataset_card.md"
    card_path.write_text(card)
    return card_path


def upload_to_zenodo(
    output_dir: Path,
    zenodo_token: str,
    sandbox: bool = True,
) -> Optional[str]:
    """Upload batch files to Zenodo via REST API. Returns deposition URL."""
    try:
        import requests
    except ImportError:
        console.print("[red]pip install requests to enable Zenodo upload[/red]")
        return None

    api_base = ZENODO_SANDBOX_API if sandbox else ZENODO_API
    headers = {"Authorization": f"Bearer {zenodo_token}"}

    # Create deposition
    meta_path = output_dir / "zenodo_metadata.json"
    with open(meta_path) as f:
        metadata = json.load(f)

    r = requests.post(f"{api_base}/deposit/depositions", json=metadata, headers=headers)
    r.raise_for_status()
    deposition = r.json()
    dep_id = deposition["id"]
    bucket_url = deposition["links"]["bucket"]
    console.print(f"Created Zenodo deposition {dep_id}")

    # Upload files
    for fpath in sorted(output_dir.glob("*.parquet")) + sorted(output_dir.glob("*.csv")):
        with open(fpath, "rb") as f:
            r = requests.put(f"{bucket_url}/{fpath.name}", data=f, headers=headers)
            r.raise_for_status()
        console.print(f"  Uploaded {fpath.name}")

    dep_url = deposition["links"]["html"]
    console.print(f"[green]Deposition ready: {dep_url}[/green]")
    return dep_url


def write_changelog_entry(
    batch_label: str,
    counts: dict[str, int],
    dep_url: Optional[str],
    changelog_path: Path = Path("CHANGELOG.md"),
) -> None:
    entry = f"""
## Batch {batch_label} — {date.today()}

- **Compounds**: {counts.get("compounds", 0):,}
- **Result rows**: {counts.get("results", 0):,}
- **Dose-response summaries**: {counts.get("dose_response", 0):,}
- **Zenodo**: {dep_url or "pending"}

"""
    if changelog_path.exists():
        existing = changelog_path.read_text()
        changelog_path.write_text(entry + existing)
    else:
        changelog_path.write_text(f"# OligoToxDB Changelog\n{entry}")


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.command()
@click.option("--batch-label", required=True, help="Batch identifier (A, B, C, D, E)")
@click.option("--db", default="oligotoxdb.duckdb", help="Path to OligoToxDB database")
@click.option("--output", required=True, help="Output directory for release package")
@click.option("--zenodo-token", default=None, envvar="ZENODO_TOKEN", help="Zenodo API token")
@click.option("--sandbox/--production", default=True, help="Use Zenodo sandbox (default) or production")
@click.option("--no-upload", is_flag=True, help="Package files only, do not upload to Zenodo")
def cli(
    batch_label: str,
    db: str,
    output: str,
    zenodo_token: Optional[str],
    sandbox: bool,
    no_upload: bool,
) -> None:
    """Package and release a completed OligoToxDB batch."""
    output_dir = Path(output)
    console.print(f"[bold]Packaging Batch {batch_label}...[/bold]")

    with OligoToxDB(db) as odb:
        counts = export_batch(odb, batch_label, output_dir)

    table = Table(title=f"Batch {batch_label} Contents")
    table.add_column("Table"); table.add_column("Rows")
    for k, v in counts.items():
        table.add_row(k, f"{v:,}")
    console.print(table)

    write_checksums(output_dir)
    write_zenodo_metadata(output_dir, batch_label, counts)
    write_dataset_card(output_dir, batch_label, counts)
    console.print(f"[green]Package ready: {output_dir}[/green]")

    dep_url = None
    if not no_upload and zenodo_token:
        env_label = "sandbox" if sandbox else "production"
        console.print(f"[bold]Uploading to Zenodo ({env_label})...[/bold]")
        dep_url = upload_to_zenodo(output_dir, zenodo_token, sandbox=sandbox)
    elif not no_upload:
        console.print("[yellow]No ZENODO_TOKEN found — skipping upload. Use --no-upload to suppress this warning.[/yellow]")

    write_changelog_entry(batch_label, counts, dep_url)
    console.print("[bold green]Release complete.[/bold green]")


if __name__ == "__main__":
    cli()
