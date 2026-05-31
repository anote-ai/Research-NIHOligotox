"""
OligoTox Minimum Reporting Standard (OTMRS) Validator

External labs contributing data to OligoToxDB must comply with OTMRS.
This module validates submission files and generates a compliance report.

Specification published alongside the dataset (see Appendix C of Phase 2 submission).

Usage:
    python -m oligotoxdb.otmrs_validator compounds.csv results.csv
    oligotox-validate compounds.csv results.csv --output report.html
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
import click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

# ─── Required fields (OTMRS spec) ────────────────────────────────────────────

REQUIRED_COMPOUND_FIELDS = {
    "oligo_id": "Unique compound identifier",
    "sequence": "Nucleotide sequence (5'→3')",
    "backbone_class": "Backbone chemistry (PS/PO/PMO/LNA_mix/PNA/2F_ANA/morpholino)",
    "sugar_mod": "Sugar modification (DNA/2OMe/2F/LNA/mixed)",
    "conjugate": "Conjugate type (none/GalNAc/cholesterol/lipid/antibody)",
    "synthesis_vendor": "Oligonucleotide synthesis vendor",
    "purity_percent": "Purity by HPLC (%)",
    "mass_confirmed": "Mass confirmation by MS (True/False)",
    "endotoxin_eu_ml": "Endotoxin level (EU/mL)",
}

REQUIRED_RESULT_FIELDS = {
    "oligo_id": "Compound identifier (must match compound registry)",
    "assay_system": "Assay platform (PHH/KidneyOrganoid/PBMC/Platelet/MPS)",
    "endpoint": "Measured endpoint",
    "concentration_um": "Compound concentration (µM)",
    "timepoint_h": "Incubation duration (hours)",
    "replicate_bio": "Biological replicate number (≥1)",
    "replicate_tech": "Technical replicate number (≥1)",
    "batch_id": "Experimental batch identifier",
    "plate_id": "Plate identifier",
    "cell_donor_id": "Anonymized cell donor identifier",
    "value_raw": "Raw measurement value",
    "unit": "Measurement unit",
}

VALID_VALUES = {
    "backbone_class": {"PS", "PO", "PMO", "LNA_mix", "PNA", "2F_ANA", "morpholino"},
    "sugar_mod": {"DNA", "2OMe", "2F", "LNA", "mixed"},
    "conjugate": {"none", "GalNAc", "cholesterol", "lipid", "antibody"},
    "assay_system": {"PHH", "KidneyOrganoid", "PBMC", "Platelet", "MPS"},
}

VALID_NUCLEOTIDES = set("ACGTUacgtu")
MIN_PURITY = 85.0
MAX_ENDOTOXIN = 10.0  # EU/mL (strict), 1.0 for Tier 1
MIN_BIO_REPLICATES = 2
MIN_CONCENTRATIONS = 3


# ─── Validation result types ─────────────────────────────────────────────────

@dataclass
class ValidationError:
    field: str
    row: Optional[int]
    value: Optional[str]
    message: str
    severity: str = "error"  # "error" | "warning"


@dataclass
class OTMRSReport:
    compound_errors: list[ValidationError] = field(default_factory=list)
    result_errors: list[ValidationError] = field(default_factory=list)
    n_compounds: int = 0
    n_results: int = 0
    n_endpoints: int = 0
    n_bio_replicates: int = 0
    n_concentrations_per_compound: float = 0.0
    passed: bool = False

    @property
    def n_errors(self) -> int:
        return sum(1 for e in self.compound_errors + self.result_errors if e.severity == "error")

    @property
    def n_warnings(self) -> int:
        return sum(1 for e in self.compound_errors + self.result_errors if e.severity == "warning")


# ─── Compound validation ─────────────────────────────────────────────────────

def validate_compounds(df: pd.DataFrame) -> list[ValidationError]:
    errors: list[ValidationError] = []

    # Required fields
    for field_name, desc in REQUIRED_COMPOUND_FIELDS.items():
        if field_name not in df.columns:
            errors.append(ValidationError(field_name, None, None, f"Required column missing: {field_name} ({desc})", "error"))

    if errors:
        return errors  # can't validate further without columns

    # Row-level validation
    for idx, row in df.iterrows():
        i = int(idx)

        # Sequence validation
        seq = str(row.get("sequence", ""))
        if not seq or not all(c in VALID_NUCLEOTIDES for c in seq):
            errors.append(ValidationError("sequence", i, seq[:30], "Invalid nucleotides or empty sequence", "error"))

        if len(seq) < 6 or len(seq) > 50:
            errors.append(ValidationError("sequence", i, str(len(seq)), "Sequence length outside 6-50 nt range", "warning"))

        # Controlled vocabulary fields
        for field_name, valid_set in VALID_VALUES.items():
            if field_name not in df.columns:
                continue
            val = str(row.get(field_name, ""))
            if val not in valid_set:
                errors.append(ValidationError(
                    field_name, i, val,
                    f"Value '{val}' not in allowed set: {sorted(valid_set)}",
                    "error"
                ))

        # Purity
        try:
            purity = float(row.get("purity_percent", 0))
            if purity < MIN_PURITY:
                errors.append(ValidationError(
                    "purity_percent", i, str(purity),
                    f"Purity {purity:.1f}% below minimum {MIN_PURITY}%",
                    "error" if purity < 80 else "warning"
                ))
        except (ValueError, TypeError):
            errors.append(ValidationError("purity_percent", i, str(row.get("purity_percent")), "Must be numeric", "error"))

        # Endotoxin
        try:
            endo = float(row.get("endotoxin_eu_ml", np.inf))
            if endo > MAX_ENDOTOXIN:
                errors.append(ValidationError(
                    "endotoxin_eu_ml", i, str(endo),
                    f"Endotoxin {endo:.2f} EU/mL exceeds maximum {MAX_ENDOTOXIN} EU/mL",
                    "warning"
                ))
        except (ValueError, TypeError):
            errors.append(ValidationError("endotoxin_eu_ml", i, str(row.get("endotoxin_eu_ml")), "Must be numeric", "error"))

        # Mass confirmation
        mass = str(row.get("mass_confirmed", "")).lower()
        if mass not in ("true", "false", "1", "0", "yes", "no"):
            errors.append(ValidationError("mass_confirmed", i, mass, "Must be True/False", "warning"))

        # Duplicate oligo_id check
    dupes = df["oligo_id"].duplicated()
    if dupes.any():
        for idx in df[dupes].index:
            errors.append(ValidationError(
                "oligo_id", int(idx), str(df.loc[idx, "oligo_id"]),
                "Duplicate oligo_id", "error"
            ))

    return errors


# ─── Results validation ───────────────────────────────────────────────────────

def validate_results(
    results_df: pd.DataFrame,
    compound_ids: Optional[set[str]] = None,
) -> list[ValidationError]:
    errors: list[ValidationError] = []

    # Required fields
    for field_name in REQUIRED_RESULT_FIELDS:
        if field_name not in results_df.columns:
            errors.append(ValidationError(field_name, None, None, f"Required column missing: {field_name}", "error"))

    if errors:
        return errors

    # Cross-reference with compound registry
    if compound_ids is not None:
        unknown = set(results_df["oligo_id"].unique()) - compound_ids
        for uid in unknown:
            errors.append(ValidationError(
                "oligo_id", None, uid,
                f"oligo_id '{uid}' not found in compound registry",
                "error"
            ))

    # Assay system validation
    for idx, row in results_df.iterrows():
        i = int(idx)
        system = str(row.get("assay_system", ""))
        if system not in VALID_VALUES["assay_system"]:
            errors.append(ValidationError("assay_system", i, system, f"Invalid assay_system: {system}", "error"))

        conc = row.get("concentration_um")
        try:
            if float(conc) < 0:
                errors.append(ValidationError("concentration_um", i, str(conc), "Must be ≥ 0", "error"))
        except (ValueError, TypeError):
            errors.append(ValidationError("concentration_um", i, str(conc), "Must be numeric", "error"))

        val = row.get("value_raw")
        try:
            float(val)
        except (ValueError, TypeError):
            errors.append(ValidationError("value_raw", i, str(val), "Must be numeric", "error"))

    # Replicate sufficiency check
    bio_rep_counts = results_df.groupby(["oligo_id", "endpoint", "concentration_um"])["replicate_bio"].nunique()
    insufficient = bio_rep_counts[bio_rep_counts < MIN_BIO_REPLICATES]
    for (oligo_id, endpoint, conc), n in insufficient.items():
        errors.append(ValidationError(
            "replicate_bio", None, str(n),
            f"{oligo_id}/{endpoint}/{conc}µM: only {n} bio replicate(s) (min {MIN_BIO_REPLICATES})",
            "warning"
        ))

    # Concentration coverage
    conc_counts = results_df.groupby("oligo_id")["concentration_um"].nunique()
    insufficient_conc = conc_counts[conc_counts < MIN_CONCENTRATIONS]
    for oligo_id, n in insufficient_conc.items():
        errors.append(ValidationError(
            "concentration_um", None, str(n),
            f"{oligo_id}: only {n} concentration(s) (min {MIN_CONCENTRATIONS} for dose-response)",
            "warning"
        ))

    return errors


# ─── Report generation ────────────────────────────────────────────────────────

def validate_submission(
    compounds_csv: Path,
    results_csv: Path,
) -> OTMRSReport:
    """Run full OTMRS validation on a submission package."""
    report = OTMRSReport()

    compounds_df = pd.read_csv(compounds_csv)
    results_df = pd.read_csv(results_csv)

    report.n_compounds = len(compounds_df)
    report.n_results = len(results_df)
    report.n_endpoints = results_df["endpoint"].nunique() if "endpoint" in results_df.columns else 0
    report.compound_errors = validate_compounds(compounds_df)

    compound_ids = set(compounds_df["oligo_id"].dropna().astype(str)) if "oligo_id" in compounds_df.columns else None
    report.result_errors = validate_results(results_df, compound_ids)

    if "replicate_bio" in results_df.columns:
        report.n_bio_replicates = int(results_df["replicate_bio"].max())
    if "concentration_um" in results_df.columns and "oligo_id" in results_df.columns:
        report.n_concentrations_per_compound = float(
            results_df.groupby("oligo_id")["concentration_um"].nunique().mean()
        )

    report.passed = report.n_errors == 0
    return report


def print_report(report: OTMRSReport, title: str = "OTMRS Validation Report") -> None:
    status = "[bold green]PASSED[/bold green]" if report.passed else "[bold red]FAILED[/bold red]"
    console.print(Panel(f"[bold]{title}[/bold]\nStatus: {status}", expand=False))

    # Summary table
    table = Table(title="Dataset Summary")
    table.add_column("Metric"); table.add_column("Value")
    table.add_row("Compounds", str(report.n_compounds))
    table.add_row("Result rows", str(report.n_results))
    table.add_row("Endpoints", str(report.n_endpoints))
    table.add_row("Max bio replicates", str(report.n_bio_replicates))
    table.add_row("Mean concentrations/compound", f"{report.n_concentrations_per_compound:.1f}")
    table.add_row("[red]Errors[/red]", str(report.n_errors))
    table.add_row("[yellow]Warnings[/yellow]", str(report.n_warnings))
    console.print(table)

    if report.compound_errors:
        console.print(f"\n[bold]Compound Errors/Warnings ({len(report.compound_errors)}):[/bold]")
        for e in report.compound_errors[:20]:
            color = "red" if e.severity == "error" else "yellow"
            row_str = f" (row {e.row})" if e.row is not None else ""
            console.print(f"  [{color}]{e.severity.upper()}[/{color}] [{e.field}]{row_str}: {e.message}")
        if len(report.compound_errors) > 20:
            console.print(f"  ... and {len(report.compound_errors) - 20} more")

    if report.result_errors:
        console.print(f"\n[bold]Result Errors/Warnings ({len(report.result_errors)}):[/bold]")
        for e in report.result_errors[:20]:
            color = "red" if e.severity == "error" else "yellow"
            row_str = f" (row {e.row})" if e.row is not None else ""
            console.print(f"  [{color}]{e.severity.upper()}[/{color}] [{e.field}]{row_str}: {e.message}")
        if len(report.result_errors) > 20:
            console.print(f"  ... and {len(report.result_errors) - 20} more")


def generate_html_report(report: OTMRSReport, output_path: Path) -> None:
    """Write a self-contained HTML validation report."""
    status_color = "#2ecc71" if report.passed else "#e74c3c"
    status_text = "PASSED" if report.passed else "FAILED"

    all_errors = report.compound_errors + report.result_errors
    error_rows = "\n".join(
        f"<tr class='{e.severity}'>"
        f"<td>{e.severity.upper()}</td>"
        f"<td>{e.field}</td>"
        f"<td>{e.row if e.row is not None else 'N/A'}</td>"
        f"<td>{e.value or ''}</td>"
        f"<td>{e.message}</td>"
        f"</tr>"
        for e in all_errors
    )

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>OTMRS Validation Report</title>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1000px; margin: 40px auto; padding: 0 20px; }}
  h1 {{ color: #2c3e50; }}
  .status {{ font-size: 2em; font-weight: bold; color: {status_color}; }}
  table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; font-size: 0.9em; }}
  th {{ background: #f5f5f5; }}
  tr.error {{ background: #fdecea; }}
  tr.warning {{ background: #fff8e1; }}
  .summary td:first-child {{ font-weight: bold; }}
</style>
</head><body>
<h1>OligoTox Minimum Reporting Standard (OTMRS) Validation Report</h1>
<p class="status">{status_text}</p>
<h2>Dataset Summary</h2>
<table class="summary">
  <tr><td>Compounds</td><td>{report.n_compounds}</td></tr>
  <tr><td>Result rows</td><td>{report.n_results}</td></tr>
  <tr><td>Endpoints</td><td>{report.n_endpoints}</td></tr>
  <tr><td>Max biological replicates</td><td>{report.n_bio_replicates}</td></tr>
  <tr><td>Mean concentrations/compound</td><td>{report.n_concentrations_per_compound:.1f}</td></tr>
  <tr><td>Errors</td><td style="color: #e74c3c; font-weight: bold">{report.n_errors}</td></tr>
  <tr><td>Warnings</td><td style="color: #f39c12; font-weight: bold">{report.n_warnings}</td></tr>
</table>
<h2>Validation Findings</h2>
<table>
  <tr><th>Severity</th><th>Field</th><th>Row</th><th>Value</th><th>Message</th></tr>
  {error_rows if error_rows else "<tr><td colspan='5' style='color:green'>No issues found</td></tr>"}
</table>
<hr><p style="color:#888;font-size:0.8em">
Generated by OligoToxDB OTMRS Validator | anote-ai/nih-oligotox
</p>
</body></html>"""

    output_path = Path(output_path)
    output_path.write_text(html)
    console.print(f"[green]HTML report saved → {output_path}[/green]")


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.command()
@click.argument("compounds_csv", type=click.Path(exists=True))
@click.argument("results_csv", type=click.Path(exists=True))
@click.option("--output", "-o", default=None, help="Optional HTML report output path")
@click.option("--strict", is_flag=True, help="Treat warnings as errors")
def cli(compounds_csv: str, results_csv: str, output: Optional[str], strict: bool) -> None:
    """Validate a submission package against the OTMRS standard."""
    console.print(f"[bold]Validating OTMRS submission...[/bold]")
    report = validate_submission(Path(compounds_csv), Path(results_csv))

    if strict:
        report.passed = (report.n_errors + report.n_warnings) == 0

    print_report(report)

    if output:
        generate_html_report(report, Path(output))

    sys.exit(0 if report.passed else 1)


if __name__ == "__main__":
    cli()
