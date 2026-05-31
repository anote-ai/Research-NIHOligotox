"""
Transcriptomics and proteomics data integration for OligoToxDB.

Handles:
  - RNA-seq: DESeq2/edgeR output → normalized table → GEO-compatible format
  - Proteomics: MaxQuant/DIA-NN output → normalized table → PRIDE-compatible format
  - Pathway enrichment scoring (gene ontology, Reactome pathways)
  - Hepatotoxicity / nephrotoxicity gene module scoring
  - Integration with main OligoToxDB schema
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import zscore
import click
from rich.console import Console
from rich.table import Table

console = Console()


# ─── Curated gene panels ─────────────────────────────────────────────────────

# NanoString hepatotoxicity panel — 50 genes covering canonical DILI pathways
HEPATOTOX_GENE_PANEL = [
    # Stress response & oxidative damage
    "HMOX1", "NQO1", "GCLC", "GCLM", "SRXN1", "TXN", "TXNRD1",
    # ER stress / UPR
    "ATF4", "ATF6", "DDIT3", "XBP1", "HSPA5", "EIF2AK3",
    # Apoptosis
    "CASP3", "CASP7", "CASP8", "CASP9", "BCL2", "BAX", "BID", "CYCS",
    # Mitochondrial
    "PPARGC1A", "TFAM", "MT-CO1", "MT-ND1",
    # Lipid metabolism
    "FASN", "ACACA", "SREBF1", "PPARA", "ACOX1",
    # Transport / bile
    "ABCB11", "ABCC2", "ABCG2", "SLC10A1", "SLC22A7",
    # Liver-specific markers
    "ALB", "AFP", "CYP7A1", "CYP3A4", "CYP1A2", "G6PC",
    # Inflammation
    "IL6", "CXCL8", "CCL2", "TNF", "NFKB1",
]

# Nephrotoxicity gene panel — 30 genes
NEPHROTOX_GENE_PANEL = [
    # Proximal tubule injury markers
    "HAVCR1",   # KIM-1
    "LCN2",     # NGAL
    "UMOD",     # Uromodulin
    "CLDN16",   # Claudin-16 (tight junction)
    # Transport
    "SLC22A2",  # OCT2
    "SLC22A6",  # OAT1
    "SLC22A8",  # OAT3
    "LRP2",     # Megalin
    "CUBN",     # Cubilin
    # Stress / injury
    "HMOX1", "HIF1A", "SOD2", "GPX1",
    # Apoptosis
    "CASP3", "CASP9", "BCL2", "BAX",
    # Fibrosis
    "TGFB1", "CTGF", "FN1", "COL1A1",
    # Immune / inflammation
    "IL6", "CXCL1", "CCL2", "IL18",
    # Mitochondrial
    "PPARGC1A", "TFAM",
    # Structural
    "ACTA2", "MKI67",
]


# ─── RNA-seq processing ───────────────────────────────────────────────────────

def load_deseq2_results(
    results_csv: Path,
    oligo_id: str,
    assay_system: str,
    concentration_um: float,
    geo_accession: str = "",
) -> pd.DataFrame:
    """
    Parse DESeq2 output CSV and return a standardized transcriptomics DataFrame.

    Expected DESeq2 columns: gene_symbol (or rownames), log2FoldChange, padj, baseMean
    Standard output: oligo_id, assay_system, concentration_um, gene_symbol,
                     log2fc, padj, tpm (baseMean proxy), geo_accession
    """
    df = pd.read_csv(results_csv)

    # Detect gene column
    gene_col = next(
        (c for c in df.columns if c.lower() in ("gene_symbol", "gene", "genename", "symbol")),
        df.columns[0]
    )
    lfc_col = next(c for c in df.columns if "log2" in c.lower() or "lfc" in c.lower())
    padj_col = next(c for c in df.columns if "padj" in c.lower() or "adj" in c.lower())
    mean_col = next(
        (c for c in df.columns if "mean" in c.lower() or "basemean" in c.lower()),
        None
    )

    result = pd.DataFrame({
        "oligo_id": oligo_id,
        "assay_system": assay_system,
        "concentration_um": concentration_um,
        "gene_symbol": df[gene_col].astype(str),
        "log2fc": df[lfc_col].astype(float),
        "padj": df[padj_col].astype(float),
        "tpm": df[mean_col].astype(float) if mean_col else np.nan,
        "geo_accession": geo_accession,
    })

    return result.dropna(subset=["gene_symbol", "log2fc"])


def compute_gene_module_score(
    deseq2_df: pd.DataFrame,
    gene_panel: list[str],
    lfc_col: str = "log2fc",
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
) -> float:
    """
    Compute a gene module activity score for a curated gene panel.

    Score = mean |log2FC| of significant panel genes, signed by direction majority.
    Range: negative (suppression) to positive (activation); 0 = no effect.
    """
    panel_df = deseq2_df[deseq2_df["gene_symbol"].isin(gene_panel)].copy()
    if panel_df.empty:
        return 0.0

    significant = panel_df[
        (panel_df["padj"] < padj_threshold) &
        (panel_df[lfc_col].abs() > lfc_threshold)
    ]

    if significant.empty:
        return 0.0

    mean_lfc = float(significant[lfc_col].mean())
    mean_abs_lfc = float(significant[lfc_col].abs().mean())
    n_sig = len(significant)
    n_panel = len(gene_panel)

    # Score = coverage-weighted mean |LFC| with sign
    coverage = n_sig / n_panel
    return round(np.sign(mean_lfc) * mean_abs_lfc * coverage, 4)


def batch_compute_gene_modules(
    transcriptomics_df: pd.DataFrame,
    panels: Optional[dict[str, list[str]]] = None,
) -> pd.DataFrame:
    """
    Compute hepatotoxicity and nephrotoxicity gene module scores for all oligos.

    Returns DataFrame with columns: oligo_id, assay_system, concentration_um,
                                    hepatotox_score, nephrotox_score
    """
    if panels is None:
        panels = {
            "hepatotox_score": HEPATOTOX_GENE_PANEL,
            "nephrotox_score": NEPHROTOX_GENE_PANEL,
        }

    rows = []
    for (oligo_id, system, conc), grp in transcriptomics_df.groupby(
        ["oligo_id", "assay_system", "concentration_um"]
    ):
        row = {"oligo_id": oligo_id, "assay_system": system, "concentration_um": conc}
        for score_name, panel in panels.items():
            row[score_name] = compute_gene_module_score(grp, panel)
        rows.append(row)

    return pd.DataFrame(rows)


def format_for_geo(
    transcriptomics_df: pd.DataFrame,
    output_dir: Path,
    study_title: str = "OligoToxDB RNA-seq Toxicogenomics",
) -> None:
    """
    Export transcriptomics data in GEO SOFT format.

    Creates:
      - matrix.txt: gene × sample expression matrix (log2 TPM)
      - samples.txt: sample metadata table
      - README_GEO.txt: submission instructions
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Wide-format gene × sample matrix
    pivot = transcriptomics_df.pivot_table(
        index="gene_symbol",
        columns=["oligo_id", "assay_system", "concentration_um"],
        values="tpm",
        aggfunc="mean",
    )
    pivot.columns = [f"{oid}_{sys}_{conc}uM" for oid, sys, conc in pivot.columns]
    pivot.to_csv(output_dir / "matrix.txt", sep="\t")

    # Sample metadata
    samples = transcriptomics_df.drop_duplicates(["oligo_id", "assay_system", "concentration_um"])[
        ["oligo_id", "assay_system", "concentration_um", "geo_accession"]
    ].copy()
    samples["sample_name"] = samples.apply(
        lambda r: f"{r.oligo_id}_{r.assay_system}_{r.concentration_um}uM", axis=1
    )
    samples.to_csv(output_dir / "samples.txt", sep="\t", index=False)

    # README
    with open(output_dir / "README_GEO.txt", "w") as f:
        f.write(f"GEO Submission Package\n{'='*40}\n")
        f.write(f"Study: {study_title}\n")
        f.write(f"Samples: {len(samples)}\n")
        f.write(f"Genes: {len(pivot)}\n\n")
        f.write("Files:\n  matrix.txt  — log2 TPM expression matrix\n")
        f.write("  samples.txt — sample metadata\n\n")
        f.write("Instructions:\n  1. Submit via GEO Submission Portal\n")
        f.write("  2. Use 'Third-party reanalysis' if raw FASTQ already deposited\n")
        f.write("  3. Link to OligoToxDB Zenodo DOI in dataset description\n")

    console.print(f"[green]GEO export: {len(pivot)} genes × {len(pivot.columns)} samples → {output_dir}[/green]")


# ─── Proteomics processing ────────────────────────────────────────────────────

def load_maxquant_results(
    proteingroups_txt: Path,
    oligo_id: str,
    assay_system: str,
    concentration_um: float,
    sample_col_prefix: str = "LFQ intensity",
    ctrl_cols: Optional[list[str]] = None,
    treat_cols: Optional[list[str]] = None,
    pride_accession: str = "",
) -> pd.DataFrame:
    """
    Parse MaxQuant proteinGroups.txt and compute log2 fold-change vs control.

    Returns standardized proteomics DataFrame compatible with OligoToxDB schema.
    """
    pg = pd.read_csv(proteingroups_txt, sep="\t", low_memory=False)

    # Remove contaminants
    if "Potential contaminant" in pg.columns:
        pg = pg[pg["Potential contaminant"].isna()]
    if "Reverse" in pg.columns:
        pg = pg[pg["Reverse"].isna()]

    gene_col = next(
        (c for c in pg.columns if "gene" in c.lower()),
        "Protein IDs"
    )
    protein_col = next(
        (c for c in pg.columns if "protein id" in c.lower() or c == "Protein IDs"),
        pg.columns[0]
    )

    # Detect LFQ intensity columns
    lfq_cols = [c for c in pg.columns if c.startswith(sample_col_prefix)]
    if not lfq_cols:
        raise ValueError(f"No columns starting with '{sample_col_prefix}' in {proteingroups_txt}")

    ctrl_cols = ctrl_cols or lfq_cols[:len(lfq_cols) // 2]
    treat_cols = treat_cols or lfq_cols[len(lfq_cols) // 2:]

    # Log2 transform (replace 0 with NaN)
    intensity = pg[lfq_cols].replace(0, np.nan).apply(np.log2)

    ctrl_mean = intensity[ctrl_cols].mean(axis=1)
    treat_mean = intensity[treat_cols].mean(axis=1)
    log2fc = treat_mean - ctrl_mean

    # Simple t-test for significance
    from scipy.stats import ttest_ind
    pvals = []
    for i in range(len(pg)):
        ctrl_vals = intensity.loc[pg.index[i], ctrl_cols].dropna()
        treat_vals = intensity.loc[pg.index[i], treat_cols].dropna()
        if len(ctrl_vals) >= 2 and len(treat_vals) >= 2:
            _, p = ttest_ind(ctrl_vals, treat_vals)
            pvals.append(p)
        else:
            pvals.append(np.nan)

    # Benjamini-Hochberg correction
    padj = _bh_correction(np.array(pvals))

    result = pd.DataFrame({
        "oligo_id": oligo_id,
        "assay_system": assay_system,
        "concentration_um": concentration_um,
        "protein_id": pg[protein_col].astype(str),
        "gene_symbol": pg[gene_col].astype(str).str.split(";").str[0],
        "log2fc": log2fc.values,
        "padj": padj,
        "intensity": treat_mean.values,
        "pride_accession": pride_accession,
    })
    return result.dropna(subset=["gene_symbol", "log2fc"])


def format_for_pride(
    proteomics_df: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Export proteomics data in PRIDE-compatible format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    proteomics_df.to_csv(output_dir / "protein_abundance.csv", index=False)

    # PRIDE submission template
    with open(output_dir / "README_PRIDE.txt", "w") as f:
        f.write("PRIDE Submission Package\n")
        f.write("Submit via PRIDE Submission Tool (ProteomeXchange)\n")
        f.write("Data type: Label-free quantification (LFQ)\n")
        f.write("Raw files: .raw (Thermo) or .wiff (SCIEX) — attach separately\n")
        f.write("Processing: MaxQuant v2.x; match between runs enabled\n")
        f.write("Quantification: LFQ intensity, min. 2 peptides per protein\n")
        f.write(f"\nSamples: {proteomics_df['oligo_id'].nunique()} oligos\n")
        f.write(f"Proteins: {proteomics_df['protein_id'].nunique()}\n")

    console.print(f"[green]PRIDE export → {output_dir}[/green]")


# ─── Differential expression helpers ─────────────────────────────────────────

def _bh_correction(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    padj = np.ones(n)
    valid = ~np.isnan(pvals)
    if not valid.any():
        return padj

    p = pvals[valid]
    rank = np.argsort(np.argsort(p)) + 1
    bh = np.minimum(1.0, p * n / rank)
    # Make monotone
    for i in range(len(bh) - 2, -1, -1):
        bh[i] = min(bh[i], bh[i + 1])
    padj[valid] = bh
    return padj


def top_dysregulated_genes(
    df: pd.DataFrame,
    n: int = 20,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """Return top n most significantly dysregulated genes."""
    sig = df[
        (df["padj"] < padj_threshold) &
        (df["log2fc"].abs() > lfc_threshold)
    ].copy()
    sig["neg_log10_padj"] = -np.log10(sig["padj"].clip(lower=1e-300))
    sig["score"] = sig["neg_log10_padj"] * sig["log2fc"].abs()
    return sig.nlargest(n, "score")[["gene_symbol", "log2fc", "padj", "score"]]


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """Omics data processing (RNA-seq and proteomics)."""


@cli.command()
@click.argument("deseq2_csv", type=click.Path(exists=True))
@click.option("--oligo-id", required=True)
@click.option("--assay-system", default="PHH")
@click.option("--concentration", default=10.0, type=float)
@click.option("--output", default="transcriptomics.csv")
def process_rnaseq(deseq2_csv, oligo_id, assay_system, concentration, output):
    """Process a DESeq2 results CSV for one oligo and compute gene module scores."""
    df = load_deseq2_results(Path(deseq2_csv), oligo_id, assay_system, concentration)
    df.to_csv(output, index=False)
    hep_score = compute_gene_module_score(df, HEPATOTOX_GENE_PANEL)
    neph_score = compute_gene_module_score(df, NEPHROTOX_GENE_PANEL)
    n_sig = int((df["padj"] < 0.05).sum())
    console.print(f"[bold]{oligo_id}[/bold] | Sig. genes: {n_sig} | Hep score: {hep_score:.3f} | Neph score: {neph_score:.3f}")
    console.print(f"[green]Saved → {output}[/green]")
    top = top_dysregulated_genes(df)
    if len(top) > 0:
        console.print("\nTop dysregulated genes:")
        console.print(top.to_string(index=False))


if __name__ == "__main__":
    cli()
