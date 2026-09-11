# Raw Data Access Instructions — OligoTox Phase 2 Submission

> **STATUS: PLACEHOLDER — no real dataset has been deposited yet.** Per the challenge rules,
> the dataset must be provided either as an included data file (Excel/CSV or similar) or as
> instructions for the challenge sponsors on how to access and download the raw data. Fill in
> the section below that applies once real data exists, and delete the other.

## Option A — Data included directly in the submission package

If the final dataset is small enough to attach directly:

- `compounds.csv` — Compound Registry (Table 1)
- `features.csv` — Computed Features (Table 2)
- `results.csv` — Experimental Results, raw + normalized (Table 3)
- `summary_stats.csv` — Derived Summary Statistics (Table 4)
- `plate_controls.csv` — Per-plate QC/control performance records

*(Replace with actual real-data file names/paths once generated; current repository placeholders
live at `oligotoxdb/data/synthetic/` but are SYNTHETIC and must not be submitted as-is.)*

## Option B — Hosted repository with download instructions

If deposited externally (recommended for large files, e.g. raw FASTQ/mzML):

| Data type | Repository | DOI / URL | Format |
|---|---|---|---|
| Compound registry + assay results | Zenodo | *[DOI once deposited]* | CSV, Parquet |
| Raw transcriptomics | GEO (NCBI) | *[GEO accession once deposited]* | FASTQ, TPM matrix |
| Raw proteomics | PRIDE Archive (EBI) | *[PXD accession once deposited]* | mzML, abundance CSV |

Access is open (no login required) and licensed CC BY 4.0 — see
`../03_public_access_and_dissemination_plan.md`.

## Checksums

*[Add SHA-256 checksums for each delivered file once final, so NCATS can verify integrity of
downloaded/attached data.]*

## Contact for Data Access Issues

Natan Vidra, Anote, Inc. — nvidra@anote.ai
