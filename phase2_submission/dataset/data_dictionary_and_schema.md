# OligoTox Open Data Challenge — Phase 2 Dataset: Data Dictionary and Schema

**Submitter:** Anote, Inc. | No page limit per challenge rules, but must be complete and precise.

> **STATUS: SCHEMA DEFINED, POPULATED WITH SYNTHETIC PLACEHOLDER DATA ONLY.**
> The schema below is implemented and exercised end-to-end by the `oligotoxdb` codebase in this
> repository (see `oligotoxdb/oligotoxdb/`, `oligotoxdb/scripts/generate_synthetic_data.py`).
> However, the data currently produced by that pipeline (`oligotoxdb/data/synthetic/*.csv`) is
> **explicitly synthetic** — generated for pipeline testing, not from real wet-lab experiments.
> **This file must be repointed at real experimental output before submission.** Submitting the
> synthetic files as the Phase 2 dataset would misrepresent fabricated data as real results.

## 1. Required Contents Per Challenge Rules

The dataset must include, at minimum:
- Sequences of all oligos tested
- Location of all chemical modifications in each oligo
- Purity and characterization data for each oligo
- Any additional metadata
- A data dictionary and schema documenting all of the above
- Access to the raw data (file included directly, or documented access/download instructions)

## 2. Schema

### Table 1 — Compound Registry

| Field | Type | Description |
|---|---|---|
| `OligoTox_ID` | string | Unique identifier, format `OT-XXXXX` |
| `Sequence_5to3` | string | Full nucleotide sequence, 5'→3' |
| `Length` | integer | Sequence length (nt) |
| `Backbone_class` | enum | PS, PO, PMO, LNA, PNA, 2F-ANA, morpholino |
| `Sugar_mod` | enum | DNA, 2'OMe, 2'F, LNA, mixed |
| `Modification_positions` | string | Position-indexed map of every chemical modification along the sequence (**required by challenge rules** — must specify *location*, not just an overall class label) |
| `Gapmer_design` | string | e.g. "5-10-5 LNA gapmer"; blank if not a gapmer |
| `Conjugate` | enum | none, GalNAc, cholesterol, lipid, antibody |
| `Synthesis_vendor` | string | CRO that synthesized the compound |
| `Purity_percent` | float | RP-HPLC purity (%) |
| `Mass_confirmed` | boolean | LC-MS identity confirmation pass/fail |
| `Mass_observed_Da` | float | Observed mass from MS |
| `Mass_expected_Da` | float | Expected mass from sequence + chemistry |
| `Endotoxin_EU_mL` | float | Endotoxin level in working solution |
| `Control_type` | enum | none, positive, negative, reference | 

### Table 2 — Computed Features (predictor variables)

| Field | Type | Description |
|---|---|---|
| `OligoTox_ID` | string | Foreign key |
| `GC_content` | float | % GC |
| `CpG_count` | integer | Count of CpG dinucleotides |
| `MFE_kcal_mol` | float | Minimum free energy (RNAfold) |
| `Tm_target` | float | Melting temp vs. intended target (°C) |
| `Tm_offtarget_max` | float | Melting temp vs. top off-target BLAST hit (°C) |
| `G4_score` | float | G-quadruplex propensity (QGRS) |
| `Self_comp_index` | float | Self-complementarity index |
| `Hybridization_dG` | float | Predicted hybridization ΔG (kcal/mol) |
| `Offtarget_score` | float | log10 BLAST E-value, top off-target |
| `MW_Da` | float | Molecular weight |
| `Net_charge_pH7` | integer | Net charge at physiological pH |

### Table 3 — Experimental Results (indicators, raw + normalized)

| Field | Type | Description |
|---|---|---|
| `OligoTox_ID` | string | Foreign key |
| `Assay_system` | enum | PHH, KidneyOrganoid, PBMC, Platelet, MPS |
| `Endpoint` | string | Readout name (controlled vocabulary; see endpoint registry) |
| `Concentration_uM` | float | Dosing concentration |
| `Timepoint_h` | integer | Exposure duration (hours) |
| `Replicate_bio` | integer | Biological replicate index |
| `Replicate_tech` | integer | Technical replicate index |
| `Batch_ID` | string | Experimental batch identifier |
| `Cell_donor_ID` | string | Anonymized donor/lot identifier |
| `Value_raw` | float | Raw instrument reading |
| `Value_normalized` | float | % of vehicle control |
| `Unit` | string | Unit of raw value |
| `QC_flag` | enum | pass, marginal, fail |

### Table 4 — Derived Summary Statistics

| Field | Type | Description |
|---|---|---|
| `OligoTox_ID` | string | Foreign key |
| `Assay_system` | string | Foreign key |
| `Endpoint` | string | Foreign key |
| `IC50_uM` | float | Fitted IC50 |
| `IC50_CI_low` / `IC50_CI_high` | float | 95% CI bounds |
| `Emax_percent` | float | Maximum effect observed |
| `Hill_coefficient` | float | Curve slope/cooperativity |
| `NOEC_uM` | float | No observed effect concentration |
| `AUDRC` | float | Area under dose-response curve |
| `N_bio_reps` | integer | Number of biological replicates |
| `Toxicity_class` | enum | inactive, low, moderate, high |

## 3. Controlled Vocabularies / Ontologies

- Chemical entities: ChEBI (+ oligo-specific extensions)
- Cell types: Cell Ontology (CL)
- Assay types: Bioassay Ontology (BAO)
- Toxicity endpoints: NCI Thesaurus (NCIt)
- Genes/transcripts: HGNC symbols, Ensembl IDs

## 4. File Formats Provided

| Content | Format | License |
|---|---|---|
| Compound registry, features, results, summary stats | CSV / Parquet | CC BY 4.0 |
| Raw transcriptomics (if applicable) | FASTQ + TPM matrix | CC BY 4.0 |
| Raw proteomics (if applicable) | mzML + abundance CSV | CC BY 4.0 |

## 5. Access to Raw Data

See `README_raw_data_access.md` in this folder for exact download links/DOIs (to be filled in
once the real dataset is deposited) or, if provided as attached files, the file names/checksums
listed there.

## 6. Reference Implementation

The schema above is implemented in code at:
- `oligotoxdb/oligotoxdb/database.py`, `oligotoxdb/oligotoxdb/endpoints.py` (schema + endpoint
  registry)
- `oligotoxdb/oligotoxdb/qc.py` (QC flag computation)
- `oligotoxdb/oligotoxdb/ingestion.py` (loading raw instrument files into this schema)

**Before submission:** point ingestion at real instrument/LIMS output (replacing
`oligotoxdb/data/synthetic/`) and regenerate `compounds.csv`, `features.csv`, `results.csv`, and
`plate_controls.csv` (or equivalent) from actual experimental data.
