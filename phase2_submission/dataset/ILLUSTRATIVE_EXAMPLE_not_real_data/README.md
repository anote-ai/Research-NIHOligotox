# ⚠️ SYNTHETIC ILLUSTRATIVE EXAMPLE — NOT REAL DATA — DO NOT SUBMIT ⚠️

Every file in this folder is **fabricated by code**
(`oligotoxdb/scripts/generate_synthetic_data.py`, `oligotoxdb/oligotoxdb/features.py`,
`oligotoxdb/oligotoxdb/qc.py`), not measured from any real oligonucleotide, cell, or assay. No
compound described here was ever synthesized. No cell was ever dosed. Every numeric value —
purity, mass, cytokine level, IC50 — is a randomly generated number shaped to look plausible.

**Every row in every file carries a `DATA_STATUS = SYNTHETIC_ILLUSTRATIVE_NOT_REAL` column so
this cannot be silently mixed into a real dataset later.** Do not delete that column. Do not
strip this README from the folder. Do not attach these files to the NIH submission.

## What this is for

This is a small (15-compound), fully-worked example showing that the OligoToxDB pipeline —
compound registry → computed features → raw assay results → QC → dose-response fitting —
actually runs end-to-end and produces output in exactly the shape described in
`../data_dictionary_and_schema.md`. It exists so you (and anyone reviewing this repo) can see
the intended file formats before real data exists, and so that once real data is collected, the
existing ingestion code (`oligotoxdb/oligotoxdb/ingestion.py`) can be pointed at it with no
schema changes.

## Files

| File | Rows | Corresponds to |
|---|---|---|
| `compounds.csv` | 15 | Data dictionary Table 1 — Compound Registry |
| `features.csv` | 15 | Data dictionary Table 2 — Computed Features |
| `results.csv` | 16,200 | Data dictionary Table 3 — Experimental Results (raw) |
| `plate_controls.csv` | 6,912 | Per-plate positive/negative control readings (QC input) |
| `plate_qc.csv` | 130 | Per-plate QC pass/fail (Z'-factor, control CV) |
| `dose_response_summary.csv` | 450 | Data dictionary Table 4 — Derived Summary Statistics (IC50, Emax, Hill, NOEC, AUDRC) |

## How it was generated (for full transparency)

```
python oligotoxdb/scripts/generate_synthetic_data.py \
    --n-oligos 15 --tier 3 --seed 7 \
    --output phase2_submission/dataset/ILLUSTRATIVE_EXAMPLE_not_real_data
# + oligotoxdb.features.compute_features_batch()
# + oligotoxdb.qc.run_qc_pipeline()
```

Regenerate at any time with a different `--seed`/`--n-oligos` — it's deterministic code, not a
one-off file, precisely so nobody is tempted to treat a specific run as precious/authoritative.

## What replaces this before submission

Real data covering the same four tables, produced by actually synthesizing oligonucleotides
(with real HPLC/MS characterization) and actually running them through a real in vitro assay —
see `../../02_methodology_document.md` for the protocol this pipeline expects real data to match.
