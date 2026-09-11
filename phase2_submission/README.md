# Phase 2 Submission Package — Status

This folder splits the previously combined `OligoTox_Phase2_Submission_Anote.md` into the four
parts NCATS actually requires, formats them to NIH's PDF spec, and tracks what's real vs.
placeholder.

| # | Required part | Source | Formatted PDF | Status |
|---|---|---|---|---|
| 1 | Narrative document (≤12pp) | `01_narrative_document.md` | `pdf/01_narrative_document.pdf` (6pp) | **Finalized:** Sections 3, 5, 6 (process, gap analysis, modeling plan). **Pending real data:** Sections 1, 2, 4 — not submittable until filled from actual results |
| 2 | Methodology document (≤5pp) | `02_methodology_document.md` | `pdf/02_methodology_document.pdf` (3pp) | **Finalized** — concrete protocol (vendors, concentrations, timepoints, replicates); only needs a check against actual execution if anything is dropped/changed in practice |
| 3 | Public Access & Dissemination Plan (≤5pp) | `03_public_access_and_dissemination_plan.md` | `pdf/03_public_access_and_dissemination_plan.pdf` (3pp) | **Finalized** — submittable as-is |
| 4 | Dataset (schema + raw data, no page limit) | `dataset/data_dictionary_and_schema.md`, `dataset/README_raw_data_access.md` | — (data files, not a PDF) | Schema fully defined and pipeline-tested; **real experimental data does not exist yet — the only remaining placeholder** |
| — | Registration form | Download from Challenge.gov "Resources" tab | — | Not started — external form, not generated here |
| — | Submission email | `email_draft.md` | — | Template ready with `nvidra@anote.ai`; **do not send until the dataset is real** |

## PDF format compliance

`pdf/*.pdf` were generated with `tools/build_pdfs.py` (reportlab) to meet the announcement's
format rules: 8.5"×11" page size, 1" margins, 11pt font, line spacing >1.0, no HHS/NIH/NCATS
logos. One deviation: the announcement specifies **Arial**; no Arial font file was available in
this build environment, so **Helvetica** (a metrically near-identical substitute) was used
instead. If exact-font compliance matters, re-export via Word/Google Docs with Arial, or add a
real Arial `.ttf` and register it in `tools/build_pdfs.py` (`FONT`/`FONT_BOLD` constants) before
rebuilding.

Regenerate PDFs any time after editing the `.md` sources:
```
python3 phase2_submission/tools/build_pdfs.py
```

## Example dataset

`dataset/ILLUSTRATIVE_EXAMPLE_not_real_data/` is a small (15-compound), fully-worked example
showing the OligoToxDB pipeline running end-to-end (compound registry → computed features → raw
assay results → QC → dose-response fit) in exactly the file shapes the data dictionary
describes. **It is synthetic, fabricated by code, and watermarked in every row
(`DATA_STATUS = SYNTHETIC_ILLUSTRATIVE_NOT_REAL`) — it is not real data and must not be
submitted.** See that folder's own README for full detail on why it exists and what replaces it.

## The one thing left — real data

Every document (narrative process/gap/modeling sections, methodology, PADP) is now finalized.
What's left is exactly what you said you'd handle on your side: the dataset itself. NIH judges
"the demonstrated dataset as it is submitted" — real sequences actually synthesized, real
purity/MS characterization, and real toxicity readouts from real human in vitro experiments. As
of this draft, no wet-lab data generation has occurred.

**Once real data exists:**
1. Ingest real instrument/LIMS output through the existing `oligotoxdb` pipeline
   (`oligotoxdb/oligotoxdb/ingestion.py`), replacing `ILLUSTRATIVE_EXAMPLE_not_real_data/`.
2. Rewrite `01_narrative_document.md` Sections 1, 2, and 4 with real numbers and findings
   (Sections 3, 5, 6 should only need light edits if actual execution deviated from the plan).
3. Spot-check `02_methodology_document.md` against what was actually done (vendor used, any
   dropped system, any changed concentration/timepoint).
4. Deposit the real dataset (Zenodo/GEO/PRIDE or equivalent) and fill in
   `dataset/README_raw_data_access.md`, or attach the real files directly.
5. Rebuild the PDFs (`tools/build_pdfs.py`) and re-verify page counts/format.
6. Complete the Challenge.gov registration form and export as PDF.
7. Use `email_draft.md`, with its pre-send checklist completed, to actually send the package.

## Where and how to submit

Per the Challenge Announcement ("How to Enter"): email the complete package to
**ncatsoligotox@mail.nih.gov** before the Phase 2 deadline (December 31, 2026, 11:59 PM ET).
There is no web upload portal — Challenge.gov hosts the listing and registration form only, not
submission intake. See `email_draft.md` for the ready-to-use template.
