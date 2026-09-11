# OligoTox Open Data Challenge — Phase 2 Narrative Document

**Submitter:** Anote, Inc. (Entity) | **Point of Contact:** Natan Vidra
**Submission target:** Single PDF, ≤12 pages, 8.5"×11", ≥1" margins, ≥11pt Arial, line spacing ≥1.0
**Challenge:** NIH NCATS Oligonucleotide Toxicity (OligoTox) Open Data Challenge — Phase 2 (Data Generation)

> **STATUS: DRAFT / PLACEHOLDER.** This document is structured to match the six required
> narrative sections exactly (challenge announcement, "Phase 2 Data Generation Phase," items
> 1–6). Sections 1–2 and 4–6 currently describe *planned* work because no wet-lab data has been
> generated yet. Per the requirement, Phase 2 is judged on "the demonstrated dataset as it is
> submitted" — this document must be rewritten in the past tense with real numbers, real
> distributions, and real findings before submission. Do not submit in its current form.

---

## 1. Executive Summary of the Dataset(s) Generated, and Positive/Negative Controls Included

*[PLACEHOLDER — fill in once Batch A–E data collection is complete or substantially underway.]*

Target content once real data exists:
- Final oligo count actually synthesized and screened (library design targets ~2,800; report actual N)
- Which toxicity domains were actually covered (hepatotoxicity, nephrotoxicity, immunotoxicity,
  complement activation, coagulopathy, thrombocytopenia — per the challenge's priority list)
- Which in vitro human systems were actually run (of the 5 proposed: PHH, kidney proximal
  tubule organoids, PBMC co-culture, platelet-rich plasma, liver-kidney MPS)
- The actual positive/negative controls included and their observed responses (see Appendix A
  of the design plan for the proposed control list — CpG ODN 2006, Poly-G ODN, ISIS 104838,
  Nusinersen, Mipomersen, scrambled PS-ASO, etc.) — report whether each behaved as expected

## 2. Summary of the Main Findings and Conclusions

*[PLACEHOLDER — requires real assay results.]* This section must report actual outcomes, e.g.:
- Which chemistries/motifs showed the strongest toxicity signal in which system
- Whether control compounds validated the assay panel (Z'-factor, CV, expected direction of effect)
- Any unexpected findings or negative results (both are valuable for an open dataset)
- Summary statistics (e.g., % of library flagged toxic per endpoint, dose-response quality metrics)

## 3. Description of How Data Were Produced

This section can largely be finalized now, since it describes process rather than results.

**Experimental design:** Tiered screening across five human in vitro systems (see
`02_methodology_document.md` for full protocol detail). Library of oligonucleotides spanning
multiple backbone chemistries (PS, PO, PMO, LNA-mix, PNA, 2'-F-ANA, morpholino), sugar
modifications (DNA, 2'-OMe, 2'-F, LNA, mixed), and structural classes (gapmers, siRNA,
aptamers, splice-switching oligos, GalNAc/lipid conjugates), plus systematic CpG-content,
GC-content, and length series, alongside known positive/negative controls.

**Data acquisition:** Automated liquid handling (acoustic dispensing for compound transfer,
programmable handlers for cell/reagent dispensing) with plate-reader-based endpoint capture
(luminescence, fluorescence, absorbance, multiplex immunoassay) and, for a subset, RNA-seq
and proteomics.

**Computational processing:** Raw plate reader output → normalization to vehicle control →
per-plate QC (Z'-factor ≥0.5, control CV ≤15%) → dose-response curve fitting (4-parameter
logistic) to derive IC50/EC50, Emax, Hill coefficient, NOEC, and AUDRC → batch-level drift
correction via repeated anchor compounds.

*[Once real data exists, replace the above generic description with actual batch counts, actual
QC pass/fail rates, and any deviations from the planned protocol.]*

## 4. Description of How Indicators and Predictor Variables Were Measured, Their Distribution, and the Distribution of Predictor Variables Amongst Tested Oligos

*[PLACEHOLDER — requires the actual assembled dataset.]*

Indicators (toxicity readouts) and predictor variables (sequence/chemistry features) are defined
in `dataset/data_dictionary_and_schema.md`. Once real data exists, this section must report:
- Distribution of each indicator across the tested library (e.g., histogram/summary stats of
  viability, cytokine levels, IC50 values, by assay system and endpoint)
- Distribution of predictor variables amongst tested oligos (e.g., how many oligos per backbone
  class, sugar modification, length bucket, GC-content bin, CpG count) — i.e., does the library
  actually achieve the chemical diversity it was designed for, or are there gaps/imbalances
- Missingness/coverage: which oligo × endpoint combinations have data vs. do not

## 5. Discussion of How the Results Address a Gap in the Publicly Available Data

Argument that can be drafted now, to be tightened once real scale/composition is known:
- Existing public oligo toxicity data is dominated by small (<200 compound), single-endpoint,
  animal-derived, PS-ASO-heavy datasets (see literature review in Phase 1 submission)
- This dataset's contribution: human in vitro origin, multi-endpoint (hepatotox, nephrotox,
  immunotox, complement, coagulopathy, thrombocytopenia in one coordinated panel), and
  chemically diverse (multiple backbone/sugar/conjugate classes) coverage
- *[Once real data exists: quantify the actual gap closed — e.g., "adds N chemically-diverse
  human-system data points, a Xx increase over the largest comparable public dataset."]*

## 6. Discussion of How the Data Could Be Used to Develop a Predictive Model

- Structured, multi-endpoint labels support both single-task and multi-task model training
- Chemistry + sequence features (documented in the data dictionary) support both interpretable
  feature-based models (e.g., gradient-boosted trees with SHAP) and sequence foundation-model
  fine-tuning
- Dose-response structure (not just binary tox/non-tox) supports regression and uncertainty-aware
  modeling, and benchmark splits (random, chemistry-stratified, prospective) support rigorous
  out-of-distribution evaluation
- *[Once real data exists: report actual baseline model performance on the real dataset — not
  simulated numbers — if any modeling has been done as part of the submission.]*

---

*Prepared by Natan Vidra, Anote, Inc. Contact: oligotoxdb@anote.ai*
