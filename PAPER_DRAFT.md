# OligoTox: Toward Multi-Modal Prediction of Oligonucleotide Toxicity
### Paper Draft Skeleton (early stage)

**Status:** Skeleton with abstract/intro/methods populated from `DESIGN_DOC.md` and the
Phase 2 NIH NCATS submission (`OligoTox_Phase2_Submission_Anote.md`). Results sections
contain only what has actually been computed in this repository so far; all
DESIGN_DOC.md target numbers are explicitly marked **(projected, pending full
experiment run)**. No fabricated measurements are included.

---

## Abstract

Oligonucleotide therapeutics (ASOs, siRNAs, aptamers) are a fast-growing drug class
whose clinical development is frequently limited by toxicity rather than efficacy.
Existing computational toxicity models are largely chemistry-only (e.g., ECFP4 +
QSAR) and ignore sequence-specific mechanisms such as off-target hybridization and
CpG-driven immune activation. We outline **OligoTox**, a research program to (1)
curate a multi-endpoint oligonucleotide toxicity benchmark ("OligoTox-Bench" /
"OligoToxDB"), (2) build a combined sequence+chemistry representation ("SeqChem"),
and (3) evaluate whether multi-task learning and off-target binding features improve
toxicity prediction over chemistry-only baselines, particularly under
scaffold/chemistry-split out-of-distribution evaluation. This draft documents the
current implementation status: a synthetic-data toolkit and rule-based/logistic
baseline are implemented and runnable end-to-end; the curated real-world dataset,
SeqChem representation, multi-task model, and off-target integration described in
the design document are not yet implemented. **(projected, pending full experiment
run)** is used throughout to flag target numbers that have not been measured.

## 1. Introduction

Oligonucleotide drugs (15+ approved, hundreds in clinical trials) commonly fail late
in development due to hepatotoxicity, nephrotoxicity, immunostimulation, complement
activation, thrombocytopenia, and injection-site reactions (see `DESIGN_DOC.md`,
"Problem Statement & Novelty"). Four gaps motivate this work:

1. Chemistry-only QSAR models miss sequence-specific mechanisms (off-target
   hybridization, CpG/TLR9 activation).
2. Toxicity data is sparse and siloed across in vitro, animal, and clinical sources.
3. Toxicity endpoints are modeled independently despite shared biological mechanisms.
4. Public benchmarks suffer from scaffold/chemistry leakage between train and test
   sets, inflating reported performance.

This repository's companion document, `OligoTox_Phase2_Submission_Anote.md`, describes
a parallel, complementary effort: a wet-lab data-generation program (OligoToxDB,
2,800 oligos x 47 endpoints across 5 human in vitro systems) submitted to the NIH
NCATS OligoTox Open Data Challenge. That program, if funded and executed, would
directly produce the kind of curated, multi-endpoint dataset that the AI/ML
objectives in `DESIGN_DOC.md` require. As of this draft, that wet-lab data has not
yet been generated (data release schedule begins August 2026 per the submission
document), so all AI/ML results in this paper draft use synthetic data only.

## 2. Related Work

- DeepTox, AttentiveFP, Tox21 (chemistry-only / graph-based toxicity models)
- Nucleotide Transformer (Dalla-Torre et al. 2023) -- sequence foundation model,
  proposed as a fine-tuning backbone in `OligoTox_Phase2_Submission_Anote.md` Sec 8.2
- Sewing et al. (2016) -- hepatotoxicity prediction from sequence + modification
  pattern, closest prior work to the SeqChem hypothesis
- MoleculeNet / ADMET-AI benchmark conventions, referenced as the target format for
  the proposed OligoTox-RandBench / ChemSplit / ProspBench splits

*(Full literature review: future work -- see Improvement Plan.)*

## 3. Methods

### 3.1 Data

- **Implemented today:** `src/oligotox/data.py::make_dataset()` generates synthetic
  oligonucleotides with three programmatically-defined risk archetypes (PS-backbone
  hepatotoxicity-elevated, high-GC complement-risk, low-risk PO-backbone), with
  labels drawn directly from the generating rule plus noise.
- **Designed, not yet implemented:** the 3,200-point curated real dataset spanning
  NIH toxicity DB, ChEMBL, ClinicalTrials.gov, and internal industry data
  (`DESIGN_DOC.md`, "Dataset Construction"), and the larger 2,800-oligo x 47-endpoint
  OligoToxDB wet-lab dataset (`OligoTox_Phase2_Submission_Anote.md` Sec 3-6), which is
  on a rolling release schedule starting August 2026.

### 3.2 Features

`src/oligotox/core.py::extract_features()` computes 7 features per oligonucleotide:
GC content, CpG ratio, sequence length, is-PS-backbone flag, is-LNA-backbone flag,
modification count, and a hand-coded motif risk score (CCGG/GGGG/TTTT/GCGC/CGCG/AAAA
weighted k-mer matches).

This is a small fraction of the **SeqChem representation** specified in
`DESIGN_DOC.md`: 64-dim k-mer embeddings, top-50 predicted off-target binding
features, 2048-dim ECFP4 fingerprints, and 32-dim modification encoding. None of the
k-mer embedding, off-target binding predictor (TargetScan/RNAhybrid), or ECFP4/SMILES
pipeline exists in this repository yet.

### 3.3 Models

- **Implemented:** `OligotoxPipeline` (logistic regression over the 7 features above);
  a rule-based `ToxicityPredictor` baseline.
- **Designed, not yet implemented:** QSAR (Random Forest + ECFP4), DeepTox,
  AttentiveFP, SeqOnly, SeqChem v1/v2/v3 (multi-task, off-target) from
  `DESIGN_DOC.md`'s "Models Under Evaluation" table. The `oligotoxdb/models/`
  directory contains scaffolding for `OligoTox-XGB`, `OligoTox-Transformer`
  (Nucleotide Transformer fine-tune), and `OligoTox-ActiveLearn` (Gaussian Process
  active learning), aimed at the wet-lab dataset rather than the SeqChem benchmark
  described in `DESIGN_DOC.md`.

### 3.4 Evaluation

`src/oligotox/evaluate.py::aucroc()` implements AUROC from first principles
(trapezoidal rule). `experiments/exp0_baseline.py` (added alongside this draft) runs
the implemented logistic-regression baseline on synthetic data with a random 80/20
split over multiple seeds and reports the *measured* mean AUROC. Murcko scaffold
splitting (`DESIGN_DOC.md` Experiment 1) is not implemented; oligonucleotides do not
have a conventional Murcko scaffold the way small molecules do, so a domain-specific
splitting scheme (e.g., by backbone class, sequence-identity clustering, or batch) is
needed -- this is the most pressing methodological gap for OOD evaluation.

## 4. Results

| Experiment | Status | Result |
|---|---|---|
| Baseline (Exp. 0): logistic regression on synthetic data, random split | **Measured** (`experiments/exp0_baseline.py`, run at draft time) | Mean AUROC ~ 0.99 over 10 seeds, n=400 synthetic oligos/seed. **Caveat: labels are deterministic functions of the same features used for prediction, so this measures pipeline correctness, not real-world predictive validity.** |
| Baseline (Exp. 0): QSAR (RF + ECFP4) on real hepatotoxicity data | Not implemented | (projected, pending full experiment run) -- DESIGN_DOC.md target: AUROC ~ 0.79 |
| Exp. 1: Scaffold-split OOD generalization | Not implemented (no scaffold-split code, no real dataset) | (projected, pending full experiment run) -- DESIGN_DOC.md targets: -17pp (QSAR) to -8pp (SeqChem-OT) OOD drop |
| Exp. 2: SeqChem vs. chemistry-only | Not implemented (no SeqChem representation) | (projected, pending full experiment run) -- DESIGN_DOC.md target: +17pp AUROC for SeqChem over chemistry-only |
| Exp. 3: Multi-task learning across 6 endpoints | Not implemented (no multi-task model; only 5 endpoint enum values exist in code vs. 6 in design doc) | (projected, pending full experiment run) -- DESIGN_DOC.md target: +3 to +9pp per endpoint |
| Exp. 4: Off-target binding integration | Not implemented (no off-target predictor) | (projected, pending full experiment run) -- DESIGN_DOC.md target: +5pp AUROC at top-50 off-target sites |

No claim in this section beyond the single "Measured" row should be treated as an
experimental result. We explicitly avoid restating DESIGN_DOC.md's "Expected results"
tables as if they were obtained -- they are hypotheses to be tested once the
corresponding code and data exist.

## 5. Discussion

The core scientific hypothesis -- that sequence-aware (SeqChem) representations
outperform chemistry-only QSAR models on out-of-distribution oligonucleotide
toxicity prediction, especially via captured off-target hybridization effects -- is
plausible given known biology (Sewing et al. 2016; CpG/TLR9 mechanisms) but remains
untested in this codebase. The repository currently has two largely disconnected
efforts: (a) a wet-lab data-generation grant submission (OligoToxDB) aimed at NIH
NCATS, and (b) a lightweight synthetic ML toolkit (`src/oligotox`) that does not yet
consume real data or implement the SeqChem/multi-task/off-target methods described
in `DESIGN_DOC.md`. Bridging these -- i.e., actually building SeqChem, scaffold-style
splits, and multi-task training, and running them against either OligoToxDB (once
released) or an interim curated public dataset -- is the critical path to obtaining
real results for Experiments 1-4.

## 6. Limitations

- No real (non-synthetic) toxicity data is used anywhere in current results.
- Synthetic data labels are generated from the same feature space used for
  prediction, which inflates AUROC and limits the informativeness of any synthetic
  benchmark number.
- No scaffold/chemistry-based OOD split exists.
- No SeqChem, multi-task, or off-target model is implemented.
- Sample sizes, in any future real-data run, will be small (per `DESIGN_DOC.md`,
  thrombocytopenia n=400, positive rate 12%), so statistical power for endpoint-level
  claims will need explicit confidence intervals, not point estimates.

## 7. Future Work (see also Improvement Plan in the audit issue)

1. Implement Murcko-style or backbone/sequence-cluster-based scaffold splitting.
2. Implement a real ECFP4/RDKit chemistry-only baseline (Experiment 0/2 chemistry arm).
3. Implement the SeqChem k-mer + chemistry concatenated representation.
4. Implement a shared-encoder multi-task model (Experiment 3).
5. Wire in TargetScan/RNAhybrid-style off-target scoring (Experiment 4).
6. Replace synthetic data with the first real OligoToxDB batch once released
   (per `OligoTox_Phase2_Submission_Anote.md` Sec 6, Batch A target August 2026), or
   with a public interim dataset (e.g., literature-curated hepatotoxic ASO panel)
   in the meantime.

---

*This draft was scaffolded as part of a research-readiness audit. Sections 1-3 and
the limitations/future-work sections are populated from existing project documents;
Section 4 contains exactly one measured result, computed by running
`experiments/exp0_baseline.py` in this repository. All other numbers are explicitly
marked as projections from `DESIGN_DOC.md`, not measurements.*
