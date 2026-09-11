# OligoTox Open Data Challenge — Phase 2 Narrative Document

**Submitter:** Anote, Inc. (Entity) | **Point of Contact:** Natan Vidra
**Submission target:** Single PDF, ≤12 pages, 8.5"×11", ≥1" margins, ≥11pt Arial, line spacing ≥1.0
**Challenge:** NIH NCATS Oligonucleotide Toxicity (OligoTox) Open Data Challenge — Phase 2 (Data Generation)

> **STATUS: DRAFT.** This document is structured to match the six required narrative sections
> exactly. Sections 1, 2, and 4 are placeholders — they require the actual generated dataset,
> which does not exist yet, and must not be filled with invented numbers. Sections 3, 5, and 6
> describe the experimental design, the gap this work addresses, and the modeling plan — none of
> that content depends on results, so it is written in full below and should only need light
> editing once real data collection begins (e.g., to reflect deviations from the plan). **Do not
> submit until Sections 1, 2, and 4 are rewritten from real experimental outcomes.**

---

## 1. Executive Summary of the Dataset(s) Generated, and Positive/Negative Controls Included

*[PLACEHOLDER — requires the completed dataset.]* When filled in, this section reports:
- The final oligo count actually synthesized and screened, broken out by chemistry class
- Which toxicity domains from the priority list were actually covered (hepatotoxicity,
  nephrotoxicity, immunotoxicity, complement activation, coagulopathy, thrombocytopenia)
- Which of the in vitro human systems described in Section 3 were actually run
- The positive/negative controls actually included (see the candidate list in
  `02_methodology_document.md` §5) and a one-line statement of whether each behaved as expected

## 2. Summary of the Main Findings and Conclusions

*[PLACEHOLDER — requires real assay results.]* When filled in, this section reports actual
outcomes: which chemistries or sequence motifs produced the strongest toxicity signal in which
system; whether control compounds validated the assay panel (Z'-factor, CV, direction of effect
matching literature expectation); any unexpected or negative findings (both are informative for
an open dataset); and summary statistics such as the fraction of the library flagged toxic per
endpoint and the proportion of dose-response curves meeting the fit-quality threshold.

## 3. Description of How Data Were Produced

### 3.1 Oligonucleotide Library

The library is designed to give the resulting dataset broad, systematic coverage of the
chemical variables believed to drive oligo toxicity, rather than being a random or
convenience sample. It spans:

- **Backbone chemistries:** phosphorothioate (PS), phosphodiester (PO), phosphorodiamidate
  morpholino (PMO), LNA-mixed backbones, peptide nucleic acid (PNA), and 2'-F-ANA
- **Sugar modifications:** unmodified DNA, 2'-O-methyl (2'-OMe), 2'-fluoro (2'-F), locked
  nucleic acid (LNA), and mixed-modification patterns
- **Structural classes:** gapmers (with systematic wing/gap length variation), siRNA/siRNA-like
  duplexes, aptamer scaffolds, splice-switching oligonucleotides, and GalNAc- or
  lipid-conjugated ASOs
- **Systematic gradients:** CpG-dinucleotide content, overall GC content, and sequence length,
  each varied independently of chemistry so their individual contribution to toxicity signal can
  be isolated statistically rather than confounded with backbone choice
- **Controls:** clinically approved/advanced reference ASOs with known safety profiles (as
  benchmarks), known toxicity-positive compounds (e.g., a CpG ODN TLR9 agonist, a
  hepatotoxic ASO analog), and scrambled/non-targeting sequences matched in chemistry (as
  chemistry-matched negative controls) — full rationale for each in
  `02_methodology_document.md`

This design choice — varying one axis (chemistry, CpG content, GC content, length) at a time
against otherwise-matched backgrounds — is what lets the resulting dataset support attribution
(which variable actually drove an observed toxicity signal), not just prediction.

### 3.2 In Vitro Human Model Systems

Toxicity is assessed using human-derived in vitro systems chosen for direct relevance to the
priority toxicity domains, rather than substitute animal models:

| System | Toxicity domain | Why this system |
|---|---|---|
| Primary human hepatocytes (PHH), sandwich culture | Hepatotoxicity | Retains CYP activity, transport protein expression, and hepatocyte polarity relevant to oligo uptake/metabolism — the field's accepted gold standard in vitro hepatotoxicity model |
| Human kidney proximal tubule organoids | Nephrotoxicity | Proximal tubule is the primary site of oligo accumulation via scavenger-receptor uptake and the dominant nephrotoxicity signal in clinical ASO experience; organoids better capture tubular transporter expression than 2D culture |
| Human PBMC co-culture | Immunotoxicity, complement activation | Captures innate immune activation (monocyte, NK, B-cell) most relevant to CpG- and PS-backbone-driven immunostimulation |
| Human platelet-rich plasma / citrated plasma | Thrombocytopenia, coagulopathy | Directly measures platelet activation/aggregation and coagulation factor interaction, the proposed mechanism behind PS-backbone-associated thrombocytopenia |

Each system's dosing range, timepoints, and full readout panel are specified in
`02_methodology_document.md`. Compounds are prioritized across systems using a tiered strategy
(a full panel for a subset of high-information compounds, a narrower panel for the rest), so
that the completed dataset trades comprehensive multi-endpoint data on some compounds against
broader single-endpoint coverage across the whole library — both are useful signal for model
training, and the schema (`dataset/data_dictionary_and_schema.md`) records exactly which
combinations were tested and which were not.

### 3.3 Data Acquisition and Computational Processing

Data acquisition uses automated liquid handling for compound dosing and plate-reader-based
endpoint capture (luminescence, fluorescence, absorbance, multiplex immunoassay), supplemented
by transcriptomic and/or proteomic profiling for a subset of high-priority compounds.

Raw instrument output is processed as follows before entering the released dataset:
1. Normalization of each raw reading to its plate's vehicle control
2. Per-plate quality control: Z'-factor (target ≥0.5) and vehicle-control coefficient of
   variation (target ≤15%) computed from positive/negative control wells; plates failing either
   threshold are flagged and repeated rather than silently included
3. Outlier detection across biological replicates (Grubbs' test)
4. Dose-response curve fitting (4-parameter logistic) on QC-passed, outlier-cleaned data, to
   derive IC50/EC50 with 95% confidence interval, Emax, Hill coefficient, no-observed-effect
   concentration (NOEC), and area-under-the-dose-response-curve (AUDRC) as a single summary
   statistic
5. Batch-level drift correction using repeated anchor compounds present in every batch

Every stage above — raw value, normalized value, QC flag, and derived statistic — is retained in
the released dataset (see `dataset/data_dictionary_and_schema.md`, Tables 3–4), not just the
final fitted values, so downstream users can re-derive or re-QC the data under their own criteria
rather than trusting ours by default.

*[Once real data exists: replace the generic description above with the actual number of
batches run, actual QC pass/fail rates, and any deviations from this plan — e.g., a system that
was dropped, a timepoint that was changed.]*

## 4. Description of How Indicators and Predictor Variables Were Measured, Their Distribution, and the Distribution of Predictor Variables Amongst Tested Oligos

*[PLACEHOLDER — requires the actual assembled dataset.]* The full definitions of every indicator
(toxicity readout) and predictor variable (sequence/chemistry feature) are fixed in advance in
`dataset/data_dictionary_and_schema.md` — Table 2 (computed sequence/chemistry features) and
Table 3 (experimental readouts). Once real data exists, this section must report: the
distribution of each indicator across the tested library (e.g., what fraction of compounds show
a dose-dependent hepatotoxicity signal, the spread of IC50 values by endpoint); the distribution
of predictor variables actually achieved (how many tested oligos fall into each backbone class,
sugar modification, length bucket, GC-content bin, and CpG-count bin — i.e., whether the library
actually achieved the systematic chemical diversity described in Section 3.1, or whether gaps
emerged during real synthesis/screening); and the coverage matrix of which oligo × assay-system ×
endpoint combinations have data versus were not tested (per the tiered strategy in Section 3.2).

## 5. Discussion of How the Results Address a Gap in the Publicly Available Data

Public data available today for training oligo toxicity models has five well-documented
limitations, each of which this dataset is designed to directly address:

1. **Small scale.** Most published oligo toxicology datasets contain well under 200 compounds
   with a single toxicity readout each, too little signal for multi-endpoint or deep-learning
   approaches. This dataset's systematic library design (Section 3.1) is sized specifically to
   exceed that scale within a single coordinated study.
2. **Chemical homogeneity.** Public data heavily over-represents phosphorothioate-backbone ASOs;
   gapmer, LNA, 2'-OMe, and PMO chemistries are comparatively underrepresented, which biases any
   model trained only on public data toward PS-backbone chemistry. The library's backbone and
   sugar-modification axes (Section 3.1) are explicitly designed to correct this imbalance.
3. **Animal origin.** The majority of available toxicology data derives from rodent in vivo
   studies, which are known to mis-predict human hepatotoxicity and immunotoxicity due to
   species differences in innate immune receptor profiles and liver biology — the exact
   motivation stated in NCATS's RFI on this challenge. Every system used here (Section 3.2) is
   human-derived.
4. **Proprietary lock-up.** The largest existing oligo toxicity datasets are held by industry
   and not publicly accessible, so the research community cannot use them to build or validate
   models regardless of their scale or quality. This dataset is released under CC BY 4.0 with
   no access restriction (see `03_public_access_and_dissemination_plan.md`).
5. **Siloed, single-endpoint measurement.** Existing datasets rarely combine orthogonal
   mechanistic readouts (e.g., transcriptomic + phenotypic) for the same compound, which
   prevents multi-task modeling that could exploit shared biology across toxicity domains. This
   dataset's multi-endpoint panel per compound (Section 3.2, tiered strategy) is designed
   specifically to enable that.

*[Once real data exists: replace the general argument above with the actual, quantified gap
closed — e.g., "this dataset adds N real human-in-vitro data points across M chemistry classes,
a Kx increase over the largest comparable public dataset in [specific toxicity domain]," citing
the specific public datasets being compared against.]*

## 6. Discussion of How the Data Could Be Used to Develop a Predictive Model

The dataset is structured, from the schema stage onward, to be directly usable for training and
evaluating predictive models of oligo toxicity, in three complementary ways:

**Feature-based models.** The computed sequence/chemistry features in
`dataset/data_dictionary_and_schema.md` Table 2 (GC content, CpG count, predicted secondary
structure, predicted melting temperature, off-target hybridization score, backbone/sugar
one-hot encodings, physicochemical properties) support interpretable models such as
gradient-boosted trees, where per-prediction feature attribution (e.g., SHAP values) lets a
researcher see which molecular features drove a given toxicity prediction — useful for
generating testable mechanistic hypotheses, not just point predictions.

**Sequence-based / foundation-model approaches.** Because the raw nucleotide sequence and
per-position chemical modification locations are retained (not just summary chemistry labels),
the dataset also supports fine-tuning of pretrained genomic sequence models, which can capture
sequence-motif effects (e.g., specific CpG or poly-G contexts) that hand-engineered features may
miss.

**Multi-task and dose-response-aware modeling.** Because each compound is measured across
multiple toxicity endpoints where the tiered protocol allows (Section 3.2) and full
dose-response curves are retained rather than collapsed to a single binary toxic/non-toxic
label (Section 3.3), the dataset supports multi-task learning that can share statistical
strength across correlated endpoints, and regression- or uncertainty-aware modeling of the full
concentration-response relationship rather than a single threshold call.

**Rigorous evaluation.** The dataset design also enables benchmark splits harder than a random
train/test split: a chemistry-stratified split (train on some backbone/sugar combinations, test
on held-out ones, to measure whether a model generalizes across chemistry rather than
memorizing it) and — since data collection proceeds in batches (see
`02_methodology_document.md`) — a prospective split (train on early batches, test on later
ones), which is a more honest test of real-world predictive utility than a random split alone.

*[Once real data exists: report actual baseline model performance on the real dataset — not
simulated numbers — if any modeling has been done as part of the submission, and note which of
the above approaches were actually tried.]*

---

*Prepared by Natan Vidra, Anote, Inc. Contact: nvidra@anote.ai*
