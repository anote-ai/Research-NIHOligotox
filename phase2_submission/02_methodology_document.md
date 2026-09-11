# OligoTox Open Data Challenge — Phase 2 Methodology Document

**Submitter:** Anote, Inc. | **Point of Contact:** Natan Vidra
**Submission target:** Single PDF, ≤5 pages, 8.5"×11", ≥1" margins, ≥11pt Arial, line spacing ≥1.0

> **STATUS: DRAFT.** Materials and methods below describe the protocol as designed. Before
> submission, confirm each method was actually followed as written (or update to reflect any
> deviations made during real data generation) — the challenge requires this document to
> describe "the materials and methods used to generate the data contained in the dataset files,"
> i.e., what was actually done, not only what was planned.

## 1. Oligonucleotide Library

Library spans 7 backbone chemistries (PS, PO, PMO, LNA-mix, PNA, 2'-F-ANA, morpholino), 4
sugar modifications (DNA, 2'-OMe, 2'-F, LNA/mixed), and multiple structural classes (gapmers,
siRNA/siRNA-like duplexes, aptamer scaffolds, splice-switching oligos, GalNAc- and
lipid-conjugated ASOs), plus systematic CpG-content, GC-content, and length gradients, known
clinical/reference compounds, and scrambled/non-targeting negative controls.

## 2. Oligonucleotide Synthesis, Purification, and Identity Characterization

*(This section directly answers the challenge's required content: "methods used to purify and
characterize oligo identity.")*

- **Synthesis:** Solid-phase phosphoramidite synthesis via CRO partner(s) [name vendor(s) actually
  used, e.g., IDT / Eurofins Genomics].
- **Purification:** Reverse-phase HPLC (RP-HPLC) purification for all oligos; target purity
  ≥95% for standard chemistries, ≥90% for heavily modified chemistries.
- **Identity confirmation:** Mass spectrometry (LC-MS) confirmation of molecular weight against
  expected mass for each synthesized sequence/chemistry combination.
- **Purity re-assessment on outliers:** Any compound producing an unexpected toxicity signal is
  re-confirmed by MS before inclusion in the final dataset.
- **Endotoxin testing:** LAL or equivalent assay; acceptance threshold <1 EU/mL in working
  solution, to rule out endotoxin contamination as a confound in immunotoxicity/PBMC readouts.
- **Formulation:** Sterile filtration, lyophilization, and reconstitution in RNase-free PBS prior
  to dosing.
- **Reporting per compound:** vendor, synthesis scale, HPLC purity (%), MS-confirmed mass
  (pass/fail), endotoxin level (EU/mL) — all recorded in the Compound Registry table of the
  dataset (see `dataset/data_dictionary_and_schema.md`).

## 3. In Vitro Human Model Systems

*[List only the systems actually used to generate the submitted data; delete unused systems from
the planned five below or annotate as "not yet run."]*

1. **Primary human hepatocytes (PHH), sandwich culture** — hepatotoxicity panel.
2. **Human kidney proximal tubule organoids** — nephrotoxicity panel.
3. **Human PBMC co-culture** — immunotoxicity and complement activation panel.
4. **Human platelet-rich plasma / citrated plasma** — thrombocytopenia and coagulopathy panel.
5. **Liver-kidney microphysiological system (organ-on-chip)** — multi-organ / chronic exposure
   panel.

For each system used, report: cell/tissue source and donor count, dosing concentrations and
timepoints actually tested, and replicate structure (biological × technical).

## 4. Toxicology Readouts

Each system's readouts (viability, membrane integrity, apoptosis markers, organ-injury
biomarkers, cytokine/complement multiplex, platelet activation/aggregation, coagulation timing
assays, and — for subsets — transcriptomic/proteomic profiling) are enumerated in the endpoint
registry referenced in the dataset files. Report which of the planned endpoints were actually
measured for the submitted data.

## 5. Dosing and Replication

Report actual concentrations/doses tested per system, actual number of biological replicates
(donors/preparations) and technical replicates per condition, and any anchor/reference compounds
repeated across batches for drift monitoring.

## 6. Quality Control

- **Plate-level:** Z'-factor and vehicle-control CV computed per plate; plates failing threshold
  repeated.
- **Compound-level:** Outlier detection (e.g., Grubbs test) across replicates; dose-response
  curve fit quality (R²) flagged when below threshold.
- **Batch-level:** Anchor compounds repeated across batches; inter-batch CV monitored, with
  normalization applied if drift detected.
- Report actual QC pass rates and any compounds/plates excluded, with justification.

## 7. Data Processing and Curve Fitting

Raw instrument output normalized to vehicle control; dose-response curves fit via 4-parameter
logistic regression to derive IC50/EC50 (with 95% CI), Emax, Hill coefficient, NOEC, and
AUDRC. Software/package and version used for fitting should be named here for reproducibility.

---

*Prepared by Natan Vidra, Anote, Inc. Contact: oligotoxdb@anote.ai*
