# OligoTox Open Data Challenge — Phase 2 Methodology Document

**Submitter:** Anote, Inc. | **Point of Contact:** Natan Vidra
**Submission target:** Single PDF, ≤5 pages, 8.5"×11", ≥1" margins, ≥11pt Arial, line spacing ≥1.0

## 1. Oligonucleotide Library

The library spans 7 backbone chemistries (PS, PO, PMO, LNA-mix, PNA, 2'-F-ANA, morpholino), 4
sugar modifications (DNA, 2'-OMe, 2'-F, LNA/mixed), and multiple structural classes (gapmers,
siRNA/siRNA-like duplexes, aptamer scaffolds, splice-switching oligos, GalNAc- and
lipid-conjugated ASOs), plus systematic CpG-content, GC-content, and length gradients, known
clinical/reference compounds, and scrambled/non-targeting negative controls.

## 2. Oligonucleotide Synthesis, Purification, and Identity Characterization

*(This section directly answers the challenge's required content: "methods used to purify and
characterize oligo identity.")*

- **Synthesis:** Solid-phase phosphoramidite synthesis via Integrated DNA Technologies (IDT) and
  Eurofins Genomics, split by chemistry class (standard PS/PO chemistries via either vendor;
  heavily modified chemistries — LNA, PMO, PNA, 2'-F-ANA — routed to whichever vendor's
  catalog covers that modification).
- **Purification:** Reverse-phase HPLC (RP-HPLC) purification for all oligos; target purity
  ≥95% for standard chemistries, ≥90% for heavily modified chemistries.
- **Identity confirmation:** Mass spectrometry (LC-MS) confirmation of molecular weight against
  expected mass for each synthesized sequence/chemistry combination.
- **Purity re-assessment on outliers:** Any compound producing an unexpected toxicity signal is
  re-confirmed by MS before inclusion in the final dataset.
- **Endotoxin testing:** LAL (Limulus Amebocyte Lysate) assay; acceptance threshold <1 EU/mL in
  working solution, to rule out endotoxin contamination as a confound in immunotoxicity/PBMC
  readouts.
- **Formulation:** Sterile filtration, lyophilization, and reconstitution in RNase-free PBS prior
  to dosing.
- **Reporting per compound:** vendor, synthesis scale, HPLC purity (%), MS-confirmed mass
  (pass/fail), endotoxin level (EU/mL) — all recorded in the Compound Registry table of the
  dataset (see `dataset/data_dictionary_and_schema.md`).

## 3. In Vitro Human Model Systems

| System | Source | Dosing | Timepoints |
|---|---|---|---|
| Primary human hepatocytes (PHH), sandwich culture | Commercial human hepatocyte supplier, 3 donors/batch (n=3 biological replicates) | 6 concentrations: 0.01, 0.1, 1, 10, 50, 100 µM | 24h and 72h |
| Human kidney proximal tubule organoids | Certified human kidney organoids, differentiated to proximal tubule identity | 5 concentrations: 0.1, 1, 10, 50, 100 µM | 48h and 96h |
| Human PBMC co-culture | Fresh PBMCs, 5 healthy donors (3 per run, rotating) | 4 concentrations: 0.1, 1, 10, 50 µM | 24h |
| Human platelet-rich plasma / citrated plasma | Fresh citrated plasma, 5 donors; platelets isolated by differential centrifugation | 3 concentrations: 1, 10, 50 µM | 30 min and 2h |

*(If the final protocol drops a system or changes a parameter above, update this table to match
what was actually run — the dataset's own QC records, not this document, are the source of
truth for what was actually executed.)*

## 4. Toxicology Readouts

- **PHH (hepatotoxicity):** cell viability (ATP-lite), LDH release, ALT/AST secretion, caspase
  3/7 activity, HMGB1 release, ROS (CellROX), lipid accumulation (Oil Red O), mitochondrial
  membrane potential (JC-1), bile transport (CDFDA efflux); RNA-seq for a high-priority subset.
- **Kidney organoids (nephrotoxicity):** organoid viability (CellTiter-Glo 3D), KIM-1 and NGAL
  secretion, caspase 3/7, lysosomal accumulation (LysoTracker), tight-junction integrity
  (ZO-1 imaging), mitochondrial stress (MitoTracker).
- **PBMC (immunotoxicity, complement):** cytokine multiplex (IFN-α, IFN-γ, IL-1β, IL-6, IL-8,
  IL-10, IL-12p70, TNF-α, MCP-1, IP-10, GM-CSF, IL-4/5/13, IL-17A), complement C3a/C5a, TLR9
  pathway activation (NF-κB reporter), NK-cell activation (CD69).
- **Platelet/plasma (thrombocytopenia, coagulopathy):** platelet aggregation (light transmission
  aggregometry), P-selectin (CD62P) and CD63 activation markers, aPTT and PT clotting times,
  Factor Xa and thrombin inhibition (chromogenic assays).

The full endpoint registry (all readouts, units, and controlled-vocabulary names) is implemented
in `oligotoxdb/oligotoxdb/endpoints.py` and documented in
`dataset/data_dictionary_and_schema.md`.

## 5. Dosing and Replication

- **Biological replicates:** n=3 independent donors/preparations for PHH, PBMC, and kidney
  organoid assays; n=5 donors for platelet/plasma work (rotated across runs).
- **Technical replicates:** n=2 per biological replicate for all primary endpoints.
- **Anchor compounds:** a fixed panel of reference compounds is repeated in every batch to
  monitor and correct for inter-batch drift.
- Concentrations and timepoints per system are fixed in Section 3 above and are not adjusted
  compound-by-compound, so results are directly comparable across the library.

## 6. Quality Control

- **Plate-level:** Z'-factor (target ≥0.5) and vehicle-control coefficient of variation (target
  ≤15%) computed from positive/negative control wells on every plate; plates failing either
  threshold are flagged and repeated, not silently included.
- **Compound-level:** outlier detection across biological replicates (Grubbs' test, α=0.05);
  dose-response curve fit quality (R² ≥0.85 target) flagged when below threshold.
- **Batch-level:** anchor compounds (Section 5) monitored for inter-batch CV; drift beyond
  threshold triggers normalization (e.g., ComBat) before release.
- Final QC pass/fail counts, and any compounds or plates excluded with justification, are
  reported in the dataset's own QC table (`plate_qc.csv` in the released data), not duplicated
  here.

## 7. Data Processing and Curve Fitting

Raw instrument output is normalized to each plate's vehicle control; dose-response curves are
fit via 4-parameter logistic regression (`oligotoxdb/oligotoxdb/qc.py`, using `scipy.optimize`)
to derive IC50/EC50 with 95% confidence interval, Emax, Hill coefficient, NOEC, and AUDRC. The
fitting code, exact package versions, and QC thresholds are version-controlled in this
repository so the derivation is fully reproducible from raw values.

---

*Prepared by Natan Vidra, Anote, Inc. Contact: nvidra@anote.ai*
