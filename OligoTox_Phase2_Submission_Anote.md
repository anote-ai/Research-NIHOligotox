# OligoTox Open Data Challenge – Phase 2 Submission
## Anote, Inc. | Data Generation & Public Dissemination Plan

**Lead Innovator:** Natan Vidra, Anote, Inc.  
**Submission Phase:** Phase 2 – Data Generation  
**Submission Date:** June 2026  
**Challenge Sponsor:** NIH National Center for Advancing Translational Sciences (NCATS)  

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Scientific Background & Rationale](#2-scientific-background--rationale)
3. [Oligonucleotide Library Design](#3-oligonucleotide-library-design)
4. [Experimental Systems & Model Selection](#4-experimental-systems--model-selection)
5. [Toxicology Assay Panel](#5-toxicology-assay-panel)
6. [Data Generation Plan](#6-data-generation-plan)
7. [Data Schema & Annotation Standards](#7-data-schema--annotation-standards)
8. [AI/ML Framework: OligoTox-Predict](#8-aiml-framework-oligotox-predict)
9. [Quality Assurance & Reproducibility](#9-quality-assurance--reproducibility)
10. [Public Access & Dissemination Plan (PADP)](#10-public-access--dissemination-plan-padp)
11. [Team & Collaborators](#11-team--collaborators)
12. [Budget Overview](#12-budget-overview)
13. [Project Timeline](#13-project-timeline)
14. [Anticipated Impact](#14-anticipated-impact)

---

## 1. Executive Summary

Anote, Inc. proposes to generate and publicly release **OligoToxDB** — the largest, most comprehensively annotated open dataset of human in vitro oligonucleotide toxicity profiles ever assembled. Building on our Phase 1 ideation framework, we will experimentally screen **2,800 structurally diverse oligonucleotides** across **six clinically critical toxicity domains** (hepatotoxicity, nephrotoxicity, immunotoxicity, complement activation, coagulopathy, and thrombocytopenia) using **five human-relevant in vitro model systems**. Each oligo will be characterized by a rich multimodal feature profile: sequence, chemical modifications, secondary structure, predicted hybridization thermodynamics, and up to **47 experimental toxicity readouts** across dose-response conditions.

Uniquely, we integrate an **AI-in-the-loop experimental design strategy** — our OligoTox-Predict machine learning system actively guides which compounds to screen next, maximizing information gain per experiment and accelerating discovery while reducing cost. The resulting dataset, associated ML models, analysis tools, and interactive data portal will all be released under permissive open licenses, immediately available to the global research community.

This work directly addresses a critical gap: the oligonucleotide therapeutics field lacks the large-scale, standardized, human-relevant toxicity datasets needed to train predictive AI/ML models, currently forcing reliance on costly, often poorly predictive animal studies. OligoToxDB will change that.

**Total oligos screened:** 2,800  
**Toxicity endpoints:** 47 readouts across 6 organ/mechanism domains  
**Model systems:** 5 human in vitro platforms  
**Data release:** Fully open, FAIR-compliant, by November 30, 2026  
**ML models released:** 3 model architectures, all weights and training code open-sourced  

---

## 2. Scientific Background & Rationale

### 2.1 The Oligonucleotide Toxicity Problem

Oligonucleotide (oligo) therapeutics — antisense oligonucleotides (ASOs), siRNA, miRNA mimics/inhibitors, aptamers, and splice-switching oligonucleotides — represent one of the fastest-growing drug modalities, with over 15 approved drugs and hundreds in clinical trials. Yet toxicity remains the primary cause of clinical attrition. Key toxicities include:

- **Hepatotoxicity**: The liver is the primary organ of accumulation for many backbone chemistries. Phosphorothioate (PS) modifications, high GC content, and specific sequence motifs (CpG, poly-G) cause hepatocellular injury through both sequence-dependent (hybridization) and sequence-independent (protein binding) mechanisms.
- **Nephrotoxicity**: Kidney proximal tubule cells actively uptake and accumulate oligos through scavenger receptor pathways, leading to tubular injury — the dominant safety signal in early clinical ASO trials.
- **Immunotoxicity & Complement Activation**: Unmethylated CpG motifs activate TLR9; backbone chemistry activates innate immune pathways; certain sequence patterns trigger complement cascade through both classical and alternative pathways.
- **Coagulopathy & Thrombocytopenia**: PS-backbone oligos bind coagulation factors and platelet surface proteins, extending clotting times and reducing platelet counts — an underappreciated and incompletely mechanistically understood toxicity.

### 2.2 Why Existing Data Is Insufficient

Current public data on oligo toxicity suffers from five critical deficiencies:

1. **Small scale**: Most published datasets contain <200 compounds with single-endpoint readouts.
2. **Chemical homogeneity**: Datasets heavily over-represent PS-backbone ASOs; gapmer, LNA, 2'-OMe, and phosphorodiamidate morpholino (PMO) chemistries are underrepresented.
3. **Animal-based**: The majority of available toxicology data comes from rodent in vivo studies, which systematically mis-predict human hepatotoxicity and immunotoxicity due to species differences in innate immune receptor profiles, complement regulation, and liver biology.
4. **Proprietary**: Industry-generated toxicity data, typically the most extensive, is not publicly accessible.
5. **Siloed**: No dataset integrates orthogonal mechanistic readouts (transcriptomics + proteomics + phenotypic endpoints) that would enable multi-task AI/ML models.

### 2.3 Our Approach Fills This Gap

OligoToxDB addresses each deficiency: large-scale (2,800 compounds), chemically diverse (7 backbone classes, 4 sugar modifications, 3 conjugation types), fully human in vitro, publicly released, and multimodal. The resulting dataset is purpose-built for AI/ML training — structured, annotated, and accompanied by curated ML benchmarks.

---

## 3. Oligonucleotide Library Design

### 3.1 Library Composition (n = 2,800)

Our library is designed to maximize structural and chemical diversity while ensuring systematic coverage of known and hypothesized toxicity drivers. Composition:

| Category | Count | Description |
|---|---|---|
| **Clinically approved or advanced ASOs** | 120 | Reference standards; known safety profiles |
| **Systematic backbone variation** | 350 | Identical sequences across 7 backbone chemistries |
| **Systematic sugar modification** | 280 | Fixed backbone; 4 sugar modifications (2'-OMe, 2'-F, LNA, DNA) |
| **CpG motif series** | 150 | Graded CpG content; methylation variants |
| **Poly-G / G-quadruplex series** | 100 | Known complement and coagulation activators |
| **GC content gradient** | 200 | 20% → 80% GC; length-matched |
| **Length series** | 120 | 8-mer → 25-mer; matched chemistry |
| **Gapmer design space** | 300 | Systematic gap/wing length variation |
| **siRNA / siRNA-like** | 200 | dsRNA oligos; RISC-mediated targeting variants |
| **Aptamer scaffolds** | 100 | Thrombin, PDGF, and randomized aptamers |
| **Splice-switching oligonucleotides** | 150 | Exon-skipping designs |
| **GalNAc-conjugated ASOs** | 200 | Liver-targeted; conjugate variation |
| **Lipid-conjugated oligos** | 100 | Membrane affinity series |
| **Scrambled / negative controls** | 150 | Matched chemistry; non-targeting sequences |
| **Known toxic positive controls** | 80 | CpG ODN 2006, ISIS 104838 hepatotoxic variants |
| **Active learning-selected** | 200 | Iteratively selected by OligoTox-Predict (see §8) |
| **Novel/exploratory** | 300 | Novel chemical space candidates |

**Total: 2,800 oligonucleotides**

### 3.2 Chemical Metadata Captured Per Compound

Each oligo in OligoToxDB is annotated with the following features, which serve as both experimental metadata and AI/ML model input features:

**Sequence & Structure Features:**
- Full nucleotide sequence (5'→3')
- Length (nt)
- GC content (%)
- CpG dinucleotide count and position
- Poly-purine / poly-pyrimidine run length
- Predicted secondary structure (minimum free energy, MFE) via RNAfold/Vienna RNA
- Predicted melting temperature (Tm) against matched target and common off-targets
- G-quadruplex propensity score (QGRS Mapper)
- Self-complementarity index
- BLAST off-target binding score (human transcriptome)
- Predicted hybridization ΔG to intended target

**Chemical Modification Features:**
- Backbone chemistry (PS, PO, PMO, PNA, LNA-mix, 2'-F-ANA, morpholino)
- Sugar modification pattern (DNA, 2'-OMe, 2'-F, LNA, mixed)
- Gapmer architecture (wing lengths, gap length)
- 5'/3' end modifications (phosphate, cap, etc.)
- Conjugation type (unconjugated, GalNAc, cholesterol, lipid, antibody)
- Formulation/delivery vehicle (if applicable)

**Physicochemical Features:**
- Molecular weight
- Net charge at physiological pH
- Predicted logP / hydrophobicity index
- Predicted plasma protein binding affinity (in silico)

### 3.3 Synthesis Quality Standards

All oligonucleotides will be synthesized by two validated CRO partners (Integrated DNA Technologies and Eurofins Genomics) and subjected to:
- RP-HPLC purity assessment (target: >95% for standard oligos, >90% for heavily modified)
- Mass spectrometry confirmation of molecular identity
- Endotoxin testing (<1 EU/mL in working solution)
- Sterile filtration
- Lyophilization and reconstitution in RNase-free PBS

---

## 4. Experimental Systems & Model Selection

We selected five complementary human in vitro model systems based on four criteria: (1) relevance to the target organ/mechanism, (2) established use in oligo toxicology, (3) scalability for high-throughput screening, and (4) availability from qualified suppliers with batch consistency guarantees.

### 4.1 System 1: Primary Human Hepatocytes (PHH) in Sandwich Culture

**Relevance:** Gold standard for hepatotoxicity assessment; retains CYP enzyme activity, transport protein expression, and hepatocyte polarity critical for oligo uptake and metabolism.  
**Source:** BioIVT (3 independent donors per batch run; n=3 biological replicates)  
**Dosing:** 6 concentrations (0.01, 0.1, 1, 10, 50, 100 µM); 24h and 72h timepoints  
**Readouts:** Cell viability (ATP-lite), LDH release, ALT/AST secretion, caspase 3/7 activity, ROS generation, lipid accumulation (Oil Red O), bile acid transport (CDFDA), transcriptomics (3' RNA-seq at 1 µM, 10 µM for cytotoxic compounds)

### 4.2 System 2: Human Kidney Proximal Tubule Organoids

**Relevance:** Proximal tubule is the primary site of oligo-induced nephrotoxicity. Organoids capture tubular morphology and transporter expression better than standard 2D cultures.  
**Source:** Hubrecht Organoid Technology (HUB) certified human kidney organoids, differentiated to proximal tubule identity (SLC22A2+, LRP2+)  
**Dosing:** 5 concentrations (0.1, 1, 10, 50, 100 µM); 48h and 96h timepoints  
**Readouts:** Organoid viability (CellTiter-Glo 3D), KIM-1 secretion (ELISA), NGAL secretion, caspase activation, tubular injury morphology (brightfield + AI image analysis), lysosomal accumulation (LysoTracker), transcriptomics subset (n=300 highest priority oligos)

### 4.3 System 3: Human PBMC Co-culture (Immunotoxicity & Complement Panel)

**Relevance:** PBMCs capture innate immune activation (monocytes, NK cells, B cells) and T cell responses; the most relevant system for CpG- and PS-backbone-driven immunostimulation.  
**Source:** Fresh PBMCs from 5 healthy donors (Hemacare); 3 donors per experiment run (rotation)  
**Dosing:** 4 concentrations (0.1, 1, 10, 50 µM); 24h timepoint  
**Readouts:** Cytokine multiplex (28-plex Luminex: IFN-α, IFN-γ, IL-1β, IL-6, IL-8, IL-10, IL-12p70, TNF-α, MCP-1, IP-10, and 18 others), complement C3a and C5a (ELISA), complement C1q binding (direct binding assay), TLR9 activation (NF-κB reporter in HEK-Blue hTLR9 cells, run in parallel), NK cell activation (CD69 flow cytometry)

### 4.4 System 4: Human Platelet-Rich Plasma (Thrombocytopenia & Coagulopathy)

**Relevance:** Direct assessment of platelet aggregation, activation, and coagulation factor interaction — the mechanisms underlying PS-backbone-induced thrombocytopenia and coagulopathy.  
**Source:** Fresh human citrated plasma (5 donors); platelets isolated by differential centrifugation  
**Dosing:** 3 concentrations (1, 10, 50 µM); 30 min and 2h incubations  
**Readouts:** Platelet aggregation (light transmission aggregometry, ADP-induced), platelet activation markers (CD62P/P-selectin, CD63 by flow cytometry), aPTT and PT clotting times (automated coagulometer), Factor Xa and thrombin inhibition (chromogenic substrate assays), complement activation in plasma (C3a, sC5b-9)

### 4.5 System 5: Liver-Kidney Microphysiological System (MPS/Organ-on-Chip)

**Relevance:** Multi-organ MPS captures systemic distribution and sequential organ toxicity better than single-organ models, most closely mimicking in vivo pharmacokinetics relevant to chronic dosing.  
**Source:** Emulate, Inc. Liver-Chip and Kidney-Chip platforms (n=3 chips per condition)  
**Dosing:** 3 concentrations (1, 10, 50 µM); 72h and 7-day timepoints  
**Readouts:** Effluent sampling (ALT, AST, albumin, KIM-1, creatinine), chip morphology imaging, barrier integrity (TEER), transcriptomics from chip lysates, proteomics from effluent (nano-LC-MS/MS, 300 most impactful oligos)

### 4.6 Tiered Screening Strategy

Not all 2,800 oligos will be run through all 5 systems at full resolution. We use a tiered strategy to maximize scientific value within budget:

| Tier | N oligos | Systems | Depth |
|---|---|---|---|
| Tier 1 (full panel) | 600 | All 5 systems | Full dose-response, all readouts |
| Tier 2 (priority panel) | 1,200 | PHH + PBMC + Platelet | Dose-response, standard readouts |
| Tier 3 (high-throughput) | 1,000 | PHH viability + PBMC cytokines | 3-point dose, key readouts |

Active learning directs which Tier 3 compounds graduate to Tier 2 or Tier 1 based on predicted information value (see §8.3).

---

## 5. Toxicology Assay Panel

### 5.1 Complete Readout Inventory (47 endpoints)

**Hepatotoxicity Panel (PHH)**
1. Cell viability — ATP-lite luminescence
2. Membrane integrity — LDH release
3. ALT secretion (hepatocyte-specific injury marker)
4. AST secretion
5. Apoptosis — Caspase 3/7 activity (Caspase-Glo)
6. Necrosis marker — HMGB1 release
7. Oxidative stress — CellROX fluorescence
8. Lipid accumulation — Oil Red O staining (image quantified)
9. Mitochondrial membrane potential — JC-1 ratio
10. Bile transport function — CDFDA efflux assay
11. Gene expression module score — hepatotoxicity gene panel (50 gene NanoString panel)
12. Transcriptome-wide RNA-seq (subset: n=600 oligos at 1 and 10 µM)

**Nephrotoxicity Panel (Kidney Organoids)**
13. Organoid viability — CellTiter-Glo 3D
14. KIM-1 secretion (tubular injury biomarker)
15. NGAL secretion (tubular injury biomarker)
16. Apoptosis — Caspase 3/7
17. Lysosomal accumulation — LysoTracker fluorescence
18. Tight junction integrity — ZO-1 immunofluorescence (imaging)
19. Mitochondrial stress — Mitotracker signal
20. Tubular injury transcriptomics — targeted 30-gene panel (NanoString)

**Immunotoxicity & Complement Panel (PBMC)**
21. IFN-α (innate antiviral response)
22. IFN-γ
23. IL-1β
24. IL-6
25. IL-8 (CXCL8)
26. IL-10
27. IL-12p70
28. TNF-α
29. MCP-1 (CCL2)
30. IP-10 (CXCL10)
31. GM-CSF
32. IL-4, IL-5, IL-13 (Th2)
33. IL-17A (Th17)
34. TLR9 pathway activation (NF-κB reporter, HEK-hTLR9)
35. NK cell activation — CD69 surface expression (flow)
36. Complement C3a (anaphylatoxin)
37. Complement C5a (anaphylatoxin)

**Thrombocytopenia & Coagulopathy Panel (Plasma/Platelets)**
38. Platelet aggregation (% max aggregation, ADP-induced)
39. P-selectin expression (CD62P) — platelet activation
40. CD63 expression — platelet degranulation
41. Activated Partial Thromboplastin Time (aPTT)
42. Prothrombin Time (PT)
43. Factor Xa inhibition (chromogenic)
44. Thrombin inhibition (chromogenic)

**Multi-organ MPS Panel (Organ-on-Chip)**
45. Chip effluent ALT + AST (hepatotoxicity under flow)
46. Chip effluent KIM-1 + creatinine (nephrotoxicity under flow)
47. Multi-omic profiling (proteomics/transcriptomics, subset n=300)

### 5.2 Dose-Response Data Structure

For Tier 1 and Tier 2 compounds, each readout is collected as a complete dose-response curve, enabling calculation of:
- **IC50 / EC50** with 95% CI (4-parameter logistic fit)
- **Maximum effect (Emax)** at highest tested concentration
- **Hill coefficient** (slope/cooperativity)
- **No Observed Effect Concentration (NOEC)**
- **Area Under the Dose-Response Curve (AUDRC)** — a single summary statistic suitable for regression modeling

---

## 6. Data Generation Plan

### 6.1 Experimental Timeline

Data generation is organized into five sequential batches (Batches A–E), each processing ~560 oligos. Batch structure allows rolling data release, early model training, and active learning iteration between batches.

| Batch | Oligos | Period | Key Content |
|---|---|---|---|
| A | 560 | Jun–Jul 2026 | Controls, backbone series, CpG series |
| B | 560 | Jul–Aug 2026 | Sugar modification series, length series |
| C | 560 | Aug–Sep 2026 | GalNAc conjugates, siRNA panel |
| D | 560 | Sep–Oct 2026 | Active learning-selected compounds |
| E | 560 | Oct–Nov 2026 | Exploratory + active learning round 2 |

### 6.2 Workflow Summary

1. **Compound management**: Oligos received, QC-tested, plated into acoustic dispensing-compatible format (Echo 555)
2. **Dosing**: Acoustic liquid handler (Labcyte Echo) for consistent nanoliter dispensing
3. **Cell plating**: Automated liquid handling (Hamilton STAR) for PHH, PBMC, platelet, and organoid assays
4. **Incubation & measurement**: Automated plate readers (Perkin Elmer EnVision, BioTek Cytation 5 for imaging)
5. **Data capture**: Automated LIMS integration → structured data tables per assay per compound
6. **QC pipeline**: Automated Z'-factor, positive/negative control performance checks per plate
7. **Upload**: Weekly data uploads to OligoToxDB repository with provenance metadata

### 6.3 Replication Strategy

- **Biological replicates**: n=3 (independent cell donors or organoid preparations) for PHH, PBMC, and kidney organoid assays
- **Technical replicates**: n=2 per biological replicate for all Tier 1 endpoints
- **Inter-batch reproducibility**: 20 reference compounds repeated in each batch to monitor drift
- **Inter-lab reproducibility**: 50-compound subset shared with one external partner lab for independent replication

---

## 7. Data Schema & Annotation Standards

### 7.1 OligoToxDB Data Model

Each database record is keyed by a unique **OligoTox ID** (format: OT-XXXXX). The schema is organized into four linked tables:

**Table 1: Compound Registry**
```
OligoTox_ID      | Unique identifier
Sequence_5to3    | Full nucleotide sequence
Length           | Integer (nt)
Backbone_class   | Enumerated (PS, PO, PMO, LNA, PNA, etc.)
Sugar_mod        | Enumerated (DNA, 2'OMe, 2'F, LNA, mixed)
Gapmer_design    | String (e.g., "5-10-5 LNA gapmer")
Conjugate        | Enumerated (none, GalNAc, cholesterol, lipid)
Synthesis_vendor | String
Purity_percent   | Float
Mass_confirmed   | Boolean
Endotoxin_EU_mL  | Float
```

**Table 2: Computed Features**
```
OligoTox_ID      | Foreign key
GC_content       | Float (%)
CpG_count        | Integer
MFE_kcal_mol     | Float (RNAfold)
Tm_target        | Float (°C)
Tm_offtarget_max | Float (°C, top BLAST hit)
G4_score         | Float (QGRS)
Self_comp_index  | Float
Hybridization_dG | Float (kcal/mol)
Offtarget_score  | Float (log10 BLAST E-value)
MW_Da            | Float
Net_charge_pH7   | Integer
```

**Table 3: Experimental Results**
```
OligoTox_ID        | Foreign key
Assay_system       | Enumerated (PHH, KidneyOrganoid, PBMC, Platelet, MPS)
Endpoint           | String (e.g., "Cell_viability_ATPLite")
Concentration_uM   | Float
Timepoint_h        | Integer
Replicate_bio      | Integer (1, 2, 3)
Replicate_tech     | Integer (1, 2)
Batch_ID           | String
Cell_donor_ID      | Anonymized donor identifier
Value_raw          | Float
Value_normalized   | Float (% of vehicle control)
Unit               | String
QC_flag            | Enumerated (pass, marginal, fail)
```

**Table 4: Derived Summary Statistics**
```
OligoTox_ID      | Foreign key
Assay_system     | Foreign key
Endpoint         | Foreign key
IC50_uM          | Float
IC50_CI_low      | Float
IC50_CI_high     | Float
Emax_percent     | Float
Hill_coefficient | Float
NOEC_uM          | Float
AUDRC            | Float
N_bio_reps       | Integer
Toxicity_class   | Enumerated (inactive / low / moderate / high)
```

### 7.2 Controlled Vocabularies & Ontologies

All fields use controlled vocabularies aligned with established ontologies:
- **Chemical entities**: ChEBI ontology; oligo-specific extensions
- **Cell types**: Cell Ontology (CL)
- **Assay types**: Bioassay Ontology (BAO)
- **Organisms**: NCBI Taxonomy
- **Toxicity endpoints**: NCI Thesaurus (NCIt)
- **Genes/transcripts**: HGNC approved symbols; Ensembl IDs
- **FAIR compliance**: Each dataset assigned a persistent DOI via Zenodo; metadata in Dublin Core + DataCite schema

### 7.3 Metadata Standards

All experimental data is accompanied by minimum metadata following:
- **MIAME** (Minimum Information About Microarray Experiments) for transcriptomics
- **MIAPE** (Minimum Information About a Proteomics Experiment) for proteomics
- **CDISC SEND** (Standard for Exchange of Nonclinical Data) for toxicology study metadata
- Custom **OligoTox Minimum Reporting Standard (OTMRS)** — a lightweight standard we will publish as part of this work — specifying required fields for oligo identity, cell model, dosing conditions, and readout annotation

---

## 8. AI/ML Framework: OligoTox-Predict

### 8.1 Overview

OligoTox-Predict is a suite of open-source machine learning models trained on OligoToxDB. We release three model architectures alongside the dataset, each optimized for a different use case: discovery screening, mechanistic understanding, and uncertainty-aware lead optimization.

All models, training code, hyperparameter configurations, and evaluation benchmarks are released on GitHub under Apache 2.0 license. Trained model weights are released on HuggingFace Hub.

### 8.2 Model 1: OligoTox-Transformer (Sequence Foundation Model Fine-Tune)

**Architecture:** We fine-tune the Nucleotide Transformer (NT-v2, 500M parameters; Dalla-Torre et al. 2023) on OligoToxDB toxicity labels using task-specific regression heads.

**Why this approach:** Nucleotide Transformer was pre-trained on 3,000+ genomes and encodes deep representations of nucleotide co-occurrence patterns, secondary structure correlates, and evolutionary conservation. Fine-tuning on toxicity labels transfers these representations to the oligo safety prediction task, enabling strong performance even for chemical classes with limited training data.

**Input:** Raw nucleotide sequence (tokenized as 6-mers)  
**Output:** 47-dimensional toxicity vector (one regression value per endpoint) + per-endpoint uncertainty estimate (MC-Dropout, n=50 forward passes)

**Modifications:**
- Chemical modification tokens appended as a separate encoding track, fused at the embedding layer via cross-attention
- Multi-task training head: all 47 endpoints trained jointly with task-specific loss weighting (inverse-frequency weighting for rare toxic compounds)
- Uncertainty head: Monte Carlo Dropout for epistemic uncertainty; predicted variance for aleatoric uncertainty

**Interpretability:** Attention rollout and Integrated Gradients to identify nucleotide positions most responsible for toxicity predictions — generating testable mechanistic hypotheses (e.g., "positions 3-5 in the gap of this gapmer class drive hepatotoxicity signal").

### 8.3 Model 2: OligoTox-XGB (Gradient Boosted Feature Model)

**Architecture:** XGBoost with hand-engineered features from Table 2 (Compound Registry + Computed Features).

**Why this approach:** Interpretable, computationally lightweight, and deployable without GPU infrastructure. Suitable for rapid in-silico prescreening during early compound design, before synthesis. Performs best for interpolation within well-sampled chemical space.

**Input:** 45 computed features (sequence, modification, physicochemical)  
**Output:** Predicted IC50 per endpoint (6 primary organ/mechanism predictions); toxicity class probability

**Key features:** SHAP values computed for every prediction — researchers can see exactly which molecular features drove the prediction, with confidence intervals.

### 8.4 Model 3: OligoTox-ActiveLearn (Active Learning Loop)

**Architecture:** Gaussian Process (GP) ensemble with deep kernel learning (DKL-GP), wrapped in a Bayesian optimization acquisition function.

**Purpose:** This is not primarily a standalone predictor — it is the engine that guides experimental design across Batches D and E. Given the current dataset, OligoTox-ActiveLearn selects the next 200 compounds to synthesize and screen by maximizing **Expected Information Gain (EIG)** — prioritizing regions of chemical space where the model is most uncertain and where new data points would most reduce overall prediction error.

**Acquisition function:** We use a multi-objective Bayesian optimization (BO) with a Pareto front across:
1. Prediction uncertainty (epistemic) — explore unknown space
2. Predicted toxicity diversity — ensure coverage of active and inactive regions
3. Chemical space distance from existing data — avoid redundant screening

**Impact:** Simulation studies on pilot data suggest active learning can achieve equivalent model performance to random screening with 30-40% fewer compounds — saving ~$120,000 in experimental costs while enabling broader chemical space exploration.

### 8.5 Benchmark & Evaluation

We will release three standardized benchmark splits for community model evaluation:

| Benchmark | Split type | Task | Metric |
|---|---|---|---|
| OligoTox-HepBench | Random 80/10/10 | Hepatotoxicity classification | AUROC, F1 |
| OligoTox-ChemSplit | By backbone class | Cross-chemistry generalization | AUROC per class |
| OligoTox-ProspBench | Time-split (Batches A-C train, D-E test) | Prospective prediction | Spearman r, RMSE |

Benchmarks follow the MoleculeNet/ADMET-AI standard format and will be hosted on the HuggingFace Datasets Hub with a public leaderboard.

### 8.6 Interactive Prediction Portal

We will deploy **OligoTox-Web** — a publicly accessible web application (hosted on Hugging Face Spaces) where any researcher can:

1. **Query the database**: Search OligoToxDB by sequence, chemistry, or toxicity profile
2. **Run predictions**: Input any oligo sequence + modification scheme → receive predicted toxicity profile across all 47 endpoints with uncertainty estimates
3. **Visualize results**: Dose-response curve viewer, chemical space PCA plot, SHAP feature importance visualization
4. **Download data**: Export filtered datasets as CSV, SDF, or HDF5

The portal will not require login for basic access; researchers can create accounts to save predictions and contribute to a federated annotation layer.

---

## 9. Quality Assurance & Reproducibility

### 9.1 Experimental QC Standards

**Plate-level QC (applied to every plate):**
- Z'-factor ≥ 0.5 (calculated from positive/negative controls)
- Coefficient of Variation (CV) of vehicle controls ≤ 15%
- Positive control response within 2 SD of historical mean
- Plates failing QC are repeated; chronic failure triggers assay revalidation

**Compound-level QC:**
- Outlier detection: Grubbs test across biological replicates (α = 0.05)
- Dose-response curve quality: R² ≥ 0.85 for fitted curves; flagged otherwise
- Purity re-check: Any compound with unexpected toxicity profile re-confirmed by MS

**Batch-level QC:**
- 20 anchor compounds repeated in each batch
- Inter-batch CV < 20% for all anchor compounds; drift corrected via ComBat normalization if needed
- Transcriptomics: RNA integrity (RIN) ≥ 8.0 for all RNA-seq samples

### 9.2 Reproducibility Package

For each published dataset batch, we release:
- **Experimental protocols**: Full SOPs as PDF (also deposited in protocols.io under CC-BY license)
- **Raw data**: Unprocessed plate reader files in native format
- **Processing code**: Python scripts for all data normalization, QC, and curve fitting
- **Analysis notebooks**: Jupyter notebooks reproducing all figures in data releases
- **Container image**: Docker image with all software dependencies for full reproducibility

### 9.3 External Validation

A 50-compound blind validation subset will be sent to one external laboratory (Leiden University's Division of Systems Biomedicine and Pharmacology, a leader in oligo toxicology) for independent replication. Results will be published alongside the primary dataset release, with inter-lab agreement statistics reported.

---

## 10. Public Access & Dissemination Plan (PADP)

### 10.1 Data Release Philosophy

OligoToxDB will be maximally open. There are no proprietary restrictions, no embargo periods beyond the rolling release schedule below, and no registration barriers for data download.

### 10.2 Data Repositories & Formats

| Data type | Repository | Format | License |
|---|---|---|---|
| Compound registry + assay results | Zenodo (DOI-assigned) | CSV, HDF5, Parquet | CC BY 4.0 |
| Transcriptomics (RNA-seq) | GEO (NCBI) | FASTQ (raw), TPM matrix | CC BY 4.0 |
| Proteomics (MS/MS) | PRIDE Archive (EBI) | mzML (raw), protein abundance CSV | CC BY 4.0 |
| Cytokine profiles | Zenodo | CSV | CC BY 4.0 |
| ML models & weights | HuggingFace Hub | PyTorch .pt + ONNX | Apache 2.0 |
| Training code | GitHub (anote-ai/oligotox-predict) | Python | Apache 2.0 |
| Analysis notebooks | GitHub + Binder | Jupyter | Apache 2.0 |
| SOPs & protocols | protocols.io | PDF | CC BY 4.0 |

**Data formats are chosen for maximum interoperability:**
- **CSV/Parquet**: Universal; readable in Python, R, Excel without specialized tools
- **HDF5**: Efficient storage for large multimodal arrays; readable via h5py, rhdf5
- **FASTQ/mzML**: Community-standard raw data formats for transcriptomics/proteomics
- **ONNX**: Hardware-agnostic model format; runs in Python, JavaScript, C++, Java

### 10.3 Rolling Release Schedule

Data will be released in rolling batches — not held until project completion:

| Milestone | Release date | Content |
|---|---|---|
| Batch A complete | August 15, 2026 | ~560 compounds, full assay panel |
| Batch B complete | September 15, 2026 | Additional ~560 compounds |
| Batch C complete | October 15, 2026 | Additional ~560 compounds + RNA-seq |
| Batch D complete | November 1, 2026 | Active learning batch + proteomics |
| Batch E complete | November 30, 2026 | Final batch + all models + portal |
| Final submission | December 15, 2026 | Complete integrated dataset, paper preprint |

### 10.4 FAIR Data Principles Implementation

**Findable:**
- All datasets assigned persistent DOIs via Zenodo
- Metadata deposited in DataCite and registered with re3data.org
- OligoToxDB indexed in FAIRsharing.org registry
- Cross-referenced in PubChem, ChEMBL (compound registry)

**Accessible:**
- All data downloadable without login at primary repositories
- API access for programmatic queries (Zenodo REST API; custom OligoToxDB API)
- Mirror copies on OSF (Open Science Framework) for redundancy
- Persistent access guaranteed through Zenodo's long-term preservation commitment

**Interoperable:**
- Controlled vocabularies and ontologies for all fields (§7.2)
- Standard formats (CSV, HDF5, FASTQ, mzML)
- Linked to PubChem CIDs and ChEMBL IDs for cross-database queries
- SPARQL endpoint for RDF-serialized metadata

**Reusable:**
- CC BY 4.0 license on all datasets: attribution required, all uses permitted including commercial
- Apache 2.0 on all code: permissive, patent-grant included
- Detailed data dictionaries, codebooks, and usage examples included
- Contact point for data questions: nvidra@anote.ai

### 10.5 Community Engagement & Dissemination

**Scientific publications:**
- Data descriptor paper submitted to *Scientific Data* (Nature) upon Batch A release
- Methods paper submitted to *Nucleic Acid Therapeutics* or *Molecular Therapy – Nucleic Acids*
- Machine learning paper submitted to *Nature Machine Intelligence* or *PLOS Computational Biology*
- All papers submitted to bioRxiv as preprints simultaneously

**Conference presentations:**
- Oligonucleotide Therapeutics Society (OTS) Annual Meeting, Fall 2026 — poster + oral presentation
- Society of Toxicology (SOT) Annual Meeting, 2027 — oral presentation
- NeurIPS 2026 ML for Drug Discovery workshop — poster submission

**Community workshops:**
- Free virtual workshop (Q4 2026): "Using OligoToxDB for oligo safety AI" — hands-on tutorial with Jupyter notebooks
- Office hours: Monthly 1-hour Q&A sessions for researchers using the dataset

**Direct outreach:**
- Direct notification to >50 oligonucleotide therapeutics labs globally (mailing list maintained)
- Partnership with Oligonucleotide Therapeutics Society to announce dataset to membership
- Direct outreach to pharma companies with active oligo programs (AstraZeneca, Ionis, Alnylam, Sarepta, Novartis Gene Therapies) encouraging dataset use and potential contribution of additional data

### 10.6 Long-term Sustainability

OligoToxDB is designed for long-term utility beyond the grant period:
- **Zenodo/EMBL-EBI long-term archiving**: Both platforms commit to 20+ year data preservation
- **Community contribution layer**: Researchers can submit additional oligo toxicity data following OTMRS standard; curated and integrated with versioned releases
- **Anote commitment**: Anote will maintain the prediction portal and GitHub repositories for a minimum of 3 years post-project. The portal uses infrastructure-as-code (Terraform) for easy handoff to community stewardship if needed.
- **Open governance**: An OligoToxDB Advisory Board (3-5 members from academia, industry, and regulatory) will guide long-term curation decisions

---

## 11. Team & Collaborators

### Core Team (Anote, Inc.)

| Role | Expertise |
|---|---|
| **Natan Vidra** (PI, Lead Innovator) | AI/ML, data platform, project coordination |
| **ML Engineer** (to be hired Q2 2026) | Deep learning, sequence modeling, model training |
| **Bioinformatician** (to be hired Q2 2026) | Transcriptomics/proteomics analysis, data pipeline |
| **Project Manager** | Timeline, vendor management, QA oversight |

### Experimental Partners (CRO)

| Partner | Role |
|---|---|
| **BioIVT** | Primary human hepatocyte supply and PHH sandwich culture assays |
| **Integrated DNA Technologies (IDT)** | Oligonucleotide synthesis (standard chemistries) |
| **Eurofins Genomics** | Oligonucleotide synthesis (modified chemistries, QC) |
| **Emulate, Inc.** | Liver-Chip and Kidney-Chip MPS assays |
| **HUB Organoids** | Kidney proximal tubule organoids |
| **Novogene** | RNA-seq library prep and sequencing |

### Scientific Advisory & Validation

| Partner | Role |
|---|---|
| **Leiden University (Dr. Manoharan lab)** | External validation replication; oligo chemistry expertise |
| **Karolinska Institutet (DILI-PREDICT consortium)** | Hepatotoxicity endpoint validation, PHH protocol QA |

### Regulatory Liaison

We will engage with FDA Center for Drug Evaluation and Research (CDER) staff scientists (Division of Applied Regulatory Science) to ensure OligoToxDB endpoints and data formats are maximally useful for regulatory submissions and IND-enabling studies. Two informal consultations planned: project kickoff (June 2026) and data release (November 2026).

---

## 12. Budget Overview

| Category | Estimated Cost |
|---|---|
| Oligonucleotide synthesis (2,800 compounds) | $280,000 |
| PHH assays (all batches) | $180,000 |
| Kidney organoid assays | $120,000 |
| PBMC / immunotoxicity assays | $90,000 |
| Platelet / coagulation assays | $60,000 |
| MPS / organ-on-chip assays | $150,000 |
| RNA-seq (600 samples) | $90,000 |
| Proteomics (300 samples) | $120,000 |
| Data infrastructure & portal | $40,000 |
| Personnel (ML engineer, bioinformatician) | $200,000 |
| Publications & dissemination | $20,000 |
| Project management & QA | $50,000 |
| **Total** | **~$1,400,000** |

*Note: This project is designed as a multi-funder effort. The requested Phase 2 prize contribution ($100,000–$200,000) will be combined with matched funding from Anote, Inc. and prospective co-funding from industry partners (discussions ongoing with two oligo therapeutic companies). The budget above reflects the full project scope; Phase 2 prize funds will be allocated to data generation and ML model development.*

---

## 13. Project Timeline

```
May 2026
├── Project kickoff, team hiring finalized
├── Vendor contracts executed (IDT, Eurofins, BioIVT, Emulate)
├── LIMS setup and data pipeline testing
└── Oligo library design finalized, orders placed (Batch A)

June 2026
├── Batch A synthesis complete; QC and plating
├── PHH Batch A assays begin
├── OligoTox-XGB trained on pilot data (n=120 reference compounds)
└── GitHub repository and Zenodo community space set up

July 2026
├── Batch A PHH + PBMC + Platelet assays complete
├── Batch B synthesis complete
├── Batch B assays begin
└── Batch A data release (Zenodo, CC BY 4.0)

August 2026
├── Batch B assays complete
├── Batch C synthesis complete
├── First RNA-seq run (Batch A subset, n=200)
├── OligoTox-Transformer fine-tuning begins (on Batch A+B data)
└── Batch B data release

September 2026
├── Batch C assays complete
├── Active learning round 1: OligoTox-ActiveLearn selects Batch D compounds
├── Batch D synthesis ordered
└── Batch C data release

October 2026
├── Batch D assays complete
├── Proteomics run (MPS effluent, n=300 subset)
├── OligoTox-Transformer training complete; model validation
├── Active learning round 2: Batch E compounds selected
└── Batch D data release

November 2026
├── Batch E assays complete
├── Full dataset integration, QC, harmonization
├── External validation results received from Leiden University
├── OligoTox-Web portal beta launch
└── All data released (Zenodo, GEO, PRIDE)

December 2026
├── Final submission to NIH OligoTox challenge
├── Scientific Data paper submitted
├── OTS conference presentation materials prepared
└── ML paper preprint posted to bioRxiv
```

---

## 14. Anticipated Impact

### For the Oligonucleotide Therapeutics Field

OligoToxDB will be the definitive public resource for oligo safety data. It directly enables:

- **In silico safety screening**: Companies and academics can pre-filter compound libraries before synthesis, reducing late-stage safety attrition
- **IND-enabling studies**: Regulators can use OligoToxDB as a benchmark for interpreting company-submitted oligo safety data
- **Mechanistic understanding**: The multimodal, multiendpoint nature enables causal analysis — why does a compound cause hepatotoxicity but not immunotoxicity? The data can answer this.
- **Reduce animal testing**: A well-validated predictive model trained on human in vitro data can substitute rodent toxicology screens early in discovery, aligning with the NIH 3Rs mandate

### For the AI/ML Community

OligoToxDB is purpose-built as an ML training resource. It provides:
- **Scale**: 2,800 compounds × 47 endpoints = 131,600 labeled data points (accounting for dose-response structure, orders of magnitude more)
- **Structure**: Clean, standardized schema; no missing data imputation required for primary endpoints
- **Benchmarks**: Three competition-ready benchmark splits for community model evaluation
- **Baselines**: Three released baseline models for comparison

### For NIH's Scientific Mission

This work directly advances:
- **NCATS mandate**: Accelerating translation of therapeutics by de-risking safety early
- **NIH 3Rs policy**: Replacing animal safety models with human-relevant in vitro systems
- **Open science priorities**: Fully FAIR, open-licensed, community-accessible data
- **Precision medicine**: Enables safety-aware oligo design for patient-specific gene therapy applications

### Projected Usage

Based on the current OTS membership (600+ organizations), PubMed citation rates for comparable public toxicology datasets (e.g., Tox21: >500 citations/year), and the growing oligonucleotide therapeutics field (14% CAGR), we project:
- **Year 1**: 50+ lab groups download and use OligoToxDB
- **Year 2**: 10+ published studies use OligoToxDB as training data or validation resource
- **Year 3**: OligoTox-Predict integrated into ≥2 commercial oligo design platforms

---

## Appendices

### Appendix A: Positive and Negative Control List

| Control | Type | Expected Response | Endpoint(s) |
|---|---|---|---|
| CpG ODN 2006 (TLR9 agonist) | Positive | High IFN-α, IL-6, TNF-α | Immunotoxicity |
| Poly-G ODN (G-quartet) | Positive | Complement C3a ↑, platelet aggregation ↑ | Complement, coagulopathy |
| ISIS 104838 | Positive | ALT/AST ↑, caspase ↑ (hepatotoxic ASO) | Hepatotoxicity |
| Nusinersen (Spinraza) | Negative | Minimal toxicity across all endpoints | All |
| Mipomersen (Kynamro) | Moderate positive | Mild hepatotoxicity signal | Hepatotoxicity |
| GalNAc-siRNA (inclisiran analog) | Negative | Minimal toxicity | Hepatotoxicity |
| Scrambled PS-ASO (non-targeting) | Negative/vehicle | Backbone-level baseline only | All |
| RNase-treated scrambled | Negative | No effect | All |

### Appendix B: Key References

1. Crooke, S.T. et al. (2021). Molecular Mechanisms of Antisense Oligonucleotides. *Nucleic Acid Therapeutics*, 31(1), 1-20.
2. Geary, R.S. et al. (2015). Pharmacokinetic Properties of 2'-O-(2-Methoxyethyl)-Modified Oligonucleotide Analogs in Rats. *J Pharmacol Exp Ther*, 355(3), 460-467.
3. Dalla-Torre, H. et al. (2023). The Nucleotide Transformer: Building and Evaluating Robust Foundation Models for Human Genomics. *bioRxiv*, 2023.01.11.523679.
4. Hagedorn, P.H. et al. (2018). Managing the Potential and Pitfalls During Translation of Sense and Antisense Oligonucleotides. *Nucleic Acid Therapeutics*, 28(3), 1-9.
5. Sewing, S. et al. (2016). Hepatotoxic Potential of Therapeutic Oligonucleotides Can Be Predicted from Their Sequence and Modification Pattern. *Nucleic Acid Therapeutics*, 26(5), 324-334.
6. Bhatt, D.K. et al. (2017). Plasma Protein Binding of Highly Bound Drugs as a Potential Confounding Factor in Basic Pharmacokinetic/Pharmacodynamic Modeling. *J Pharm Sci*, 106(6), 1395-1400.
7. Wilkinson, M.D. et al. (2016). The FAIR Guiding Principles for Scientific Data Management and Stewardship. *Scientific Data*, 3, 160018.

### Appendix C: OligoTox Minimum Reporting Standard (OTMRS) — Draft

The following fields are required for any experimental oligo toxicity record to be included in OligoToxDB or any dataset claiming OTMRS compliance:

**Compound identity (required):**
- Full sequence (5'→3')
- Backbone chemistry (controlled vocabulary)
- Sugar modification pattern
- Synthesis purity (%)
- Endotoxin level (EU/mL)

**Experimental conditions (required):**
- Cell model (cell type, source, passage or lot)
- Dosing concentrations (list all, µM)
- Vehicle/solvent
- Timepoint(s)
- Number of biological replicates (minimum n=2)

**Readout (required):**
- Endpoint name (controlled vocabulary)
- Normalization method
- Positive control response (mean ± SD)
- Negative control response (mean ± SD)
- QC acceptance criteria and pass/fail status

**Derived values (if reported):**
- IC50/EC50 with confidence interval
- Fitting method and software version
- R² of fitted curve

---

*Submission prepared by Natan Vidra, Anote, Inc. | nvidra@anote.ai*  
*Building on Phase 1 Ideation Submission, awarded by NIH NCATS OligoTox Open Data Challenge*
