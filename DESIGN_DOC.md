# Research Design Document: NIH Oligotoxicology

## Vision Statement

Build the first **multi-modal, multi-task deep learning framework** for predicting oligonucleotide toxicity that jointly models chemical structure, target sequence, and off-target binding — demonstrating that sequence+chemistry representation improves AUROC by ≥15 points over chemistry-only models and establishing **OligoTox-Bench** as the standard evaluation suite for AI-driven oligonucleotide safety prediction.

---

## Problem Statement & Novelty

Oligonucleotides (ASOs, siRNAs, aptamers) are a rapidly growing drug class, but their toxicity prediction remains a critical bottleneck in drug development:

1. **Chemistry-only models miss sequence-specific effects**: Hepatotoxicity from ASOs is driven by both chemical modifications (backbone, 2' sugar) AND sequence-specific hybridization with off-target transcripts — but existing QSAR models use chemistry features only.
2. **Sparse labeled data**: Human clinical toxicity data is limited; models must leverage in vitro assays, animal studies, and cross-species transfer.
3. **Multi-endpoint prediction gap**: Each toxicity endpoint (hepatotox, nephrotox, immunostimulation, complement activation) is typically modeled separately, ignoring shared biological mechanisms.
4. **Scaffold bias in benchmarks**: Existing toxicity benchmarks have severe train/test scaffold overlap, leading to optimistic evaluation.

### Novel Contributions

| Contribution | Description |
|---|---|
| **OligoTox-Bench** | 3,200 curated oligonucleotide-toxicity data points across 6 endpoints, with scaffold splits |
| **SeqChem representation** | Combined sequence (k-mer embeddings) + chemistry (ECFP4 + modification features) representation |
| **Multi-task OligoTox model** | Jointly predicts 6 toxicity endpoints with shared representation |
| **Off-target predictor integration** | First model incorporating predicted off-target binding profiles into toxicity prediction |
| **Scaffold-split evaluation** | Rigorous OOD evaluation using Murcko scaffold splits to prevent data leakage |

### SeqChem Representation

```python
# Combined representation
def seqchem_representation(oligonucleotide):
    # Sequence component
    seq_features = kmer_embedding(oligonucleotide.sequence, k=3)  # 64-dim
    offtarget_features = predict_offtarget_binding(oligonucleotide)  # top-50 targets
    
    # Chemistry component  
    chem_features = ECFP4(oligonucleotide.smiles, radius=2)  # 2048-dim
    mod_features = modification_encoding(oligonucleotide.backbone_mods)  # 32-dim
    
    return concat([seq_features, offtarget_features, chem_features, mod_features])
```

---

## Research Objectives

1. Demonstrate that **SeqChem representation** improves AUROC by ≥15 pp over chemistry-only QSAR on scaffold-split evaluation.
2. Show that **multi-task learning** across 6 toxicity endpoints improves per-endpoint AUROC by ≥5 pp over single-task models.
3. Quantify the **off-target binding contribution**: how much does including predicted off-target profiles improve hepatotoxicity prediction?
4. Characterize **out-of-distribution generalization**: do improvements hold on scaffold-split OOD evaluation?
5. Build **OligoTox-Bench** as the standard evaluation dataset with rigorous train/test splits.

---

## Dataset Construction

### Data Sources

| Source | Data Type | Count | Endpoints |
|---|---|---|---|
| NIH toxicity database | In vivo (rat/mouse) | 800 | Hepatotox, nephrotox |
| ChEMBL bioassays | In vitro | 1,200 | Immunostim, complement |
| ClinicalTrials.gov | Human clinical | 400 | All endpoints (sparse) |
| Internal industry data | Mixed | 800 | All endpoints |

**Total: ~3,200 curated data points**

### Toxicity Endpoints

| Endpoint | Samples | Positive Rate | Prediction Difficulty |
|---|---|---|---|
| Hepatotoxicity | 1,100 | 28% | High |
| Nephrotoxicity | 900 | 19% | High |
| Immunostimulation (TLR9) | 800 | 35% | Medium |
| Complement activation | 700 | 22% | Medium |
| Thrombocytopenia | 400 | 12% | Very High |
| Injection site reactions | 300 | 31% | Medium |

### Data Curation Protocol

```
1. Deduplicate by canonical SMILES + sequence
2. Standardize SMILES (RDKit)
3. Curate endpoint labels: binary (toxic/non-toxic) + continuous (dose)
4. Scaffold split: Murcko scaffolds, 70/15/15 train/val/test
5. Verify: no scaffold overlap between train and test
6. Human expert review for data quality (random 10% sample)
```

---

## Models Under Evaluation

| Model | Representation | Multi-task | Off-target |
|---|---|---|---|
| QSAR (Random Forest) | Chemistry only (ECFP4) | No | No |
| DeepTox | Chemistry (CNN) | No | No |
| AttentiveFP | Molecular graph | No | No |
| SeqOnly | Sequence (k-mer) | No | No |
| SeqChem (ours, v1) | SeqChem | No | No |
| SeqChem-MT (ours, v2) | SeqChem | Yes | No |
| SeqChem-OT (ours, v3) | SeqChem + off-target | Yes | Yes |

---

## Experimental Design

### Baseline Experiment (Experiment 0)
**Protocol**: Train QSAR (Random Forest + ECFP4) on hepatotoxicity endpoint with random split. Compute AUROC.

**Expected result**: AUROC ≈ 0.79 with random split. This replicates published QSAR baseline and confirms our data pipeline.

---

### Experiment 1: Scaffold-Split Evaluation (OOD Generalization)
**Hypothesis**: All models show ≥10 pp AUROC drop under scaffold-split vs. random split; SeqChem models maintain better OOD performance than chemistry-only.

**Protocol**:
1. Train all models with both random and Murcko scaffold splits.
2. Compare AUROC: random vs. scaffold split per model.
3. Compute OOD degradation = random_AUROC - scaffold_AUROC.

**Expected results**:

| Model | Random Split AUROC | Scaffold Split AUROC | OOD Drop |
|---|---|---|---|
| QSAR (RF + ECFP4) | 0.79 | 0.62 | −17 pp |
| AttentiveFP | 0.83 | 0.67 | −16 pp |
| SeqChem (ours, v1) | 0.84 | 0.74 | −10 pp |
| SeqChem-OT (ours, v3) | 0.87 | 0.79 | −8 pp |

- Sequence features provide OOD robustness because hybridization patterns generalize across chemical scaffolds.

---

### Experiment 2: SeqChem vs. Chemistry-Only
**Hypothesis**: SeqChem representation improves AUROC by ≥15 pp over chemistry-only QSAR on scaffold-split evaluation for hepatotoxicity.

**Protocol**:
1. Systematic ablation: (a) Chemistry-only, (b) Sequence-only, (c) SeqChem combined.
2. Evaluate under scaffold split (strict OOD evaluation).
3. Statistical test: paired t-test across 5 random seeds.

**Expected results**:
- Chemistry-only AUROC: 0.62 (scaffold split)
- Sequence-only AUROC: 0.69 (scaffold split)
- SeqChem AUROC: 0.79 (scaffold split) — +17 pp over chemistry-only
- Both components contribute; combination is superadditive
- Key mechanism: sequence features capture off-target hybridization patterns that chemistry alone cannot

---

### Experiment 3: Multi-Task Learning Benefits
**Hypothesis**: Multi-task prediction across 6 endpoints improves mean per-endpoint AUROC by ≥5 pp vs. single-task, with largest improvements on data-sparse endpoints (thrombocytopenia, injection site).

**Protocol**:
1. Train: (a) 6 single-task models, (b) multi-task model (shared encoder).
2. Compare per-endpoint AUROC on scaffold-split test set.
3. Analysis: which endpoint pairs share the most representation (cosine similarity of learned features).

**Expected results**:

| Endpoint | Single-task AUROC | Multi-task AUROC | Improvement |
|---|---|---|---|
| Hepatotoxicity | 0.79 | 0.83 | +4 pp |
| Nephrotoxicity | 0.74 | 0.79 | +5 pp |
| Immunostimulation | 0.81 | 0.84 | +3 pp |
| Complement | 0.77 | 0.82 | +5 pp |
| Thrombocytopenia | 0.68 | 0.77 | +9 pp |
| Injection site | 0.71 | 0.79 | +8 pp |

- Sparse endpoints benefit most: shared representation transfers knowledge from data-rich endpoints.

---

### Experiment 4: Off-Target Binding Integration
**Hypothesis**: Adding predicted off-target binding profiles (top-50 transcripts with predicted binding affinity) improves hepatotoxicity AUROC by ≥5 pp.

**Protocol**:
1. Use TargetScan/RNAhybrid to predict off-target binding for each oligonucleotide.
2. Add off-target features to SeqChem-MT model.
3. Evaluate on hepatotoxicity (where off-target mechanisms are best characterized).
4. Ablation: how many off-target sites are needed (top-10, top-20, top-50)?

**Expected results**:
- SeqChem-MT without off-target: AUROC 0.83
- SeqChem-MT + top-10 off-target: AUROC 0.85 (+2 pp)
- SeqChem-MT + top-50 off-target: AUROC 0.88 (+5 pp)
- Saturates at ~50 off-target sites
- Key finding: off-target contribution is additive to sequence features (not redundant)

---

## Expected Results Summary

| Metric | Baseline (QSAR) | Our Best (SeqChem-OT) | Improvement |
|---|---|---|---|
| Hepatotox AUROC (scaffold) | 0.62 | 0.79 | +17 pp |
| Mean AUROC across 6 endpoints | 0.67 | 0.83 | +16 pp |
| Thrombocytopenia AUROC | 0.58 | 0.77 | +19 pp |
| OOD degradation | 17 pp | 8 pp | 2× better OOD robustness |

**Primary claim**: SeqChem representation combining sequence k-mer embeddings with chemistry features improves oligonucleotide toxicity prediction by 15–19 AUROC points on scaffold-split OOD evaluation, driven by sequence features capturing off-target hybridization mechanisms invisible to chemistry-only models.

---

## Why This Matters

**For researchers**: OligoTox-Bench provides the first rigorous benchmark for oligonucleotide toxicity prediction with proper scaffold splits.

**For drug development**: A 15 pp AUROC improvement in hepatotoxicity prediction could eliminate toxic compound classes earlier in development, reducing attrition and development costs.

**Market**: Oligonucleotide therapeutics is a $10B+ market growing rapidly; toxicity prediction tools are critical for pipeline safety.

**NIH/FDA relevance**: Aligns with NIH's push for in silico safety tools to reduce animal testing.

---

## Implementation Plan

```
research-niholigotox/
├── data/
│   ├── raw/             # Raw data from NIH, ChEMBL, etc.
│   ├── curated/         # 3,200 curated data points
│   └── splits/          # Scaffold train/val/test splits
├── representations/
│   ├── seqchem.py       # SeqChem combined representation
│   ├── offtarget.py     # Off-target binding predictor
│   └── chemistry.py     # ECFP4, modification features
├── models/
│   ├── single_task.py
│   ├── multi_task.py
│   └── baselines/       # QSAR, AttentiveFP, DeepTox
├── experiments/
│   ├── exp0_baseline.py
│   ├── exp1_scaffold_split.py
│   ├── exp2_seqchem.py
│   ├── exp3_multitask.py
│   └── exp4_offtarget.py
```

---

## Timeline

| Phase | Duration | Deliverable |
|---|---|---|
| Data curation | 6 weeks | 3,200 curated data points |
| Representation development | 4 weeks | SeqChem + off-target features |
| Model training + baseline | 4 weeks | All models trained |
| Experiments | 5 weeks | All results |
| Paper writing | 4 weeks | J. Chem. Info. submission |

**Target venues**: J. Chemical Information & Modeling, NeurIPS ML4H Workshop, or J. Cheminformatics

---

## Open Questions & Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Sparse positive labels for rare endpoints | High | Oversampling + cost-sensitive loss |
| Off-target predictor accuracy | Medium | Use ensemble of TargetScan + RNAhybrid |
| Cross-species extrapolation | High | Analyze species-specific subsets separately |
| Scaffold split reduces effective training size | Medium | Pre-training on broader toxicity data |

---

## Related Issues

- NIH grant alignment: FDA in silico safety initiative
- Reproducibility: scaffold split protocol
- Ethics: animal study data use policy
- Related work audit: DeepTox, Tox21, AttentiveFP
