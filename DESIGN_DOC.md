# NIHOligoTox — Research Design Document

## Goal

Build and evaluate a machine learning model that predicts hepatotoxicity, nephrotoxicity, and immunotoxicity of oligonucleotide therapeutics from sequence and chemistry features — enabling early-stage computational screening to prioritize candidates before expensive wet-lab toxicity testing.

## Objective

1. Assemble the largest publicly-available dataset of oligonucleotide-toxicity pairs from ChEMBL, published literature, and public regulatory submissions
2. Train and evaluate a graph neural network / transformer hybrid model that represents both nucleotide sequence and backbone chemistry modifications
3. Demonstrate multi-task learning across toxicity endpoints outperforms single-endpoint models, and that the model generalizes to novel chemistries

## Background / Motivation

Oligonucleotide therapeutics are one of the fastest-growing drug classes (12 approved in the US since 2016, 100+ in clinical trials). The biggest clinical failure mode is unexpected toxicity, particularly hepatotoxicity. Computational toxicity prediction for small molecules has advanced significantly (DeepTox, MPNN), but oligonucleotides have unique properties that invalidate small-molecule approaches: sequence-based features and backbone chemistry modifications are critical and not captured by existing tools.

## Experimental Design

### Baseline Experiment

**Evaluate 3 molecular ML baselines (ECFP fingerprint + XGBoost, DeepTox, standard GNN with SMILES) on a held-out test set of 200 oligonucleotides with known hepatotoxicity labels**

- Metric: AUROC, AUPRC, sensitivity at 95% specificity
- Purpose: confirm existing small-molecule ML approaches underperform on oligonucleotides
- Expected result: AUROC ≈ 0.65–0.75 for baselines

### Test Experiment 1: Sequence + Chemistry Representation

Build a model with: (1) sequence encoder = transformer over nucleotide sequence with modification tokens, (2) chemistry encoder = GNN over backbone molecular graph, combined with cross-attention. Train on hepatotoxicity labels.

**Expected result:** combined encoder achieves AUROC ≈ 0.85+ vs. baseline 0.65–0.75 — 15+ point improvement demonstrating that backbone chemistry features contain toxicity signal SMILES representations miss

### Test Experiment 2: Multi-Task Learning Across Toxicity Endpoints

Train separate single-task models and a joint multi-task model for hepatotoxicity, nephrotoxicity, and immunotoxicity. Compare AUROCs per endpoint.

**Expected result:** multi-task model improves nephrotoxicity AUROC by 5–10 points vs. single-task (nephrotoxicity has fewer training examples; sharing with hepatotoxicity helps)

### Test Experiment 3: Out-of-Distribution Generalization

Hold out all training examples with one modification type (e.g., all LNA-containing oligonucleotides). Train on remainder. Evaluate on held-out modification type.

**Expected result:** AUROC drops from 0.85 to ~0.70 for held-out chemistry type — an honest limitation to document, with uncertainty quantification flagging predictions on novel chemistries as low-confidence

## Expected Results

1. The largest public oligonucleotide toxicity dataset
2. A trained multi-task GNN/transformer model with open-source weights
3. **Key finding:** "Backbone chemistry modifications carry as much toxicity signal as nucleotide sequence — models that ignore chemistry underperform by 15+ AUROC points"
4. Multi-task learning benefit: +5–10 AUROC on data-scarce endpoints
5. A practical screening tool with calibrated confidence and OOD detection

## Why This Matters / Why People Would Care

- **Oligonucleotide drug developers** (Ionis, Alnylam, Novo Nordisk): reducing toxicity-related failures could save $10M–$100M+ per program
- **NIH and FDA:** a public, open-source toxicity model reduces barriers for academic researchers and small biotechs
- **Computational chemists:** oligonucleotide ML is an underexplored area adjacent to the well-funded small-molecule ML space
- **Drug discovery AI field:** the sequence + chemistry multi-modal representation is generalizable to other modified nucleic acid modalities

## Timeline

| Month | Milestone |
|---|---|
| 1–2 | Dataset assembly (literature + ChEMBL + public regulatory submissions) |
| 3 | Model implementation (sequence transformer + chemistry GNN + cross-attention) |
| 4 | Baseline + test experiments |
| 5 | Calibration + OOD analysis |
| 6 | Submission to J. Chem. Info. & Modeling or NeurIPS ML4H |

## Related Issues

- Design doc GitHub issue: #23
- Target conferences: see issues labeled `conference-prep`
- Reproducibility package: see issues labeled `artifact-release`
