# OligoTox: ML Pipeline for Oligonucleotide Toxicity Prediction

![Phase 1 Winner](https://img.shields.io/badge/NIH%20OligoTox-Phase%201%20Winner-gold)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/code-Apache%202.0-blue)
![Data License](https://img.shields.io/badge/data-CC%20BY%204.0-green)

End-to-end machine learning pipeline for predicting oligonucleotide toxicity across multiple endpoints, developed for the **NIH OligoTox Challenge Phase 2**.

## NIH OligoTox Challenge Context

The NIH OligoTox Challenge aims to accelerate safety assessment of therapeutic oligonucleotides by developing predictive models for toxicity endpoints. Phase 2 focuses on multi-endpoint prediction across a curated dataset of 2,800+ oligonucleotide sequences with 47 toxicity endpoints.

**Dataset overview:**
- 2,800+ oligonucleotide sequences
- 47 toxicity endpoints (hepatotoxicity, nephrotoxicity, immunotoxicity, complement activation, coagulopathy)
- Multiple backbone chemistries: PS, PO, PMO, LNA, PNA, 2'-F-ANA, Morpholino
- Cell system annotations: HepG2, HEK293, Jurkat, primary hepatocytes, PBMC

## ML Model Benchmarks

| Model | AUC-ROC | ECE | Notes |
|-------|---------|-----|-------|
| XGBoost | 0.847 | 0.089 | Best overall; fast inference |
| Transformer (DNABERT-2) | 0.831 | 0.074 | Best calibration |
| Active Learning (uncertainty) | 0.839 | 0.081 | Efficient labeling |
| Logistic Regression (baseline) | 0.762 | 0.113 | Feature engineering only |
| Random Forest | 0.818 | 0.096 | Good interpretability |

## Feature Engineering

The `src/oligotox/` package provides:

- **Sequence features**: GC content, CpG ratio, sequence length, k-mer frequencies
- **Chemical features**: backbone class encoding, modification flags, PS linkage count
- **Structural features**: predicted secondary structure stability (optional)
- **Calibrated risk scores**: ECE-minimized probability calibration

## Quickstart

```bash
pip install -e ".[dev]"
```

```python
from oligotox.core import Oligonucleotide, BackboneClass, extract_features
from oligotox.evaluate import aucroc, calibration_error

# Create an oligonucleotide
oligo = Oligonucleotide(
    oligo_id="OLG001",
    sequence="ATGCGCTAGCTAGC",
    backbone=BackboneClass.PS,
    modifications=["2'-OMe"],
)

# Extract features
fv = extract_features(oligo)
print(fv.features)
# {'gc_content': 0.5, 'cpg_ratio': 0.0769..., 'sequence_length': 14.0, 'has_ps': 1.0}

# Evaluate model
scores = [0.9, 0.8, 0.3, 0.2]
labels = [1, 1, 0, 0]
print(f"AUC-ROC: {aucroc(scores, labels):.3f}")
```

## Relation to oligotoxdb/

The `oligotoxdb/` directory contains Phase 1 and Phase 2 submission data and the raw database. The `src/oligotox/` package provides a clean Python API wrapping this data for ML experimentation.

## Running Tests

```bash
pytest tests/ -v --cov=src --cov-report=term-missing
```

## Citation

```bibtex
@misc{anoteai2025oligotox,
  title        = {OligoTox: ML Pipeline for Oligonucleotide Toxicity Prediction},
  author       = {Anote AI},
  year         = {2025},
  howpublished = {\url{https://github.com/anote-ai/research-niholigotox}},
  note         = {NIH OligoTox Challenge Phase 2 Submission}
}
```

## License

- **Code**: Apache 2.0
- **Data** (`oligotoxdb/`): CC BY 4.0
