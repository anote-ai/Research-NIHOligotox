# OligoToxDB

**Open dataset and AI prediction toolkit for oligonucleotide toxicity**

[![CI](https://github.com/anote-ai/nih-oligotox/actions/workflows/ci.yml/badge.svg)](https://github.com/anote-ai/nih-oligotox/actions)
[![License: Apache 2.0](https://img.shields.io/badge/Code-Apache%202.0-blue)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-green)](https://creativecommons.org/licenses/by/4.0/)
[![NIH NCATS OligoTox Challenge](https://img.shields.io/badge/NIH%20NCATS-OligoTox%20Phase%202-orange)](https://oligotox.com)

Part of the [NIH NCATS OligoTox Open Data Challenge](https://oligotox.com) Phase 2 submission.  
Developed by [Anote, Inc.](https://anote.ai) | PI: Natan Vidra

---

## What Is OligoToxDB?

OligoToxDB is the largest open dataset of **human in vitro oligonucleotide toxicity profiles**:

- **2,800 oligonucleotides** — 7 backbone classes, 4 sugar modifications, systematic GC/CpG/gapmer gradients
- **47 toxicity endpoints** across 6 mechanisms (hepatotoxicity, nephrotoxicity, immunotoxicity, complement, coagulopathy, thrombocytopenia)
- **5 human in vitro model systems** — primary hepatocytes, kidney organoids, PBMCs, platelets, organ-on-chip
- **3 open-source ML models** — OligoTox-XGB, OligoTox-Transformer, OligoTox-ActiveLearn
- **FAIR-compliant** — CC BY 4.0 data, Apache 2.0 code, DOIs via Zenodo

---

## Quick Start

```bash
# Install
pip install -e ".[ml,portal]"

# Generate synthetic data for testing
python scripts/generate_synthetic_data.py --n-oligos 200 --output data/test/

# Compute sequence features
oligotox-features data/test/compounds.csv data/test/features.csv

# Run QC pipeline
oligotox-qc data/test/results.csv data/test/plate_controls.csv

# Validate a submission against OTMRS standard
python -m oligotoxdb.otmrs_validator data/test/compounds.csv data/test/results.csv

# Launch the web portal
oligotox-portal

# Run tests
python -m pytest tests/ -v
```

## Docker

```bash
docker build -t oligotoxdb:latest .

# Run tests
docker run oligotoxdb

# Launch portal
docker run -p 7860:7860 -e OLIGOTOXDB_PATH=/data/db.duckdb oligotoxdb python portal/app.py
```

---

## Project Structure

```
oligotoxdb/
├── oligotoxdb/
│   ├── endpoints.py         # All 47 toxicity endpoints (canonical registry)
│   ├── features.py          # Sequence & physicochemical feature computation
│   ├── qc.py                # Plate QC, 4PL dose-response fitting, Grubbs test
│   ├── database.py          # DuckDB interface (OTMRS schema)
│   ├── ingestion.py         # ETL pipeline with Pydantic validation
│   ├── omics.py             # RNA-seq and proteomics integration
│   └── otmrs_validator.py   # OTMRS compliance validator for external submissions
├── models/
│   ├── xgb_model.py         # OligoTox-XGB (multi-endpoint, SHAP interpretable)
│   ├── transformer_model.py # OligoTox-Transformer (Nucleotide Transformer fine-tune)
│   └── active_learning.py   # OligoTox-ActiveLearn (GP Bayesian active learning)
├── benchmarks/
│   ├── create_splits.py     # 3 standardized benchmark splits
│   └── evaluate.py          # RMSE, R², Spearman ρ, AUROC, F1 evaluation
├── portal/
│   └── app.py               # Gradio web app (predictor, explorer, DR viewer, PCA)
├── scripts/
│   ├── generate_synthetic_data.py   # Synthetic data for CI/testing
│   └── release_batch.py             # Rolling batch release (Zenodo + HuggingFace)
├── notebooks/
│   └── OligoToxDB_Analysis_Walkthrough.ipynb
├── data/
│   └── controls.csv         # Positive/negative control compound definitions
├── tests/                   # 57 unit tests (pytest)
├── Dockerfile
└── pyproject.toml
```

---

## Models

### OligoTox-XGB
Gradient-boosted multi-endpoint regressor + classifier per toxicity endpoint. SHAP values computed for every prediction.

```python
from models.xgb_model import OligoToxXGB

model = OligoToxXGB.load("models/xgb/")
predictions = model.predict(features_df)
shap_df = model.explain(features_df, endpoint="Cell_viability_ATPLite")
```

### OligoTox-Transformer
Fine-tuned [Nucleotide Transformer](https://github.com/instadeepai/nucleotide-transformer) (500M params) with chemical modification cross-attention and MC-Dropout uncertainty.

```python
from models.transformer_model import OligoToxTransformer
# See notebooks/OligoToxDB_Analysis_Walkthrough.ipynb
```

### OligoTox-ActiveLearn
Gaussian Process surrogate with multi-objective Expected Information Gain acquisition for experimental design. Selects which compounds to synthesize next to maximize information gain.

```python
from models.active_learning import select_next_batch

result = select_next_batch(
    train_df, candidate_df,
    feature_cols=feature_cols, endpoints=endpoints,
    batch_size=200, alpha=0.5, round_number=1,
)
print(result.selected_ids)
```

---

## Benchmarks

Three standardized benchmark splits for community model comparison:

| Benchmark | Split Type | Task | Primary Metric |
|---|---|---|---|
| OligoTox-RandBench | Random 80/10/10 | IC50 regression | Spearman ρ |
| OligoTox-ChemSplit | By backbone class | Cross-chemistry generalization | AUROC |
| OligoTox-ProspBench | Batch A-C → D-E | Prospective prediction | RMSE |

Community leaderboard: [GitHub Discussions](https://github.com/anote-ai/nih-oligotox/discussions)

---

## Data Release Schedule

| Batch | Compounds | Release Date | Zenodo |
|---|---|---|---|
| A | ~560 | Aug 15, 2026 | pending |
| B | ~560 | Sep 15, 2026 | pending |
| C | ~560 | Oct 15, 2026 | pending |
| D (active learning) | ~560 | Nov 1, 2026 | pending |
| E (active learning) | ~560 | Nov 30, 2026 | pending |

---

## External Contributions (OTMRS)

External labs can contribute data following the **OligoTox Minimum Reporting Standard (OTMRS)**:

```bash
# Validate your submission before contributing
python -m oligotoxdb.otmrs_validator your_compounds.csv your_results.csv --output report.html
```

See [Appendix C of the Phase 2 submission](../OligoTox_Phase2_Submission_Anote.md) for the full OTMRS specification.

---

## Citation

```bibtex
@dataset{oligotoxdb2026,
  author    = {Vidra, Natan and Anote, Inc.},
  title     = {OligoToxDB: Large-scale human in vitro oligonucleotide toxicity dataset},
  year      = {2026},
  publisher = {Zenodo},
  license   = {CC BY 4.0},
  url       = {https://github.com/anote-ai/nih-oligotox}
}
```

## License

- **Data**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — use freely, attribution required
- **Code**: [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)

## Contact

oligotoxdb@anote.ai | [GitHub Issues](https://github.com/anote-ai/nih-oligotox/issues)
