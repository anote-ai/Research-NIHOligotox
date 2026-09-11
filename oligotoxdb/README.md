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

# Generate synthetic data for testing (all 47 endpoints, Tier 1 = all 5 assay systems)
python scripts/generate_synthetic_data.py --n-oligos 200 --output data/test/ --tier 1

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
│   ├── features.py          # Sequence & physicochemical feature computation (45+ features)
│   ├── qc.py                # Plate QC (Z'-factor), 4PL dose-response fitting, Grubbs test
│   ├── database.py          # DuckDB interface (OTMRS schema, 7 tables)
│   ├── ingestion.py         # ETL pipeline with Pydantic v2 validation
│   ├── omics.py             # RNA-seq and proteomics integration (GEO/PRIDE export)
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
│   ├── generate_synthetic_data.py   # Synthetic data for CI/testing (all 47 endpoints)
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

## Toxicity Endpoint Panel

All 47 endpoints are defined in `oligotoxdb/endpoints.py` — the single source of truth used across the pipeline, models, and portal.

| Assay System | Count | Mechanism | Example Endpoints |
|---|---|---|---|
| Primary Human Hepatocytes (PHH) | 12 | Hepatotoxicity | Cell_viability_ATPLite, LDH_release, ALT_secretion, ROS_CellROX, Mitochondrial_JC1 |
| Kidney Proximal Tubule Organoids | 8 | Nephrotoxicity | KIM1_secretion, NGAL_secretion, Organoid_viability_3D, TightJunction_ZO1 |
| PBMC Co-culture | 18 | Immunotoxicity, Complement | IFNa, IL6, TNFa, C3a, C5a, TLR9_activation, NK_CD69 |
| Platelet-Rich Plasma | 7 | Coagulopathy, Thrombocytopenia | Platelet_aggregation, aPTT, PT, FactorXa_inhib |
| Liver-Kidney MPS (organ-on-chip) | 2 | Hepato/Nephrotoxicity | MPS_ALT, MPS_KIM1 |

---

## Models

### OligoTox-XGB
Gradient-boosted multi-endpoint regressor + calibrated classifier per toxicity endpoint. SHAP values for every prediction. One model per endpoint trained on log₁₀(IC50).

```python
from models.xgb_model import OligoToxXGB

model = OligoToxXGB.load("models/xgb/")
predictions = model.predict(features_df)           # IC50 + toxicity class per endpoint
shap_df = model.explain(features_df, endpoint="Cell_viability_ATPLite")
uncertainty = model.predict_with_uncertainty(features_df, n_samples=20)
```

### OligoTox-Transformer
Fine-tuned [Nucleotide Transformer](https://github.com/instadeepai/nucleotide-transformer) (500M params) with chemical modification cross-attention and MC-Dropout uncertainty quantification.

```python
from models.transformer_model import OligoToxTransformer
# See notebooks/OligoToxDB_Analysis_Walkthrough.ipynb for full training example
```

### OligoTox-ActiveLearn
Gaussian Process surrogate (Matérn ν=2.5) with multi-objective Expected Information Gain acquisition. Selects the next synthesis batch to maximize information gain across all endpoints.

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

Three standardized splits for community model comparison:

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

External labs can contribute data following the **OligoTox Minimum Reporting Standard (OTMRS)**. The validator checks required fields, controlled vocabularies, purity thresholds, replicate counts, and dose-response coverage.

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

nvidra@anote.ai | [GitHub Issues](https://github.com/anote-ai/nih-oligotox/issues)
