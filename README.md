# NIH NCATS OligoTox Open Data Challenge — Anote, Inc.

**Phase 1 Winner | Phase 2 Submission**

[![CI](https://github.com/anote-ai/nih-oligotox/actions/workflows/ci.yml/badge.svg)](https://github.com/anote-ai/nih-oligotox/actions)
[![License: Apache 2.0](https://img.shields.io/badge/Code-Apache%202.0-blue)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-green)](https://creativecommons.org/licenses/by/4.0/)
[![NIH NCATS OligoTox Challenge](https://img.shields.io/badge/NIH%20NCATS-OligoTox%20Phase%202-orange)](https://oligotox.com)

This repository contains the **Phase 2 submission** from [Anote, Inc.](https://anote.ai) for the [NIH NCATS Oligonucleotide Toxicity (OligoTox) Open Data Challenge](https://oligotox.com).
PI: Natan Vidra | Contact: oligotoxdb@anote.ai

---

## Submission Overview

Our submission — **OligoToxDB** — is a fully open dataset and AI prediction toolkit for oligonucleotide toxicity, built to directly address the challenge goals:

| Requirement | Our Approach |
|---|---|
| Human in vitro systems | 5 model systems: PHH, kidney organoids, PBMCs, platelets, organ-on-chip |
| Toxicity endpoints | 47 endpoints across hepatotoxicity, nephrotoxicity, immunotoxicity, complement, coagulopathy, thrombocytopenia |
| Oligo diversity | 2,800 oligos across 7 backbone classes, 4 sugar modifications, systematic GC/CpG/gapmer gradients |
| AI-ready dataset | DuckDB + Parquet, FAIR schema (OTMRS), DOIs via Zenodo |
| Predictive models | 3 open-source ML models with SHAP interpretability and uncertainty quantification |
| Open access | CC BY 4.0 (data), Apache 2.0 (code), rolling Zenodo releases Aug–Nov 2026 |

---

## Repository Contents

```
nih-oligotox/
├── OligoTox_Phase2_Submission_Anote.md   # Full Phase 2 narrative document
├── oligotoxdb/                            # OligoToxDB codebase
│   ├── oligotoxdb/        # Core library (endpoints, features, QC, database, ingestion, omics)
│   ├── models/            # OligoTox-XGB, OligoTox-Transformer, OligoTox-ActiveLearn
│   ├── benchmarks/        # Standardized splits + evaluation metrics
│   ├── portal/            # Gradio web portal (predictor, explorer, dose-response viewer)
│   ├── scripts/           # Synthetic data generation, batch release pipeline
│   ├── notebooks/         # Analysis walkthrough notebook
│   ├── tests/             # 57 unit tests
│   └── README.md          # OligoToxDB technical documentation
└── src/oligotox/          # Lightweight Python package (see below)
```

See **[`OligoTox_Phase2_Submission_Anote.md`](OligoTox_Phase2_Submission_Anote.md)** for the full combined narrative, methodology, OTMRS specification, and Public Access and Dissemination Plan (background/design reference).

See **[`phase2_submission/`](phase2_submission/)** for the actual submission package, split into
the four parts NCATS requires (narrative document, methodology document, PADP, dataset), with a
status tracker (`phase2_submission/README.md`) showing what's finalized vs. still pending real
experimental data.

See **[`oligotoxdb/README.md`](oligotoxdb/README.md)** for installation, quick-start, and technical documentation of the codebase.

---

## OligoToxDB at a Glance

- **2,800 oligonucleotides** spanning 7 backbone classes (PS, PO, PMO, LNA, PNA, 2′F-ANA, morpholino) and 4 sugar modifications
- **47 toxicity endpoints** measured across 5 human in vitro model systems
- **3 ML models**: gradient-boosted (XGB + SHAP), fine-tuned Nucleotide Transformer (500M), and Gaussian Process active learning
- **3 benchmark splits** for community model comparison (random, chemistry-stratified, prospective)
- **OTMRS standard** for external lab data contributions with automated validation
- **57 passing unit tests**, CI via GitHub Actions

### Quick Start

```bash
cd oligotoxdb
pip install -e ".[ml,portal]"

# Generate synthetic data (all 47 endpoints)
python scripts/generate_synthetic_data.py --n-oligos 200 --output data/test/ --tier 1

# Run tests
python -m pytest tests/ -v

# Launch interactive portal
oligotox-portal
```

---

## Data Release Schedule

Data generated under this challenge will be released publicly on Zenodo under CC BY 4.0:

| Batch | Compounds | Planned Release |
|---|---|---|
| A | ~560 | Aug 15, 2026 |
| B | ~560 | Sep 15, 2026 |
| C | ~560 | Oct 15, 2026 |
| D (active learning) | ~560 | Nov 1, 2026 |
| E (active learning) | ~560 | Nov 30, 2026 |

---

## About the Challenge

The [NIH NCATS OligoTox Open Data Challenge](https://oligotox.com) incentivizes publicly accessible, high-quality datasets from human cells to predict oligonucleotide toxicity from sequence and chemical modification. Anote, Inc. was selected as a **Phase 1 winner** (announced April 30, 2026) and is participating in Phase 2 (deadline December 31, 2026, prizes up to $400K).

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

---

## Python Package (`src/oligotox/`)

A lightweight Python package wrapping the core OligoTox ML pipeline for external use.

### Install

```bash
pip install -e ".[dev]"
```

### Quick Start

```python
from oligotox.data import make_dataset
from oligotox.core import extract_features, OligotoxPipeline
from oligotox.evaluate import aucroc

oligos, records = make_dataset(n=100)
features = [extract_features(o) for o in oligos]
labels = [int(r.value > 0.5) for r in records]

pipeline = OligotoxPipeline(model_type="xgb")
pipeline.fit(features[:80], [float(l) for l in labels[:80]])
scores = pipeline.predict_proba(features[80:])
print(f"AUC-ROC: {aucroc(scores, labels[80:]):.3f}")
```

### Run Tests

```bash
pytest tests/ -v --cov=src
```
