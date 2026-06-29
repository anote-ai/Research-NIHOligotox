"""Experiment 0 (Baseline): logistic-regression QSAR-style classifier on
synthetic OligoTox data, per DESIGN_DOC.md's "Baseline Experiment (Experiment 0)".

This script is a REAL, runnable implementation that produces a MEASURED
AUROC -- it does not reuse or restate the "expected results" tables in
DESIGN_DOC.md. It is intentionally scoped to what the current codebase can
actually support today:

  * Data: `oligotox.data.make_dataset()` synthetic generator (the project
    has no curated 3,200-point real dataset yet -- see DESIGN_DOC.md
    "Dataset Construction", which remains aspirational).
  * Features: the 7 features in `oligotox.core.extract_features()`
    (GC content, CpG ratio, length, is_ps, is_lna, n_modifications,
    motif_score). This is a much smaller feature set than the SeqChem
    representation described in DESIGN_DOC.md (k-mer embeddings +
    off-target features + ECFP4 + modification encoding); no chemistry
    SMILES/ECFP4 pipeline exists in this repo.
  * Model: sklearn LogisticRegression (the only model `OligotoxPipeline.fit`
    currently implements), used as a stand-in for the "QSAR (Random Forest
    + ECFP4)" baseline named in the design doc -- no Random Forest / ECFP4
    implementation exists yet either.
  * Split: random 80/20 (NOT the Murcko scaffold split described in the
    design doc -- no scaffold-splitting code exists in this repo, and
    oligonucleotides don't have a natural Murcko scaffold the way small
    molecules do, so this would need its own splitting logic, e.g. by
    backbone class or sequence similarity).

Because labels in `make_dataset()` are synthesized from simple,
near-deterministic rules over the same features used for prediction
(high GC -> elevated risk, PS backbone -> elevated risk, etc.), the
resulting AUROC mainly validates that the data/feature/model plumbing
works end-to-end -- it is NOT evidence for or against any of the
DESIGN_DOC.md hypotheses about real oligonucleotide toxicity. Treat the
printed number as a pipeline sanity check, not a scientific result.

Usage:
    python experiments/exp0_baseline.py --n 400 --seeds 10
"""
from __future__ import annotations

import argparse
import statistics
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oligotox.data import make_dataset  # noqa: E402
from oligotox.core import extract_features, OligotoxPipeline  # noqa: E402
from oligotox.evaluate import aucroc  # noqa: E402


def run_trial(n: int, seed: int, test_frac: float = 0.2) -> float:
    """Run one random-split train/test trial; return measured AUROC."""
    oligos, records = make_dataset(n=n, seed=seed)
    features = [extract_features(o) for o in oligos]
    labels = [int(r.value > 0.5) for r in records]

    split = int(n * (1 - test_frac))
    pipeline = OligotoxPipeline(model_type="logistic")
    pipeline.fit(features[:split], [float(l) for l in labels[:split]])

    scores = pipeline.predict_proba(features[split:])
    return aucroc(scores, labels[split:])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=400, help="Dataset size per trial")
    parser.add_argument("--seeds", type=int, default=10, help="Number of random seeds")
    args = parser.parse_args()

    aucs = [run_trial(args.n, seed) for seed in range(args.seeds)]
    mean_auc = statistics.mean(aucs)
    stdev_auc = statistics.pstdev(aucs)

    print("Experiment 0: Baseline (logistic-regression QSAR stand-in)")
    print(f"  n={args.n} synthetic oligos/trial, {args.seeds} seeds, random 80/20 split")
    print(f"  Per-seed AUROC: {[round(a, 4) for a in aucs]}")
    print(f"  Mean AUROC (measured): {mean_auc:.4f}")
    print(f"  Stdev across seeds:    {stdev_auc:.4f}")
    print()
    print("  NOTE: this is a MEASURED result on synthetic data with labels")
    print("  derived from the same features used for prediction. It confirms")
    print("  the pipeline runs end-to-end; it does NOT validate DESIGN_DOC.md's")
    print("  ~0.79 AUROC expectation on real hepatotoxicity data, which would")
    print("  require the curated 3,200-point dataset and ECFP4/RF baseline")
    print("  that do not yet exist in this repository.")


if __name__ == "__main__":
    main()
