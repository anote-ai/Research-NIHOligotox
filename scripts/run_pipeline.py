"""Demo script: generate data, fit OligotoxPipeline, evaluate."""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oligotox.data import make_dataset
from oligotox.core import extract_features, OligotoxPipeline
from oligotox.evaluate import aucroc, model_comparison, endpoint_breakdown


def main() -> None:
    print("Generating dataset...")
    oligos, records = make_dataset(n=50, seed=42)

    features = [extract_features(o) for o in oligos]
    labels = [int(r.value > 0.5) for r in records]

    split = 40
    pipeline = OligotoxPipeline(model_type="xgb")
    pipeline.fit(features[:split], [float(l) for l in labels[:split]])

    scores = pipeline.predict_proba(features[split:])
    auc = aucroc(scores, labels[split:])
    print(f"AUC-ROC (hold-out): {auc:.3f}")

    comparison = model_comparison(
        [{"model": "logistic", "aucroc": auc, "ece": 0.1}]
    )
    print(f"Model comparison: {comparison}")

    breakdown = endpoint_breakdown(records)
    print("Endpoint breakdown:")
    for ep, stats in breakdown.items():
        print(f"  {ep}: n={stats['n']}, mean_value={stats['mean_value']:.3f}")


if __name__ == "__main__":
    main()
