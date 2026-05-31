"""
Generate synthetic OligoToxDB data for development and testing.

Creates realistic but fake experimental data suitable for:
  - Testing the ingestion pipeline
  - Validating QC logic
  - Training initial ML models before real data is available
  - CI/CD pipeline verification

Usage:
    python scripts/generate_synthetic_data.py --n-oligos 200 --output data/synthetic/
"""

from __future__ import annotations

import random
import string
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

BACKBONE_CLASSES = ["PS", "PO", "PMO", "LNA_mix", "PNA", "2F_ANA", "morpholino"]
SUGAR_MODS = ["DNA", "2OMe", "2F", "LNA", "mixed"]
CONJUGATES = ["none", "none", "none", "GalNAc", "cholesterol", "lipid"]
ASSAY_SYSTEMS = ["PHH", "KidneyOrganoid", "PBMC", "Platelet", "MPS"]
ENDPOINTS = {
    "PHH": ["Cell_viability_ATPLite", "LDH_release", "ALT_secretion", "Caspase_3_7"],
    "KidneyOrganoid": ["KIM1_secretion", "NGAL_secretion"],
    "PBMC": ["IFNa", "IL6", "TNFa", "C3a", "C5a", "TLR9_activation"],
    "Platelet": ["Platelet_aggregation", "aPTT", "PT"],
    "MPS": ["Cell_viability_ATPLite", "KIM1_secretion"],
}
CONCENTRATIONS = [0.01, 0.1, 1.0, 10.0, 50.0, 100.0]
TIMEPOINTS = {
    "PHH": [24, 72], "KidneyOrganoid": [48, 96], "PBMC": [24],
    "Platelet": [2], "MPS": [72, 168],
}
VENDORS = ["IDT", "Eurofins"]


def random_sequence(length: int, cpg_prob: float = 0.05, g_rich: bool = False) -> str:
    nts = list("ACGT")
    weights = [0.25, 0.25, 0.25, 0.25]
    if g_rich:
        weights = [0.1, 0.1, 0.6, 0.2]
    seq = []
    for i in range(length):
        nt = random.choices(nts, weights=weights)[0]
        # Inject CpG
        if i > 0 and seq[-1] == "C" and random.random() < cpg_prob:
            nt = "G"
        seq.append(nt)
    return "".join(seq)


def generate_compound_library(n: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for i in range(n):
        length = int(rng.integers(8, 26))
        backbone = random.choice(BACKBONE_CLASSES)
        cpg_enriched = random.random() < 0.15
        g4_enriched = random.random() < 0.10
        seq = random_sequence(length, cpg_prob=0.3 if cpg_enriched else 0.05, g_rich=g4_enriched)

        is_gapmer = backbone in ("PS", "LNA_mix") and random.random() < 0.4
        gap = int(rng.integers(8, 13)) if is_gapmer else 0
        wing = int((length - gap) // 2) if is_gapmer else 0

        rows.append({
            "oligo_id": f"OT-{i+1:05d}",
            "sequence": seq,
            "backbone_class": backbone,
            "sugar_mod": random.choice(SUGAR_MODS),
            "conjugate": random.choice(CONJUGATES),
            "gapmer_design": f"{wing}-{gap}-{wing}" if is_gapmer else "",
            "gapmer_gap": gap,
            "gapmer_wing_5p": wing,
            "gapmer_wing_3p": wing,
            "synthesis_vendor": random.choice(VENDORS),
            "purity_percent": round(float(rng.uniform(92, 99.5)), 1),
            "mass_confirmed": True,
            "endotoxin_eu_ml": round(float(rng.exponential(0.1)), 3),
            "batch_group": ["A", "B", "C", "D", "E"][i % 5],
            "is_anchor": 1 if i < 20 else 0,
        })

    return pd.DataFrame(rows)


def _true_ic50(seq: str, backbone: str, endpoint: str, rng: np.random.Generator) -> float:
    """Simulate a biologically plausible IC50 based on sequence/chemistry rules."""
    cpg_count = seq.count("CG")
    poly_g = max((len(m) for m in __import__("re").findall(r"G+", seq)), default=0)
    gc = (seq.count("G") + seq.count("C")) / len(seq)

    base_ic50 = 20.0  # µM, default inactive

    if "IFN" in endpoint or "IL6" in endpoint or "TNF" in endpoint or "TLR9" in endpoint:
        if cpg_count >= 2:
            base_ic50 = float(rng.uniform(0.1, 5.0))
        elif backbone == "PS":
            base_ic50 = float(rng.uniform(5.0, 30.0))

    elif "C3a" in endpoint or "C5a" in endpoint or "Platelet" in endpoint:
        if poly_g >= 4:
            base_ic50 = float(rng.uniform(0.5, 10.0))
        elif backbone == "PS" and gc > 0.6:
            base_ic50 = float(rng.uniform(2.0, 20.0))

    elif "Cell_viability" in endpoint or "LDH" in endpoint or "ALT" in endpoint or "Caspase" in endpoint:
        if gc > 0.7 and backbone == "PS":
            base_ic50 = float(rng.uniform(1.0, 15.0))
        elif cpg_count >= 3:
            base_ic50 = float(rng.uniform(5.0, 50.0))

    elif "KIM1" in endpoint or "NGAL" in endpoint:
        if backbone == "PS" and len(seq) > 18:
            base_ic50 = float(rng.uniform(2.0, 25.0))

    # Add noise
    return base_ic50 * float(rng.lognormal(0, 0.3))


def generate_results(
    compounds_df: pd.DataFrame,
    rng: np.random.Generator,
    tier: str = "2",  # "1" = all systems, "2" = PHH+PBMC+Platelet, "3" = PHH+PBMC only
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate synthetic raw results and plate control data."""
    active_systems = {
        "1": ASSAY_SYSTEMS,
        "2": ["PHH", "PBMC", "Platelet"],
        "3": ["PHH", "PBMC"],
    }.get(tier, ["PHH"])

    result_rows = []
    plate_rows = []
    result_id = 1

    for _, compound in compounds_df.iterrows():
        for system in active_systems:
            for endpoint in ENDPOINTS[system]:
                ic50_true = _true_ic50(compound["sequence"], compound["backbone_class"], endpoint, rng)
                emax = float(rng.uniform(30, 100))
                hill = float(rng.uniform(0.8, 3.0))
                bottom = float(rng.uniform(0, 15))

                for conc in CONCENTRATIONS:
                    response_true = bottom + (emax - bottom) / (1 + (ic50_true / conc) ** hill)
                    plate_id = f"PLT-{compound['batch_group']}-{system[:3]}-{endpoint[:3]}-{int(conc*10):04d}"

                    for bio_rep in range(1, 4):  # 3 biological replicates
                        for tech_rep in range(1, 3):  # 2 technical replicates
                            noise = float(rng.normal(0, 3.0))
                            value_raw = max(0, response_true + noise) * float(rng.lognormal(0, 0.02))
                            result_rows.append({
                                "result_id": result_id,
                                "oligo_id": compound["oligo_id"],
                                "assay_system": system,
                                "endpoint": endpoint,
                                "concentration_um": conc,
                                "timepoint_h": TIMEPOINTS[system][0],
                                "replicate_bio": bio_rep,
                                "replicate_tech": tech_rep,
                                "batch_id": f"BATCH-{compound['batch_group']}",
                                "plate_id": plate_id,
                                "cell_donor_id": f"DONOR-{bio_rep:02d}",
                                "value_raw": round(value_raw, 3),
                                "value_normalized": round(value_raw, 3),
                                "unit": "% vehicle ctrl",
                                "plate_qc_pass": True,
                                "qc_flag": "pass",
                                "is_anchor": compound.get("is_anchor", 0),
                            })
                            result_id += 1

                # Generate plate controls for this plate/endpoint combo
                plate_id_ctrl = f"PLT-{compound['batch_group']}-{system[:3]}-{endpoint[:3]}-ctrl"
                for _ in range(8):
                    plate_rows.append({
                        "plate_id": plate_id_ctrl,
                        "ctrl_type": "neg",
                        "value_raw": round(float(rng.normal(100, 2.0)), 2),
                    })
                for _ in range(8):
                    plate_rows.append({
                        "plate_id": plate_id_ctrl,
                        "ctrl_type": "pos",
                        "value_raw": round(float(rng.normal(10, 0.8)), 2),
                    })

    return pd.DataFrame(result_rows), pd.DataFrame(plate_rows).drop_duplicates()


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic OligoToxDB data")
    parser.add_argument("--n-oligos", type=int, default=200)
    parser.add_argument("--output", default="data/synthetic")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tier", choices=["1", "2", "3"], default="2",
                        help="Assay coverage tier (1=all, 2=PHH+PBMC+Platelet, 3=PHH+PBMC)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Generating {args.n_oligos} compounds (Tier {args.tier} assays)...")
    compounds_df = generate_compound_library(args.n_oligos, rng)
    compounds_df.to_csv(out_dir / "compounds.csv", index=False)
    print(f"  Saved {len(compounds_df)} compounds → {out_dir / 'compounds.csv'}")

    print("Generating experimental results...")
    results_df, controls_df = generate_results(compounds_df, rng, tier=args.tier)
    results_df.to_csv(out_dir / "results.csv", index=False)
    controls_df.to_csv(out_dir / "plate_controls.csv", index=False)
    print(f"  Saved {len(results_df)} result rows → {out_dir / 'results.csv'}")
    print(f"  Saved {len(controls_df)} control rows → {out_dir / 'plate_controls.csv'}")
    print("\nDone. Use these files with:")
    print(f"  oligotox-ingest load-compounds {out_dir}/compounds.csv")
    print(f"  oligotox-ingest load-results {out_dir}/results.csv {out_dir}/plate_controls.csv")


if __name__ == "__main__":
    main()
