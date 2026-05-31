"""
Sequence and physicochemical feature computation for oligonucleotides.

Each oligo is represented as a fixed-length numeric feature vector suitable
for ML model input. Features cover sequence composition, thermodynamics,
secondary structure, off-target binding, and chemical modification encoding.
"""

from __future__ import annotations

import re
import math
import itertools
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import pandas as pd
import click
from rich.console import Console
from rich.progress import track

console = Console()

# ─── Controlled vocabularies ────────────────────────────────────────────────

BACKBONE_CLASSES = ["PS", "PO", "PMO", "LNA_mix", "PNA", "2F_ANA", "morpholino"]
SUGAR_MODS = ["DNA", "2OMe", "2F", "LNA", "mixed"]
CONJUGATES = ["none", "GalNAc", "cholesterol", "lipid", "antibody"]
NUCLEOTIDES = list("ACGT")
DINUCLEOTIDES = [a + b for a, b in itertools.product("ACGT", repeat=2)]


@dataclass
class OligoFeatures:
    """All computed features for a single oligonucleotide."""

    oligo_id: str
    sequence: str

    # ── Sequence composition ──────────────────────────────────────────────
    length: int = 0
    gc_content: float = 0.0
    a_freq: float = 0.0
    c_freq: float = 0.0
    g_freq: float = 0.0
    t_freq: float = 0.0

    # Dinucleotide frequencies (16 values)
    cpg_count: int = 0
    cpg_freq: float = 0.0
    ga_freq: float = 0.0
    tc_freq: float = 0.0
    gg_freq: float = 0.0
    cc_freq: float = 0.0
    dinuc_entropy: float = 0.0  # Shannon entropy over dinucleotide distribution

    # ── Run / motif features ──────────────────────────────────────────────
    max_poly_a: int = 0
    max_poly_g: int = 0
    max_poly_c: int = 0
    max_poly_t: int = 0
    max_purine_run: int = 0
    max_pyrimidine_run: int = 0
    g4_score: float = 0.0        # G-quadruplex propensity (QGRS heuristic)
    self_comp_index: float = 0.0 # Fraction of bases that could form hairpin

    # ── Thermodynamic features ────────────────────────────────────────────
    mfe_estimate: float = 0.0    # Estimated MFE (kcal/mol); exact via RNAfold if available
    tm_estimate: float = 0.0     # Nearest-neighbor Tm estimate (°C)
    hybridization_dg: float = 0.0

    # ── Physicochemical ───────────────────────────────────────────────────
    molecular_weight: float = 0.0
    net_charge_ph7: int = 0      # Backbone-dependent
    hydrophobicity_index: float = 0.0

    # ── Modification encoding (one-hot) ───────────────────────────────────
    backbone_ps: int = 0
    backbone_po: int = 0
    backbone_pmo: int = 0
    backbone_lna_mix: int = 0
    backbone_pna: int = 0
    backbone_2f_ana: int = 0
    backbone_morpholino: int = 0

    sugar_dna: int = 0
    sugar_2ome: int = 0
    sugar_2f: int = 0
    sugar_lna: int = 0
    sugar_mixed: int = 0

    conjugate_none: int = 0
    conjugate_galnac: int = 0
    conjugate_cholesterol: int = 0
    conjugate_lipid: int = 0
    conjugate_antibody: int = 0

    # ── Gapmer architecture ───────────────────────────────────────────────
    is_gapmer: int = 0
    gap_length: int = 0
    wing_length_5p: int = 0
    wing_length_3p: int = 0

    # ── Off-target / specificity ──────────────────────────────────────────
    offtarget_score: float = 0.0   # Placeholder; populated by BLAST pipeline

    def to_vector(self) -> np.ndarray:
        """Return numeric feature vector (excludes id and sequence fields)."""
        d = asdict(self)
        excluded = {"oligo_id", "sequence"}
        return np.array([v for k, v in d.items() if k not in excluded], dtype=np.float32)

    @classmethod
    def feature_names(cls) -> list[str]:
        excluded = {"oligo_id", "sequence"}
        return [f for f in cls.__dataclass_fields__ if f not in excluded]


# ─── Feature computation ─────────────────────────────────────────────────────

def compute_features(
    oligo_id: str,
    sequence: str,
    backbone: str = "PS",
    sugar_mod: str = "DNA",
    conjugate: str = "none",
    gapmer_gap: int = 0,
    gapmer_wing_5p: int = 0,
    gapmer_wing_3p: int = 0,
) -> OligoFeatures:
    """Compute the full feature set for a single oligonucleotide."""
    seq = sequence.upper().replace("U", "T")
    n = len(seq)

    feat = OligoFeatures(oligo_id=oligo_id, sequence=seq)

    # ── Composition ──────────────────────────────────────────────────────
    feat.length = n
    counts = {nt: seq.count(nt) for nt in "ACGT"}
    feat.a_freq = counts["A"] / n
    feat.c_freq = counts["C"] / n
    feat.g_freq = counts["G"] / n
    feat.t_freq = counts["T"] / n
    feat.gc_content = (counts["G"] + counts["C"]) / n

    # Dinucleotide
    dinuc_counts = {d: 0 for d in DINUCLEOTIDES}
    for i in range(n - 1):
        di = seq[i : i + 2]
        if di in dinuc_counts:
            dinuc_counts[di] += 1
    total_di = sum(dinuc_counts.values()) or 1
    feat.cpg_count = dinuc_counts["CG"]
    feat.cpg_freq = dinuc_counts["CG"] / total_di
    feat.ga_freq = dinuc_counts["GA"] / total_di
    feat.tc_freq = dinuc_counts["TC"] / total_di
    feat.gg_freq = dinuc_counts["GG"] / total_di
    feat.cc_freq = dinuc_counts["CC"] / total_di

    di_freqs = np.array([v / total_di for v in dinuc_counts.values()])
    di_freqs = di_freqs[di_freqs > 0]
    feat.dinuc_entropy = float(-np.sum(di_freqs * np.log2(di_freqs)))

    # ── Runs ─────────────────────────────────────────────────────────────
    feat.max_poly_a = _max_run(seq, "A")
    feat.max_poly_g = _max_run(seq, "G")
    feat.max_poly_c = _max_run(seq, "C")
    feat.max_poly_t = _max_run(seq, "T")
    feat.max_purine_run = _max_run_set(seq, set("AG"))
    feat.max_pyrimidine_run = _max_run_set(seq, set("CT"))

    # G4 heuristic: count G-tracts of length ≥ 3
    g_tracts = re.findall(r"G{3,}", seq)
    feat.g4_score = float(len(g_tracts))

    # Self-complementarity: reverse complement overlap fraction
    rc = _reverse_complement(seq)
    feat.self_comp_index = _longest_common_substring_frac(seq, rc)

    # ── Thermodynamics ───────────────────────────────────────────────────
    feat.mfe_estimate = _estimate_mfe(seq)
    feat.tm_estimate = _nearest_neighbor_tm(seq)
    feat.hybridization_dg = feat.mfe_estimate * 0.85  # rough proxy

    # ── Physicochemical ───────────────────────────────────────────────────
    feat.molecular_weight = _molecular_weight(seq, backbone, sugar_mod)
    feat.net_charge_ph7 = _net_charge(n, backbone)
    feat.hydrophobicity_index = _hydrophobicity(backbone, conjugate)

    # ── One-hot encodings ─────────────────────────────────────────────────
    _set_onehot(feat, "backbone", backbone, BACKBONE_CLASSES)
    _set_onehot(feat, "sugar", sugar_mod, SUGAR_MODS)
    _set_onehot(feat, "conjugate", conjugate, CONJUGATES)

    # ── Gapmer ───────────────────────────────────────────────────────────
    feat.is_gapmer = 1 if gapmer_gap > 0 else 0
    feat.gap_length = gapmer_gap
    feat.wing_length_5p = gapmer_wing_5p
    feat.wing_length_3p = gapmer_wing_3p

    return feat


def compute_features_batch(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute features for a DataFrame of oligos.

    Expected columns: oligo_id, sequence, backbone, sugar_mod, conjugate,
                      gapmer_gap, gapmer_wing_5p, gapmer_wing_3p
    Returns: DataFrame with one row per oligo, columns = feature names.
    """
    rows = []
    for _, row in track(df.iterrows(), total=len(df), description="Computing features"):
        feat = compute_features(
            oligo_id=row["oligo_id"],
            sequence=row["sequence"],
            backbone=row.get("backbone", "PS"),
            sugar_mod=row.get("sugar_mod", "DNA"),
            conjugate=row.get("conjugate", "none"),
            gapmer_gap=int(row.get("gapmer_gap", 0)),
            gapmer_wing_5p=int(row.get("gapmer_wing_5p", 0)),
            gapmer_wing_3p=int(row.get("gapmer_wing_3p", 0)),
        )
        rows.append(asdict(feat))
    return pd.DataFrame(rows)


# ─── Helper functions ─────────────────────────────────────────────────────────

def _max_run(seq: str, nt: str) -> int:
    matches = re.findall(rf"{nt}+", seq)
    return max((len(m) for m in matches), default=0)


def _max_run_set(seq: str, nts: set) -> int:
    max_r = cur = 0
    for c in seq:
        if c in nts:
            cur += 1
            max_r = max(max_r, cur)
        else:
            cur = 0
    return max_r


_COMPLEMENT = str.maketrans("ACGT", "TGCA")

def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _longest_common_substring_frac(s1: str, s2: str) -> float:
    n = len(s1)
    if n == 0:
        return 0.0
    # dp approach, O(n²)
    dp = [[0] * (n + 1) for _ in range(n + 1)]
    best = 0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                best = max(best, dp[i][j])
    return best / n


def _estimate_mfe(seq: str) -> float:
    """
    Simplified MFE estimate based on GC content and G4 presence.
    For production, replace with subprocess call to RNAfold/Vienna RNA.
    MFE is negative (more stable = more negative).
    """
    n = len(seq)
    gc = sum(1 for c in seq if c in "GC") / n
    g4 = len(re.findall(r"G{3,}", seq))
    base_mfe = -gc * n * 1.5  # rough: GC pairs contribute ~-3 kcal/mol / 2
    g4_bonus = -g4 * 2.0
    return round(base_mfe + g4_bonus, 2)


def _nearest_neighbor_tm(seq: str, na_conc_mm: float = 50.0) -> float:
    """
    Simplified nearest-neighbor Tm for DNA:DNA duplex (SantaLucia 1998).
    For production, use BioPython's MeltingTemp module.
    Units: °C
    """
    # Nearest-neighbor parameters (ΔH kcal/mol, ΔS cal/mol/K) for DNA:DNA
    NN_PARAMS = {
        "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
        "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
        "GA": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
        "GG": (-8.0, -19.9), "AC": (-7.8, -21.0), "TC": (-8.2, -22.2),
        "AG": (-7.8, -21.0), "TG": (-8.5, -22.7), "TT": (-7.9, -22.2),
        "CC": (-8.0, -19.9),
    }
    R = 1.987  # cal/mol/K
    dH = dS = 0.0
    for i in range(len(seq) - 1):
        di = seq[i : i + 2]
        h, s = NN_PARAMS.get(di, (-8.0, -22.0))
        dH += h
        dS += s
    dH *= 1000  # to cal/mol
    dS_total = dS - R * math.log(na_conc_mm * 1e-3)
    if dS_total == 0:
        return 0.0
    tm_k = dH / dS_total
    return round(tm_k - 273.15, 1)


# Average nucleotide MW by backbone/sugar (approximate, g/mol per residue)
_NT_MW = {"A": 313.2, "C": 289.2, "G": 329.2, "T": 304.2}
_BACKBONE_MW_ADJUST = {"PS": 16.0, "PO": 0.0, "PMO": -62.0, "LNA_mix": 22.0, "PNA": -28.0}
_SUGAR_MW_ADJUST = {"DNA": 0.0, "2OMe": 14.0, "2F": 2.0, "LNA": 22.0}

def _molecular_weight(seq: str, backbone: str, sugar_mod: str) -> float:
    base_mw = sum(_NT_MW.get(nt, 310.0) for nt in seq)
    b_adj = _BACKBONE_MW_ADJUST.get(backbone, 0.0) * len(seq)
    s_adj = _SUGAR_MW_ADJUST.get(sugar_mod, 0.0) * len(seq)
    water = 18.02  # for terminal OH
    return round(base_mw + b_adj + s_adj + water, 1)


def _net_charge(length: int, backbone: str) -> int:
    if backbone in ("PS", "PO"):
        return -(length - 1)  # one negative charge per phosphodiester linkage
    if backbone in ("PMO", "PNA"):
        return 0  # neutral backbones
    return -(length - 1)


_HYDRO = {
    "none": 0.0, "GalNAc": -0.5, "cholesterol": 2.5, "lipid": 1.8, "antibody": -0.2
}
_BACKBONE_HYDRO = {"PS": 0.8, "PO": 0.0, "PMO": 0.3, "PNA": 0.5}

def _hydrophobicity(backbone: str, conjugate: str) -> float:
    return _BACKBONE_HYDRO.get(backbone, 0.0) + _HYDRO.get(conjugate, 0.0)


def _set_onehot(feat: OligoFeatures, prefix: str, value: str, vocab: list[str]) -> None:
    norm = value.lower().replace("-", "_").replace(" ", "_")
    for v in vocab:
        attr = f"{prefix}_{v.lower().replace('-', '_').replace(' ', '_')}"
        setattr(feat, attr, 1 if v.lower() == norm else 0)


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.command()
@click.argument("input_csv", type=click.Path(exists=True))
@click.argument("output_csv", type=click.Path())
def cli(input_csv: str, output_csv: str) -> None:
    """Compute features for all oligos in INPUT_CSV and write to OUTPUT_CSV."""
    df = pd.read_csv(input_csv)
    console.print(f"[bold]Computing features for {len(df)} oligos...[/bold]")
    feat_df = compute_features_batch(df)
    feat_df.to_csv(output_csv, index=False)
    console.print(f"[green]Saved {len(feat_df)} rows × {len(feat_df.columns)} columns → {output_csv}[/green]")


if __name__ == "__main__":
    cli()
