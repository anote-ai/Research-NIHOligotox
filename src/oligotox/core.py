"""Core data structures and feature engineering for OligoTox."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class BackboneClass(str, Enum):
    """Oligonucleotide backbone chemistry classifications."""

    PS = "PS"  # Phosphorothioate
    PO = "PO"  # Phosphodiester
    PMO = "PMO"  # Phosphorodiamidate Morpholino
    LNA = "LNA"  # Locked Nucleic Acid
    PNA = "PNA"  # Peptide Nucleic Acid
    F2_ANA = "F2_ANA"  # 2'-Fluoroarabino Nucleic Acid
    MORPHOLINO = "MORPHOLINO"  # Morpholino


class ToxicityEndpoint(str, Enum):
    """Toxicity endpoint categories assessed in the NIH OligoTox dataset."""

    HEPATOTOXICITY = "HEPATOTOXICITY"
    NEPHROTOXICITY = "NEPHROTOXICITY"
    IMMUNOTOXICITY = "IMMUNOTOXICITY"
    COMPLEMENT = "COMPLEMENT"
    COAGULOPATHY = "COAGULOPATHY"


@dataclass
class Oligonucleotide:
    """Represents a single oligonucleotide sequence with associated metadata."""

    oligo_id: str
    sequence: str
    backbone: BackboneClass
    modifications: list[str] = field(default_factory=list)
    gc_content: float = 0.0

    def __post_init__(self) -> None:
        if not self.sequence:
            raise ValueError("Sequence cannot be empty.")
        if self.gc_content == 0.0 and self.sequence:
            self.gc_content = compute_gc_content(self.sequence)


@dataclass
class ToxicityRecord:
    """A single toxicity measurement for an oligonucleotide."""

    oligo_id: str
    endpoint: ToxicityEndpoint
    value: float
    cell_system: str
    confidence: float = 1.0

    def __post_init__(self) -> None:
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"confidence must be in [0, 1], got {self.confidence}")


@dataclass
class FeatureVector:
    """Feature vector extracted from an Oligonucleotide for ML models."""

    oligo_id: str
    features: dict[str, float] = field(default_factory=dict)


def compute_gc_content(sequence: str) -> float:
    """Compute GC content as fraction of G and C bases in sequence.

    Args:
        sequence: Nucleotide sequence string (case-insensitive).

    Returns:
        GC fraction in [0, 1].

    Raises:
        ValueError: If sequence is empty.
    """
    if not sequence:
        raise ValueError("sequence cannot be empty")
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return gc_count / len(seq_upper)


def compute_cpg_ratio(sequence: str) -> float:
    """Compute CpG dinucleotide ratio relative to total dinucleotides.

    Args:
        sequence: Nucleotide sequence string (case-insensitive).

    Returns:
        CpG ratio in [0, 1], or 0.0 if fewer than 2 bases.
    """
    if len(sequence) < 2:
        return 0.0
    seq_upper = sequence.upper()
    total_dinuc = len(seq_upper) - 1
    cpg_count = sum(1 for i in range(total_dinuc) if seq_upper[i : i + 2] == "CG")
    return cpg_count / total_dinuc


def extract_features(oligo: Oligonucleotide) -> FeatureVector:
    """Extract a numeric feature vector from an Oligonucleotide.

    Features computed:
        - gc_content: fraction of G+C bases
        - cpg_ratio: CpG dinucleotide frequency
        - sequence_length: number of bases
        - has_ps: 1.0 if backbone is PS else 0.0
        - is_lna: 1.0 if backbone is LNA else 0.0
        - is_pmo: 1.0 if backbone is PMO else 0.0
        - modification_count: number of chemical modifications

    Args:
        oligo: Oligonucleotide instance.

    Returns:
        FeatureVector with 7 extracted numeric features.
    """
    features: dict[str, float] = {
        "gc_content": compute_gc_content(oligo.sequence),
        "cpg_ratio": compute_cpg_ratio(oligo.sequence),
        "sequence_length": float(len(oligo.sequence)),
        "has_ps": 1.0 if oligo.backbone == BackboneClass.PS else 0.0,
        "is_lna": 1.0 if oligo.backbone == BackboneClass.LNA else 0.0,
        "is_pmo": 1.0 if oligo.backbone == BackboneClass.PMO else 0.0,
        "modification_count": float(len(oligo.modifications)),
    }
    return FeatureVector(oligo_id=oligo.oligo_id, features=features)
