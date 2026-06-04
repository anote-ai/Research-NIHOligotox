"""Synthetic data generation utilities for OligoTox."""
from __future__ import annotations

import random
from typing import List, Tuple

from .core import (
    BackboneClass,
    ToxicityEndpoint,
    Oligonucleotide,
    ToxicityRecord,
)

BACKBONE_CLASSES_LIST: List[str] = [b.value for b in BackboneClass]

_BASES = ["A", "T", "G", "C"]


def generate_sequence(length: int = 20, gc_target: float = 0.5, seed: int = 42) -> str:
    """Generate a DNA sequence with approximate GC content."""
    rng = random.Random(seed)
    seq = []
    for _ in range(length):
        if rng.random() < gc_target:
            seq.append(rng.choice(["G", "C"]))
        else:
            seq.append(rng.choice(["A", "T"]))
    return "".join(seq)


def make_oligonucleotide(
    oligo_id: str = "oligo_001",
    backbone: BackboneClass = BackboneClass.PS,
    seed: int = 42,
) -> Oligonucleotide:
    """Create a synthetic Oligonucleotide."""
    seq = generate_sequence(length=20, gc_target=0.5, seed=seed)
    return Oligonucleotide(
        oligo_id=oligo_id,
        sequence=seq,
        backbone=backbone,
        modifications=[],
    )


def make_toxicity_record(
    oligo_id: str,
    endpoint: ToxicityEndpoint = ToxicityEndpoint.HEPATOTOXICITY,
    seed: int = 42,
) -> ToxicityRecord:
    """Create a synthetic ToxicityRecord."""
    rng = random.Random(seed)
    return ToxicityRecord(
        oligo_id=oligo_id,
        endpoint=endpoint,
        value=round(rng.random(), 4),
        cell_system="PHH",
        confidence=round(0.5 + rng.random() * 0.5, 4),
    )


def _make_ps_hepatotox_oligo(oligo_id: str, seed: int) -> Tuple[Oligonucleotide, ToxicityRecord]:
    """PS backbone oligo with elevated hepatotoxicity (realistic pattern)."""
    rng = random.Random(seed)
    # PS + moderate-high GC -> hepatotoxicity risk 0.55-0.85
    seq = generate_sequence(length=20, gc_target=0.60, seed=seed)
    oligo = Oligonucleotide(oligo_id=oligo_id, sequence=seq, backbone=BackboneClass.PS)
    value = round(min(0.55 + rng.random() * 0.30, 1.0), 4)
    record = ToxicityRecord(
        oligo_id=oligo_id,
        endpoint=ToxicityEndpoint.HEPATOTOXICITY,
        value=value,
        cell_system="PHH",
        confidence=round(0.75 + rng.random() * 0.20, 4),
    )
    return oligo, record


def _make_high_gc_complement_oligo(oligo_id: str, seed: int) -> Tuple[Oligonucleotide, ToxicityRecord]:
    """High-GC oligo with complement activation risk."""
    rng = random.Random(seed)
    seq = generate_sequence(length=20, gc_target=0.75, seed=seed)
    backbone = rng.choice([BackboneClass.PS, BackboneClass.LNA])
    oligo = Oligonucleotide(oligo_id=oligo_id, sequence=seq, backbone=backbone)
    value = round(min(0.50 + rng.random() * 0.40, 1.0), 4)
    record = ToxicityRecord(
        oligo_id=oligo_id,
        endpoint=ToxicityEndpoint.COMPLEMENT,
        value=value,
        cell_system="serum",
        confidence=round(0.65 + rng.random() * 0.25, 4),
    )
    return oligo, record


def _make_low_risk_oligo(oligo_id: str, seed: int) -> Tuple[Oligonucleotide, ToxicityRecord]:
    """PO backbone, low GC oligo with low toxicity across endpoints."""
    rng = random.Random(seed)
    seq = generate_sequence(length=20, gc_target=0.35, seed=seed)
    oligo = Oligonucleotide(oligo_id=oligo_id, sequence=seq, backbone=BackboneClass.PO)
    endpoint = rng.choice(list(ToxicityEndpoint))
    value = round(rng.random() * 0.25, 4)  # low risk
    record = ToxicityRecord(
        oligo_id=oligo_id,
        endpoint=endpoint,
        value=value,
        cell_system="HEK293",
        confidence=round(0.60 + rng.random() * 0.30, 4),
    )
    return oligo, record


def make_dataset(
    n: int = 50, seed: int = 42
) -> Tuple[List[Oligonucleotide], List[ToxicityRecord]]:
    """Generate a paired dataset with realistic toxicity patterns.

    One-third of oligos are PS-backbone with elevated hepatotoxicity,
    one-third are high-GC with complement activation risk, and the
    remainder are low-risk PO-backbone oligos.
    """
    rng = random.Random(seed)
    oligos: List[Oligonucleotide] = []
    records: List[ToxicityRecord] = []

    n_hepato = n // 3
    n_complement = n // 3
    n_low = n - n_hepato - n_complement

    for i in range(n_hepato):
        oid = f"oligo_{i:04d}"
        oligo, record = _make_ps_hepatotox_oligo(oid, seed=rng.randint(0, 99999))
        oligos.append(oligo)
        records.append(record)

    for i in range(n_complement):
        oid = f"oligo_{n_hepato + i:04d}"
        oligo, record = _make_high_gc_complement_oligo(oid, seed=rng.randint(0, 99999))
        oligos.append(oligo)
        records.append(record)

    for i in range(n_low):
        oid = f"oligo_{n_hepato + n_complement + i:04d}"
        oligo, record = _make_low_risk_oligo(oid, seed=rng.randint(0, 99999))
        oligos.append(oligo)
        records.append(record)

    return oligos, records
