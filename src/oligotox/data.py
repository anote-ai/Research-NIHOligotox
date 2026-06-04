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


def make_dataset(
    n: int = 50, seed: int = 42
) -> Tuple[List[Oligonucleotide], List[ToxicityRecord]]:
    """Generate a paired dataset of oligonucleotides and toxicity records."""
    rng = random.Random(seed)
    backbones = list(BackboneClass)
    endpoints = list(ToxicityEndpoint)
    oligos = []
    records = []
    for i in range(n):
        oid = f"oligo_{i:04d}"
        backbone = rng.choice(backbones)
        oligo = make_oligonucleotide(oligo_id=oid, backbone=backbone, seed=seed + i)
        endpoint = rng.choice(endpoints)
        record = make_toxicity_record(oligo_id=oid, endpoint=endpoint, seed=seed + i)
        oligos.append(oligo)
        records.append(record)
    return oligos, records
