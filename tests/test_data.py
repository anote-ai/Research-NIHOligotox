"""Tests for oligotox.data."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from oligotox.data import (
    generate_sequence, make_oligonucleotide, make_dataset, BACKBONE_CLASSES_LIST,
)
from oligotox.core import Oligonucleotide, compute_gc_content


def test_generate_sequence_length():
    seq = generate_sequence(length=30, seed=1)
    assert len(seq) == 30


def test_generate_sequence_gc_in_range():
    seq = generate_sequence(length=100, gc_target=0.5, seed=7)
    gc = compute_gc_content(seq)
    assert 0.3 <= gc <= 0.7


def test_make_oligonucleotide_type():
    oligo = make_oligonucleotide()
    assert isinstance(oligo, Oligonucleotide)


def test_make_dataset_sizes():
    oligos, records = make_dataset(n=10, seed=1)
    assert len(oligos) == 10
    assert len(records) == 10


def test_backbone_classes_list_length():
    assert len(BACKBONE_CLASSES_LIST) == 7
