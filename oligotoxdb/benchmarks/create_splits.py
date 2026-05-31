"""
Create standardized benchmark splits for OligoToxDB.

Three splits for community model evaluation:
  1. OligoTox-RandBench: Random 80/10/10 train/val/test split
  2. OligoTox-ChemSplit: Split by backbone chemistry class (cross-chemistry generalization)
  3. OligoTox-ProspBench: Time/batch split (Batches A-C train → Batches D-E test)

Output format follows MoleculeNet convention: CSV files with split column.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import click
from rich.console import Console

console = Console()


def create_random_split(
    df: pd.DataFrame,
    label_col: str = "toxicity_class",
    train_frac: float = 0.8,
    val_frac: float = 0.1,
    seed: int = 42,
) -> pd.DataFrame:
    """Stratified random split on toxicity class."""
    df = df.copy()
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(df))
    n_train = int(len(df) * train_frac)
    n_val = int(len(df) * val_frac)

    df["split"] = "test"
    df.iloc[idx[:n_train], df.columns.get_loc("split")] = "train"
    df.iloc[idx[n_train : n_train + n_val], df.columns.get_loc("split")] = "val"
    return df


def create_chemistry_split(
    df: pd.DataFrame,
    backbone_col: str = "backbone_class",
    test_backbones: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Hold out specific backbone chemistry classes as test set.
    Default: PMO and PNA as test (least represented in historical data).
    """
    df = df.copy()
    if test_backbones is None:
        test_backbones = ["PMO", "PNA"]

    val_backbones = ["2F_ANA"]  # val: one more chemistry
    df["split"] = "train"
    df.loc[df[backbone_col].isin(val_backbones), "split"] = "val"
    df.loc[df[backbone_col].isin(test_backbones), "split"] = "test"
    return df


def create_prospective_split(
    df: pd.DataFrame,
    batch_col: str = "batch_group",
    train_batches: Optional[list[str]] = None,
    val_batches: Optional[list[str]] = None,
    test_batches: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Split by experimental batch to simulate prospective prediction.
    Train: Batches A-C | Val: Batch D | Test: Batch E
    """
    df = df.copy()
    train_batches = train_batches or ["A", "B", "C"]
    val_batches = val_batches or ["D"]
    test_batches = test_batches or ["E"]

    df["split"] = "train"
    df.loc[df[batch_col].isin(val_batches), "split"] = "val"
    df.loc[df[batch_col].isin(test_batches), "split"] = "test"
    return df


def generate_all_splits(
    df: pd.DataFrame,
    output_dir: Path,
    backbone_col: str = "backbone_class",
    batch_col: str = "batch_group",
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    splits = {
        "rand_bench": create_random_split(df),
        "chem_split": create_chemistry_split(df, backbone_col=backbone_col),
        "prosp_bench": create_prospective_split(df, batch_col=batch_col),
    }

    stats = []
    for name, split_df in splits.items():
        split_df.to_csv(output_dir / f"{name}.csv", index=False)
        dist = split_df["split"].value_counts()
        stats.append({
            "benchmark": name,
            "train": int(dist.get("train", 0)),
            "val": int(dist.get("val", 0)),
            "test": int(dist.get("test", 0)),
            "total": len(split_df),
        })
        console.print(f"[green]{name}[/green]: {dict(dist)}")

    pd.DataFrame(stats).to_csv(output_dir / "split_summary.csv", index=False)
    console.print(f"\nAll splits saved to {output_dir}")


@click.command()
@click.argument("dataset_csv", type=click.Path(exists=True))
@click.option("--output-dir", default="benchmarks/splits")
@click.option("--backbone-col", default="backbone_class")
@click.option("--batch-col", default="batch_group")
def cli(dataset_csv: str, output_dir: str, backbone_col: str, batch_col: str) -> None:
    """Generate all benchmark splits from DATASET_CSV."""
    df = pd.read_csv(dataset_csv)
    console.print(f"[bold]Generating benchmark splits for {len(df)} compounds...[/bold]")
    generate_all_splits(df, Path(output_dir), backbone_col=backbone_col, batch_col=batch_col)


if __name__ == "__main__":
    cli()
