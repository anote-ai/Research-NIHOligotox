"""
OligoTox-Transformer: Fine-tune Nucleotide Transformer for multi-endpoint toxicity prediction.

Architecture:
  - Base: InstaDeepAI/nucleotide-transformer-v2-500m-multi-species (HuggingFace)
  - Chemical modification encoding: learned embedding fused via cross-attention
  - Multi-task head: 47 regression outputs + uncertainty (MC-Dropout)
  - Training: AdamW + cosine LR decay + gradient clipping

Interpretability:
  - Integrated Gradients per nucleotide position per endpoint
  - Attention rollout visualization

Requirements:
  pip install torch transformers datasets accelerate
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, AutoModel, get_cosine_schedule_with_warmup
import click
from rich.console import Console
from rich.progress import track

console = Console()

# ─── Constants ───────────────────────────────────────────────────────────────

BASE_MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
MAX_SEQ_LEN = 512  # NT tokenizer uses 6-mer tokens; 25-nt oligo → ~5 tokens
CHEMICAL_FEAT_DIM = 25   # number of chemical modification features
HIDDEN_DIM = 1024        # NT-v2-500M hidden size
MC_DROPOUT_SAMPLES = 50  # for uncertainty estimation


# ─── Dataset ─────────────────────────────────────────────────────────────────

class OligoToxDataset(Dataset):
    """
    PyTorch Dataset for OligoToxDB training.

    Each item: (sequence_tokens, chemical_features, labels, mask)
      - sequence_tokens: tokenized DNA sequence (6-mer tokens)
      - chemical_features: float32 tensor of modification/physicochemical features
      - labels: float32 tensor of log10(IC50) per endpoint (NaN → masked in loss)
      - mask: bool tensor, True where label exists
    """

    def __init__(
        self,
        df: pd.DataFrame,
        tokenizer: AutoTokenizer,
        endpoints: list[str],
        chem_feature_cols: list[str],
        max_length: int = MAX_SEQ_LEN,
    ) -> None:
        self.df = df.reset_index(drop=True)
        self.tokenizer = tokenizer
        self.endpoints = endpoints
        self.chem_feature_cols = chem_feature_cols
        self.max_length = max_length

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, idx: int) -> dict:
        row = self.df.iloc[idx]
        seq = row["sequence"].upper().replace("U", "T")

        # Tokenize sequence
        enc = self.tokenizer(
            seq,
            return_tensors="pt",
            padding="max_length",
            truncation=True,
            max_length=self.max_length,
        )

        # Chemical features
        chem = torch.tensor(
            [float(row.get(c, 0.0)) for c in self.chem_feature_cols],
            dtype=torch.float32,
        )

        # Labels
        label_vals = [float(row.get(ep, float("nan"))) for ep in self.endpoints]
        labels = torch.tensor(label_vals, dtype=torch.float32)
        mask = ~torch.isnan(labels)
        labels = torch.nan_to_num(labels, nan=0.0)

        return {
            "input_ids": enc["input_ids"].squeeze(0),
            "attention_mask": enc["attention_mask"].squeeze(0),
            "chem_features": chem,
            "labels": labels,
            "label_mask": mask,
        }


# ─── Model architecture ───────────────────────────────────────────────────────

class ChemicalModificationEncoder(nn.Module):
    """Projects chemical modification features into the same space as sequence embeddings."""

    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class MultiTaskToxicityHead(nn.Module):
    """
    Multi-task regression head: one output per toxicity endpoint.
    Includes MC-Dropout for uncertainty quantification.
    """

    def __init__(self, input_dim: int, n_endpoints: int, dropout: float = 0.2):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.LayerNorm(512),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        # Task-specific heads
        self.heads = nn.ModuleList([
            nn.Linear(512, 1) for _ in range(n_endpoints)
        ])
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        shared = self.shared(x)
        return torch.cat([head(self.dropout(shared)) for head in self.heads], dim=-1)


class OligoToxTransformer(nn.Module):
    """
    Full OligoTox-Transformer model.

    Forward pass:
      1. Encode DNA sequence via frozen/fine-tuned Nucleotide Transformer
      2. Pool sequence representation ([CLS] token)
      3. Encode chemical modification features via ChemicalModificationEncoder
      4. Fuse via learned cross-attention (sequence attends to chem context)
      5. Predict log10(IC50) for all endpoints via MultiTaskToxicityHead
    """

    def __init__(
        self,
        n_endpoints: int,
        chem_feature_dim: int = CHEMICAL_FEAT_DIM,
        dropout: float = 0.2,
        freeze_base: bool = False,
    ) -> None:
        super().__init__()
        self.n_endpoints = n_endpoints

        # Nucleotide Transformer backbone
        self.base = AutoModel.from_pretrained(BASE_MODEL_ID)
        if freeze_base:
            for param in self.base.parameters():
                param.requires_grad = False

        # Chemical modification encoder
        self.chem_encoder = ChemicalModificationEncoder(
            input_dim=chem_feature_dim,
            hidden_dim=256,
            output_dim=HIDDEN_DIM,
        )

        # Cross-attention: sequence [CLS] attends to chem context
        self.cross_attn = nn.MultiheadAttention(
            embed_dim=HIDDEN_DIM,
            num_heads=8,
            dropout=dropout,
            batch_first=True,
        )
        self.fusion_norm = nn.LayerNorm(HIDDEN_DIM)

        # Multi-task head
        self.head = MultiTaskToxicityHead(
            input_dim=HIDDEN_DIM,
            n_endpoints=n_endpoints,
            dropout=dropout,
        )

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        chem_features: torch.Tensor,
    ) -> torch.Tensor:
        # Sequence encoding
        outputs = self.base(input_ids=input_ids, attention_mask=attention_mask)
        cls_embed = outputs.last_hidden_state[:, 0, :]  # [B, H]

        # Chemical encoding
        chem_embed = self.chem_encoder(chem_features)  # [B, H]

        # Cross-attention fusion: CLS attends to chem embedding
        cls_q = cls_embed.unsqueeze(1)    # [B, 1, H]
        chem_k = chem_embed.unsqueeze(1)  # [B, 1, H]
        fused, _ = self.cross_attn(cls_q, chem_k, chem_k)
        fused = self.fusion_norm(fused.squeeze(1) + cls_embed)  # residual

        # Multi-task prediction
        return self.head(fused)  # [B, n_endpoints]

    def predict_with_uncertainty(
        self,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        chem_features: torch.Tensor,
        n_samples: int = MC_DROPOUT_SAMPLES,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """MC-Dropout uncertainty estimation. Returns (mean, std) over samples."""
        self.train()  # enable dropout
        preds = torch.stack([
            self(input_ids, attention_mask, chem_features)
            for _ in range(n_samples)
        ])  # [n_samples, B, n_endpoints]
        self.eval()
        return preds.mean(0), preds.std(0)


# ─── Training loop ────────────────────────────────────────────────────────────

@dataclass
class TrainingConfig:
    n_epochs: int = 20
    batch_size: int = 16
    learning_rate: float = 1e-4
    base_lr: float = 1e-5       # lower LR for pretrained backbone
    warmup_steps: int = 100
    weight_decay: float = 0.01
    grad_clip: float = 1.0
    dropout: float = 0.2
    freeze_base_epochs: int = 3  # freeze backbone for first N epochs
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    save_dir: str = "models/transformer"


def masked_mse_loss(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """MSE loss computed only over labeled (non-masked) entries."""
    loss = F.mse_loss(pred, target, reduction="none")
    if mask.sum() == 0:
        return torch.tensor(0.0, requires_grad=True)
    return (loss * mask.float()).sum() / mask.float().sum()


def train_transformer(
    train_df: pd.DataFrame,
    val_df: pd.DataFrame,
    endpoints: list[str],
    chem_feature_cols: list[str],
    config: TrainingConfig = TrainingConfig(),
) -> OligoToxTransformer:
    """Full training loop for OligoTox-Transformer."""
    device = torch.device(config.device)
    console.print(f"[bold]Training OligoTox-Transformer on {device}[/bold]")

    tokenizer = AutoTokenizer.from_pretrained(BASE_MODEL_ID)
    train_ds = OligoToxDataset(train_df, tokenizer, endpoints, chem_feature_cols)
    val_ds = OligoToxDataset(val_df, tokenizer, endpoints, chem_feature_cols)
    train_loader = DataLoader(train_ds, batch_size=config.batch_size, shuffle=True, num_workers=4)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size, shuffle=False, num_workers=4)

    model = OligoToxTransformer(
        n_endpoints=len(endpoints),
        chem_feature_dim=len(chem_feature_cols),
        dropout=config.dropout,
        freeze_base=True,  # start frozen
    ).to(device)

    # Optimizer with different LRs for backbone vs head
    optimizer = torch.optim.AdamW([
        {"params": model.base.parameters(), "lr": config.base_lr},
        {"params": list(model.chem_encoder.parameters()) +
                   list(model.cross_attn.parameters()) +
                   list(model.head.parameters()), "lr": config.learning_rate},
    ], weight_decay=config.weight_decay)

    total_steps = config.n_epochs * len(train_loader)
    scheduler = get_cosine_schedule_with_warmup(optimizer, config.warmup_steps, total_steps)

    best_val_loss = float("inf")
    save_dir = Path(config.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    for epoch in range(config.n_epochs):
        # Unfreeze backbone after freeze_base_epochs
        if epoch == config.freeze_base_epochs:
            for param in model.base.parameters():
                param.requires_grad = True
            console.print("[cyan]Backbone unfrozen.[/cyan]")

        # Train
        model.train()
        train_loss = 0.0
        for batch in train_loader:
            optimizer.zero_grad()
            preds = model(
                batch["input_ids"].to(device),
                batch["attention_mask"].to(device),
                batch["chem_features"].to(device),
            )
            loss = masked_mse_loss(preds, batch["labels"].to(device), batch["label_mask"].to(device))
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
            optimizer.step()
            scheduler.step()
            train_loss += loss.item()

        # Validate
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch in val_loader:
                preds = model(
                    batch["input_ids"].to(device),
                    batch["attention_mask"].to(device),
                    batch["chem_features"].to(device),
                )
                loss = masked_mse_loss(preds, batch["labels"].to(device), batch["label_mask"].to(device))
                val_loss += loss.item()

        avg_train = train_loss / len(train_loader)
        avg_val = val_loss / len(val_loader)
        console.print(f"Epoch {epoch+1}/{config.n_epochs} | train_loss={avg_train:.4f} | val_loss={avg_val:.4f}")

        if avg_val < best_val_loss:
            best_val_loss = avg_val
            torch.save(model.state_dict(), save_dir / "best_model.pt")
            console.print(f"  [green]New best model saved (val_loss={avg_val:.4f})[/green]")

    # Load best checkpoint
    model.load_state_dict(torch.load(save_dir / "best_model.pt", map_location=device))
    # Save config
    with open(save_dir / "config.json", "w") as f:
        json.dump({
            "endpoints": endpoints,
            "chem_feature_cols": chem_feature_cols,
            "n_endpoints": len(endpoints),
            "chem_feature_dim": len(chem_feature_cols),
        }, f, indent=2)

    return model


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """OligoTox-Transformer fine-tuning and inference."""


@cli.command()
@click.argument("train_csv", type=click.Path(exists=True))
@click.argument("val_csv", type=click.Path(exists=True))
@click.option("--endpoints-json", required=True, help="JSON file listing endpoint names")
@click.option("--chem-cols-json", required=True, help="JSON file listing chemical feature column names")
@click.option("--epochs", default=20, type=int)
@click.option("--batch-size", default=16, type=int)
@click.option("--save-dir", default="models/transformer")
def train(train_csv, val_csv, endpoints_json, chem_cols_json, epochs, batch_size, save_dir):
    """Fine-tune OligoTox-Transformer on TRAIN_CSV validated against VAL_CSV."""
    train_df = pd.read_csv(train_csv)
    val_df = pd.read_csv(val_csv)
    with open(endpoints_json) as f:
        endpoints = json.load(f)
    with open(chem_cols_json) as f:
        chem_cols = json.load(f)
    cfg = TrainingConfig(n_epochs=epochs, batch_size=batch_size, save_dir=save_dir)
    train_transformer(train_df, val_df, endpoints, chem_cols, cfg)
    console.print("[bold green]Training complete.[/bold green]")


if __name__ == "__main__":
    cli()
