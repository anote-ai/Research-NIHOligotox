"""
OligoTox-ActiveLearn: Bayesian active learning for experimental design.

Strategy: Given the current dataset, select the next batch of compounds to
synthesize and screen by maximizing Expected Information Gain (EIG) across
a candidate pool of virtual oligos.

Architecture:
  - Surrogate model: Gaussian Process with deep RBF kernel (scikit-learn)
  - Acquisition functions: Expected Improvement (EI), Upper Confidence Bound (UCB),
    and multi-objective EIG (primary strategy)
  - Diversity constraint: ensures selected batch spans chemical space

This module is run between experimental batches (after Batch C, after Batch D)
to select which oligos to synthesize next.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.spatial.distance import cdist
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import click
from rich.console import Console
from rich.table import Table

console = Console()


# ─── Acquisition functions ────────────────────────────────────────────────────

def expected_improvement(
    mu: np.ndarray,
    sigma: np.ndarray,
    y_best: float,
    xi: float = 0.01,
) -> np.ndarray:
    """
    Expected Improvement acquisition function.
    Maximized at points predicted to exceed the current best by at least xi.
    """
    z = (mu - y_best - xi) / (sigma + 1e-9)
    return (mu - y_best - xi) * norm.cdf(z) + sigma * norm.pdf(z)


def upper_confidence_bound(
    mu: np.ndarray,
    sigma: np.ndarray,
    beta: float = 2.0,
) -> np.ndarray:
    """UCB: balance exploitation (high mu) and exploration (high sigma)."""
    return mu + beta * sigma


def multi_objective_eig(
    mu: np.ndarray,        # [N, n_endpoints]
    sigma: np.ndarray,     # [N, n_endpoints]
    y_best: np.ndarray,    # [n_endpoints] — current best per endpoint
    diversity_scores: np.ndarray,  # [N] — distance from existing data
    alpha: float = 0.5,    # weight: uncertainty vs. diversity
) -> np.ndarray:
    """
    Multi-objective Expected Information Gain.

    Combines:
      1. Mean EI across all endpoints (exploitation + exploration in objective space)
      2. Chemical space diversity (distance from existing data)

    Returns: [N] acquisition scores (higher = more valuable to screen)
    """
    # EI per endpoint
    ei_per_endpoint = np.stack([
        expected_improvement(mu[:, i], sigma[:, i], float(y_best[i]))
        for i in range(mu.shape[1])
    ], axis=1)  # [N, n_endpoints]

    mean_ei = ei_per_endpoint.mean(axis=1)  # [N]

    # Normalize both components to [0, 1]
    mean_ei_norm = (mean_ei - mean_ei.min()) / (mean_ei.max() - mean_ei.min() + 1e-9)
    div_norm = (diversity_scores - diversity_scores.min()) / (diversity_scores.max() - diversity_scores.min() + 1e-9)

    return alpha * mean_ei_norm + (1 - alpha) * div_norm


# ─── GP surrogate ─────────────────────────────────────────────────────────────

class OligoToxGPSurrogate:
    """
    Gaussian Process surrogate model for active learning.

    One GP per toxicity endpoint. Uses Matern kernel (ν=2.5) which is
    appropriate for biological data that may not be infinitely differentiable.
    """

    def __init__(self, endpoints: list[str], n_restarts: int = 5) -> None:
        self.endpoints = endpoints
        self.n_restarts = n_restarts
        self.gps: dict[str, GaussianProcessRegressor] = {}
        self.scalers: dict[str, StandardScaler] = {}
        self._fitted = False

    def fit(self, X: np.ndarray, y_dict: dict[str, np.ndarray]) -> None:
        """Fit one GP per endpoint on labeled training data."""
        for endpoint in self.endpoints:
            if endpoint not in y_dict:
                continue
            y = y_dict[endpoint]
            valid = ~np.isnan(y)
            if valid.sum() < 10:
                continue

            X_ep, y_ep = X[valid], y[valid]
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_ep)
            self.scalers[endpoint] = scaler

            kernel = Matern(nu=2.5) + WhiteKernel(noise_level=0.1)
            gp = GaussianProcessRegressor(
                kernel=kernel,
                n_restarts_optimizer=self.n_restarts,
                normalize_y=True,
                random_state=42,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gp.fit(X_scaled, y_ep)
            self.gps[endpoint] = gp

        self._fitted = True
        console.print(f"[green]GP fitted for {len(self.gps)} endpoints.[/green]")

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Predict mean and std for all endpoints.
        Returns: (mu [N, E], sigma [N, E])
        """
        n = len(X)
        n_ep = len(self.endpoints)
        mu = np.zeros((n, n_ep))
        sigma = np.zeros((n, n_ep))

        for i, endpoint in enumerate(self.endpoints):
            if endpoint not in self.gps:
                continue
            X_scaled = self.scalers[endpoint].transform(X)
            m, s = self.gps[endpoint].predict(X_scaled, return_std=True)
            mu[:, i] = m
            sigma[:, i] = s

        return mu, sigma


# ─── Active learning loop ────────────────────────────────────────────────────

@dataclass
class ActiveLearningResult:
    selected_ids: list[str]
    acquisition_scores: pd.DataFrame
    n_candidates: int
    n_selected: int
    round_number: int


def select_next_batch(
    train_df: pd.DataFrame,           # existing labeled data
    candidate_df: pd.DataFrame,       # unlabeled candidate pool
    feature_cols: list[str],
    endpoints: list[str],
    label_cols: Optional[list[str]] = None,
    batch_size: int = 200,
    alpha: float = 0.5,               # EIG weight: uncertainty vs diversity
    diversity_method: str = "maximin",  # "maximin" | "cluster"
    round_number: int = 1,
    random_seed: int = 42,
) -> ActiveLearningResult:
    """
    Select the next batch of compounds to synthesize and screen.

    Strategy:
      1. Fit GP surrogate on current training data
      2. Compute multi-objective EIG acquisition scores for all candidates
      3. Apply diversity constraint (avoid redundant selections)
      4. Return top-ranked candidates

    Args:
        train_df: Labeled data with feature columns + endpoint columns
        candidate_df: Unlabeled compounds (features only; oligo_id required)
        feature_cols: Feature column names
        endpoints: Toxicity endpoint names
        batch_size: Number of compounds to select
        alpha: Weight between uncertainty (1.0) and diversity (0.0)
        diversity_method: How to enforce batch diversity
    """
    np.random.seed(random_seed)

    X_train = train_df[feature_cols].fillna(0).values
    X_cand = candidate_df[feature_cols].fillna(0).values

    label_cols = label_cols or [f"log10_ic50_{ep}" for ep in endpoints]
    y_dict = {
        ep: train_df[lc].values
        for ep, lc in zip(endpoints, label_cols)
        if lc in train_df.columns
    }

    # Fit GP surrogate
    surrogate = OligoToxGPSurrogate(endpoints=endpoints)
    surrogate.fit(X_train, y_dict)

    # Predict on candidates
    mu, sigma = surrogate.predict(X_cand)

    # Current best per endpoint
    y_best = np.array([
        np.nanmax(y_dict.get(ep, np.array([0.0]))) for ep in endpoints
    ])

    # Diversity scores: mean distance from training set
    dists = cdist(X_cand, X_train, metric="euclidean")
    diversity_scores = dists.min(axis=1)  # distance to nearest training point

    # Multi-objective acquisition
    acq = multi_objective_eig(mu, sigma, y_best, diversity_scores, alpha=alpha)

    # Diversity-constrained selection (greedy maximin)
    if diversity_method == "maximin":
        selected_indices = _greedy_maximin_selection(acq, X_cand, batch_size)
    else:
        selected_indices = np.argsort(-acq)[:batch_size]

    selected_ids = candidate_df.iloc[selected_indices]["oligo_id"].tolist()

    # Build output table
    acq_df = pd.DataFrame({
        "oligo_id": candidate_df["oligo_id"].values,
        "acquisition_score": acq,
        "mean_uncertainty": sigma.mean(axis=1),
        "diversity_score": diversity_scores,
        "selected": False,
    })
    acq_df.loc[selected_indices, "selected"] = True
    acq_df = acq_df.sort_values("acquisition_score", ascending=False).reset_index(drop=True)

    console.print(f"\n[bold]Active Learning Round {round_number}[/bold]")
    console.print(f"  Candidates evaluated: {len(candidate_df)}")
    console.print(f"  Compounds selected: {len(selected_ids)}")
    console.print(f"  Mean acquisition score (selected): {acq[selected_indices].mean():.4f}")
    console.print(f"  Mean uncertainty (selected): {sigma[selected_indices].mean():.4f}")

    return ActiveLearningResult(
        selected_ids=selected_ids,
        acquisition_scores=acq_df,
        n_candidates=len(candidate_df),
        n_selected=len(selected_ids),
        round_number=round_number,
    )


def _greedy_maximin_selection(acq: np.ndarray, X: np.ndarray, k: int) -> np.ndarray:
    """
    Greedy maximin selection: maximize minimum pairwise distance within selected set,
    weighted by acquisition score.
    """
    n = len(X)
    selected = []

    # Start with the highest-acquisition compound
    first = int(np.argmax(acq))
    selected.append(first)

    dist_to_selected = cdist(X, X[[first]], metric="euclidean").squeeze()

    for _ in range(k - 1):
        remaining = [i for i in range(n) if i not in selected]
        if not remaining:
            break
        # Score = acquisition score * min distance to already-selected (normalized)
        acq_rem = acq[remaining]
        dist_rem = dist_to_selected[remaining]
        dist_norm = dist_rem / (dist_rem.max() + 1e-9)
        score = acq_rem * (0.5 + 0.5 * dist_norm)
        best = remaining[int(np.argmax(score))]
        selected.append(best)
        # Update distances
        new_dists = cdist(X, X[[best]], metric="euclidean").squeeze()
        dist_to_selected = np.minimum(dist_to_selected, new_dists)

    return np.array(selected)


# ─── Simulation: active learning vs random ───────────────────────────────────

def simulate_active_learning(
    full_df: pd.DataFrame,
    feature_cols: list[str],
    endpoints: list[str],
    label_cols: list[str],
    n_initial: int = 100,
    batch_size: int = 100,
    n_rounds: int = 5,
    n_trials: int = 3,
) -> pd.DataFrame:
    """
    Simulate active learning vs random selection to quantify efficiency gain.
    Returns DataFrame with model performance per round per strategy.
    """
    from sklearn.metrics import r2_score as sklearn_r2
    results = []

    for trial in range(n_trials):
        rng = np.random.default_rng(trial)
        idx = rng.permutation(len(full_df))
        initial_idx = idx[:n_initial]
        pool_idx = idx[n_initial:]

        for strategy in ["active", "random"]:
            labeled = list(initial_idx)
            pool = list(pool_idx)
            X_all = full_df[feature_cols].fillna(0).values

            for round_n in range(n_rounds):
                # Evaluate current model
                X_labeled = X_all[labeled]
                y_dict = {
                    ep: full_df[lc].values[labeled]
                    for ep, lc in zip(endpoints, label_cols)
                    if lc in full_df.columns
                }
                # Simple GP proxy for evaluation
                surrogate = OligoToxGPSurrogate(endpoints=[endpoints[0]])
                surrogate.fit(X_labeled, {endpoints[0]: y_dict.get(endpoints[0], np.zeros(len(labeled)))})
                X_test = X_all[pool[:200]]
                mu_test, _ = surrogate.predict(X_test)
                y_test = full_df[label_cols[0]].values[pool[:200]]
                valid = ~np.isnan(y_test)
                r2 = float(sklearn_r2(y_test[valid], mu_test[valid, 0])) if valid.sum() > 5 else np.nan

                results.append({
                    "trial": trial, "strategy": strategy, "round": round_n,
                    "n_labeled": len(labeled), "r2": r2
                })

                # Select next batch
                if strategy == "active" and len(pool) > 0:
                    cand_df = full_df.iloc[pool].assign(oligo_id=pool).reset_index(drop=True)
                    train_sub = full_df.iloc[labeled].reset_index(drop=True)
                    result = select_next_batch(
                        train_sub, cand_df, feature_cols, endpoints,
                        label_cols=label_cols, batch_size=batch_size, round_number=round_n + 1
                    )
                    new_idx = [int(oid) for oid in result.selected_ids if str(oid).isdigit()]
                else:
                    new_idx = pool[:batch_size]

                labeled.extend(new_idx[:batch_size])
                pool = [i for i in pool if i not in new_idx[:batch_size]]

    return pd.DataFrame(results)


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.group()
def cli() -> None:
    """OligoTox-ActiveLearn: Bayesian active learning for compound selection."""


@cli.command()
@click.argument("labeled_csv", type=click.Path(exists=True))
@click.argument("candidate_csv", type=click.Path(exists=True))
@click.option("--feature-cols-json", required=True)
@click.option("--endpoints-json", required=True)
@click.option("--batch-size", default=200, type=int)
@click.option("--alpha", default=0.5, type=float, help="EIG weight: uncertainty vs diversity")
@click.option("--round-number", default=1, type=int)
@click.option("--output", default="selected_compounds.csv")
def select(labeled_csv, candidate_csv, feature_cols_json, endpoints_json, batch_size, alpha, round_number, output):
    """Select next experimental batch from CANDIDATE_CSV using labeled data in LABELED_CSV."""
    labeled_df = pd.read_csv(labeled_csv)
    candidate_df = pd.read_csv(candidate_csv)
    with open(feature_cols_json) as f:
        feature_cols = json.load(f)
    with open(endpoints_json) as f:
        endpoints = json.load(f)

    result = select_next_batch(
        labeled_df, candidate_df, feature_cols, endpoints,
        batch_size=batch_size, alpha=alpha, round_number=round_number
    )

    result.acquisition_scores.to_csv(output, index=False)
    console.print(f"\n[green]Selected {result.n_selected} compounds saved to {output}[/green]")

    # Print top 10
    table = Table(title=f"Top 10 Selected Compounds (Round {round_number})")
    table.add_column("oligo_id")
    table.add_column("Acq. Score")
    table.add_column("Uncertainty")
    table.add_column("Diversity")
    for _, row in result.acquisition_scores[result.acquisition_scores["selected"]].head(10).iterrows():
        table.add_row(
            str(row["oligo_id"]),
            f"{row['acquisition_score']:.4f}",
            f"{row['mean_uncertainty']:.4f}",
            f"{row['diversity_score']:.4f}",
        )
    console.print(table)


if __name__ == "__main__":
    cli()
