"""
Plate-level and compound-level QC for OligoToxDB experimental data.

Implements:
  - Z'-factor per plate
  - 4PL dose-response curve fitting with IC50/Emax/Hill/NOEC/AUDRC
  - Grubbs outlier test (applied before dose-response fitting)
  - Inter-batch drift normalization
  - Full pipeline: assess_plates → flag bad plates → grubbs filter → fit curves
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist
import click
from rich.console import Console
from rich.table import Table

console = Console()

QC_MIN_Z_PRIME = 0.5
QC_MAX_CTRL_CV = 0.15
GRUBBS_ALPHA = 0.05


# ─── Plate QC ────────────────────────────────────────────────────────────────

@dataclass
class PlateQCResult:
    plate_id: str
    z_prime: float
    pos_ctrl_mean: float
    pos_ctrl_cv: float
    neg_ctrl_mean: float
    neg_ctrl_cv: float
    signal_window: float
    passed: bool
    failure_reasons: list[str] = field(default_factory=list)


def z_prime(pos_values: np.ndarray, neg_values: np.ndarray) -> float:
    """
    Z'-factor (Zhang et al. 1999).
    Z' = 1 - 3*(σ_pos + σ_neg) / |μ_pos - μ_neg|
    ≥ 0.5 is acceptable for HTS.
    """
    sigma_pos = np.std(pos_values, ddof=1)
    sigma_neg = np.std(neg_values, ddof=1)
    delta = abs(np.mean(pos_values) - np.mean(neg_values))
    if delta == 0:
        return -np.inf
    return 1.0 - 3.0 * (sigma_pos + sigma_neg) / delta


def assess_plate(
    plate_id: str,
    pos_ctrl_values: np.ndarray,
    neg_ctrl_values: np.ndarray,
) -> PlateQCResult:
    zp = z_prime(pos_ctrl_values, neg_ctrl_values)
    pos_mean = float(np.mean(pos_ctrl_values))
    neg_mean = float(np.mean(neg_ctrl_values))
    pos_cv = float(np.std(pos_ctrl_values, ddof=1) / pos_mean) if pos_mean != 0 else np.inf
    neg_cv = float(np.std(neg_ctrl_values, ddof=1) / neg_mean) if neg_mean != 0 else np.inf
    signal_window = (pos_mean - neg_mean) / neg_mean if neg_mean != 0 else 0.0

    failures: list[str] = []
    if zp < QC_MIN_Z_PRIME:
        failures.append(f"Z'={zp:.3f} < {QC_MIN_Z_PRIME}")
    if pos_cv > QC_MAX_CTRL_CV:
        failures.append(f"pos_ctrl CV={pos_cv:.2%} > {QC_MAX_CTRL_CV:.0%}")
    if neg_cv > QC_MAX_CTRL_CV:
        failures.append(f"neg_ctrl CV={neg_cv:.2%} > {QC_MAX_CTRL_CV:.0%}")

    return PlateQCResult(
        plate_id=plate_id,
        z_prime=round(zp, 4),
        pos_ctrl_mean=round(pos_mean, 2),
        pos_ctrl_cv=round(pos_cv, 4),
        neg_ctrl_mean=round(neg_mean, 2),
        neg_ctrl_cv=round(neg_cv, 4),
        signal_window=round(signal_window, 4),
        passed=len(failures) == 0,
        failure_reasons=failures,
    )


# ─── Grubbs outlier test ──────────────────────────────────────────────────────

def grubbs_test(values: np.ndarray, alpha: float = GRUBBS_ALPHA) -> np.ndarray:
    """
    Grubbs test for a single outlier (two-sided).
    Returns boolean mask: True = outlier.
    Minimum 3 values required; returns all-False for smaller arrays.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < 3:
        return np.zeros(n, dtype=bool)

    mean = np.mean(values)
    std = np.std(values, ddof=1)
    if std == 0:
        return np.zeros(n, dtype=bool)

    G = np.abs(values - mean) / std
    G_max = np.max(G)

    t_crit = t_dist.ppf(1 - alpha / (2 * n), n - 2)
    G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))

    outliers = np.zeros(n, dtype=bool)
    if G_max > G_crit:
        outliers[np.argmax(G)] = True
    return outliers


def flag_replicate_outliers(
    df: pd.DataFrame,
    value_col: str = "value_normalized",
    group_cols: tuple[str, ...] = ("oligo_id", "endpoint", "concentration_um"),
    alpha: float = GRUBBS_ALPHA,
) -> pd.DataFrame:
    """
    Apply Grubbs test across biological replicates for each (oligo, endpoint, concentration).
    Adds a boolean column 'grubbs_outlier' to the DataFrame.
    """
    df = df.copy()
    df["grubbs_outlier"] = False

    for _, grp in df.groupby(list(group_cols)):
        if len(grp) < 3:
            continue
        vals = grp[value_col].values
        outlier_mask = grubbs_test(vals, alpha=alpha)
        df.loc[grp.index[outlier_mask], "grubbs_outlier"] = True

    n_flagged = df["grubbs_outlier"].sum()
    if n_flagged > 0:
        console.print(f"[yellow]Grubbs test flagged {n_flagged} outlier replicates ({n_flagged/len(df):.1%})[/yellow]")
    return df


# ─── Dose-response fitting ────────────────────────────────────────────────────

@dataclass
class DoseResponseResult:
    oligo_id: str
    endpoint: str
    ic50: Optional[float]
    ic50_ci_low: Optional[float]
    ic50_ci_high: Optional[float]
    emax: Optional[float]
    hill: Optional[float]
    bottom: Optional[float]
    noec: Optional[float]
    audrc: float
    r_squared: float
    n_points: int
    n_outliers_removed: int
    fit_success: bool
    qc_flag: str


def _four_pl(x: np.ndarray, bottom: float, top: float, ic50: float, hill: float) -> np.ndarray:
    """4-parameter logistic (4PL) model."""
    return bottom + (top - bottom) / (1.0 + (ic50 / np.maximum(x, 1e-12)) ** hill)


def fit_dose_response(
    oligo_id: str,
    endpoint: str,
    concentrations: np.ndarray,
    responses: np.ndarray,
    min_valid_points: int = 4,
    remove_outliers: bool = True,
) -> DoseResponseResult:
    """
    Fit 4PL to dose-response data. Optionally removes Grubbs outliers first.
    """
    conc = np.asarray(concentrations, dtype=float)
    resp = np.asarray(responses, dtype=float)

    # Remove NaN
    valid = ~(np.isnan(conc) | np.isnan(resp))
    conc, resp = conc[valid], resp[valid]

    n_outliers = 0
    if remove_outliers and len(conc) >= 3:
        outlier_mask = grubbs_test(resp)
        n_outliers = int(outlier_mask.sum())
        conc, resp = conc[~outlier_mask], resp[~outlier_mask]

    if len(conc) < min_valid_points:
        return _failed_result(oligo_id, endpoint, conc, resp, n_outliers, "insufficient data points")

    sort_idx = np.argsort(conc)
    log_conc = np.log10(conc[sort_idx] + 1e-12)
    audrc = float(np.trapezoid(resp[sort_idx], log_conc)) if len(conc) > 1 else 0.0

    p0 = [min(resp), max(resp), float(np.median(conc)), 1.0]
    bounds = ([0, 0, conc.min() * 0.01, 0.1], [150, 150, conc.max() * 100, 10.0])

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(_four_pl, conc, resp, p0=p0, bounds=bounds, maxfev=5000)
        bottom, top, ic50, hill = popt

        resp_pred = _four_pl(conc, *popt)
        ss_res = np.sum((resp - resp_pred) ** 2)
        ss_tot = np.sum((resp - np.mean(resp)) ** 2)
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

        perr = np.sqrt(np.diag(np.abs(pcov)))  # abs to avoid sqrt of tiny negatives
        ic50_se = perr[2]
        df_fit = max(len(conc) - 4, 1)
        t_crit = t_dist.ppf(0.975, df_fit)
        ic50_low = max(0.0, ic50 - t_crit * ic50_se)
        ic50_high = ic50 + t_crit * ic50_se

        noec = _estimate_noec(conc, resp)
        qc = "pass" if r2 >= 0.85 else ("marginal" if r2 >= 0.65 else "fail")

        return DoseResponseResult(
            oligo_id=oligo_id, endpoint=endpoint,
            ic50=round(float(ic50), 4), ic50_ci_low=round(float(ic50_low), 4),
            ic50_ci_high=round(float(ic50_high), 4),
            emax=round(float(top), 2), hill=round(float(hill), 3),
            bottom=round(float(bottom), 2), noec=noec,
            audrc=round(audrc, 4), r_squared=round(r2, 4),
            n_points=int(len(conc)), n_outliers_removed=n_outliers,
            fit_success=True, qc_flag=qc,
        )
    except (RuntimeError, ValueError) as e:
        return _failed_result(oligo_id, endpoint, conc, resp, n_outliers, str(e))


def _estimate_noec(conc: np.ndarray, resp: np.ndarray, threshold: float = 20.0) -> Optional[float]:
    """Highest concentration with <threshold% deviation from vehicle (100%)."""
    sort_idx = np.argsort(conc)
    for c, r in zip(conc[sort_idx], resp[sort_idx]):
        if abs(r - 100.0) > threshold:
            return round(float(c), 4)
    return None


def _failed_result(
    oligo_id: str, endpoint: str,
    conc: np.ndarray, resp: np.ndarray,
    n_outliers: int, reason: str,
) -> DoseResponseResult:
    audrc = float(np.trapezoid(resp[np.argsort(conc)], np.log10(np.sort(conc) + 1e-12))) if len(conc) > 1 else 0.0
    return DoseResponseResult(
        oligo_id=oligo_id, endpoint=endpoint,
        ic50=None, ic50_ci_low=None, ic50_ci_high=None,
        emax=None, hill=None, bottom=None, noec=None,
        audrc=round(audrc, 4), r_squared=0.0,
        n_points=int(len(conc)), n_outliers_removed=n_outliers,
        fit_success=False, qc_flag="fail",
    )


# ─── Batch drift normalization ────────────────────────────────────────────────

def normalize_batch_drift(
    df: pd.DataFrame,
    value_col: str = "value_normalized",
    batch_col: str = "batch_id",
    reference_compound_col: str = "is_anchor",
) -> pd.DataFrame:
    """
    Median-ratio normalization using anchor compounds repeated across batches.
    Scales each batch so its anchor medians match the global median.
    """
    df = df.copy()
    anchors = df[df[reference_compound_col] == 1]
    if anchors.empty:
        console.print("[yellow]No anchor compounds found; skipping batch normalization.[/yellow]")
        return df

    global_median = anchors[value_col].median()
    for batch_id, batch_med in anchors.groupby(batch_col)[value_col].median().items():
        if batch_med == 0:
            continue
        df.loc[df[batch_col] == batch_id, value_col] *= global_median / batch_med

    return df


# ─── Full QC pipeline ────────────────────────────────────────────────────────

def run_qc_pipeline(
    raw_df: pd.DataFrame,
    plate_ctrl_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full QC pipeline:
      1. Plate-level Z'-factor check
      2. Flag data from failed plates
      3. Grubbs outlier detection across biological replicates
      4. Dose-response fitting (on Grubbs-cleaned, plate-QC-passed data)

    Args:
        raw_df: Experimental results with columns:
            oligo_id, endpoint, concentration_um, timepoint_h, batch_id,
            plate_id, replicate_bio, replicate_tech, value_raw, value_normalized
        plate_ctrl_df: Control measurements with columns:
            plate_id, ctrl_type (pos/neg), value_raw

    Returns:
        (plate_qc_df, dose_response_df)
    """
    # 1. Plate QC
    plate_qc_rows = []
    for plate_id, plate_data in plate_ctrl_df.groupby("plate_id"):
        pos = plate_data[plate_data["ctrl_type"] == "pos"]["value_raw"].values
        neg = plate_data[plate_data["ctrl_type"] == "neg"]["value_raw"].values
        if len(pos) >= 2 and len(neg) >= 2:
            result = assess_plate(str(plate_id), pos, neg)
            row = result.__dict__.copy()
            row["failure_reasons"] = "; ".join(result.failure_reasons)
            plate_qc_rows.append(row)
    plate_qc_df = pd.DataFrame(plate_qc_rows)

    # 2. Flag failed plates
    raw_df = raw_df.copy()
    if len(plate_qc_df) > 0:
        failed_plates = set(plate_qc_df.loc[~plate_qc_df["passed"], "plate_id"])
        raw_df["plate_qc_pass"] = ~raw_df["plate_id"].isin(failed_plates)
    else:
        raw_df["plate_qc_pass"] = True  # no control data → assume pass

    # 3. Grubbs outlier detection on plate-QC-passed data
    passed_df = raw_df[raw_df["plate_qc_pass"]].copy()
    if len(passed_df) > 0:
        passed_df = flag_replicate_outliers(passed_df)
        clean_df = passed_df[~passed_df["grubbs_outlier"]]
    else:
        clean_df = passed_df

    # 4. Dose-response fitting (mean across surviving replicates)
    dr_rows = []
    for (oligo_id, endpoint), grp in clean_df.groupby(["oligo_id", "endpoint"]):
        conc_resp = grp.groupby("concentration_um")["value_normalized"].mean()
        result = fit_dose_response(
            oligo_id=str(oligo_id),
            endpoint=str(endpoint),
            concentrations=conc_resp.index.values,
            responses=conc_resp.values,
            remove_outliers=False,  # already done above
        )
        dr_rows.append(result.__dict__)
    dr_df = pd.DataFrame(dr_rows)

    n_pass = int(plate_qc_df["passed"].sum()) if len(plate_qc_df) > 0 else 0
    n_total = len(plate_qc_df)
    n_dr_ok = int(dr_df["fit_success"].sum()) if len(dr_df) > 0 else 0
    console.print(
        f"QC complete | plates: {n_pass}/{n_total} passed | "
        f"DR fits: {n_dr_ok}/{len(dr_df)} successful"
    )
    return plate_qc_df, dr_df


# ─── CLI ─────────────────────────────────────────────────────────────────────

@click.command()
@click.argument("results_csv", type=click.Path(exists=True))
@click.argument("controls_csv", type=click.Path(exists=True))
@click.option("--out-qc", default="plate_qc.csv")
@click.option("--out-dr", default="dose_response.csv")
def cli(results_csv: str, controls_csv: str, out_qc: str, out_dr: str) -> None:
    """Run the full QC pipeline on RESULTS_CSV using control data in CONTROLS_CSV."""
    raw_df = pd.read_csv(results_csv)
    ctrl_df = pd.read_csv(controls_csv)
    console.print(f"[bold]QC: {len(raw_df)} measurements, {raw_df['plate_id'].nunique()} plates[/bold]")
    plate_qc, dr = run_qc_pipeline(raw_df, ctrl_df)
    plate_qc.to_csv(out_qc, index=False)
    dr.to_csv(out_dr, index=False)
    table = Table(title="QC Summary")
    table.add_column("Metric"); table.add_column("Value")
    table.add_row("Plates passed", f"{plate_qc['passed'].sum()}/{len(plate_qc)}")
    table.add_row("Mean Z'", f"{plate_qc['z_prime'].mean():.3f}")
    table.add_row("DR fits passed", f"{dr['fit_success'].sum()}/{len(dr)}")
    table.add_row("Mean R²", f"{dr.loc[dr['fit_success'], 'r_squared'].mean():.3f}")
    console.print(table)


if __name__ == "__main__":
    cli()
