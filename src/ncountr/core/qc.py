"""Quality control checks for Nanostring nCounter data."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy import stats

from ncountr.experiment import NanostringExperiment


def _fov_ratio(lane_info: pd.DataFrame) -> pd.Series:
    """Compute FOV counted / FOV total ratio per sample."""
    return lane_info["FovCounted"] / lane_info["FovCount"].replace(0, np.nan)


def _pos_ctrl_r2(pos_counts: pd.DataFrame) -> pd.Series:
    """Compute R-squared of log-log positive control linearity per sample.

    Expected concentrations are extracted from control names, e.g.
    ``POS_A(128)`` → 128.
    """
    conc_map: dict[str, float] = {}
    for name in pos_counts.index:
        match = re.search(r"\(([\d.]+)\)", name)
        if match:
            conc_map[name] = float(match.group(1))

    names_with_conc = [n for n in pos_counts.index if n in conc_map]
    concs = np.array([conc_map[n] for n in names_with_conc])

    r2_values: dict[str, float] = {}
    for sid in pos_counts.columns:
        counts = np.array([pos_counts.loc[n, sid] for n in names_with_conc])
        mask = (concs > 0) & (counts > 0)
        if mask.sum() >= 3:
            r, _p = stats.pearsonr(np.log10(concs[mask]), np.log10(counts[mask]))
            r2_values[sid] = r ** 2
        else:
            r2_values[sid] = np.nan
    return pd.Series(r2_values, name="PosCtrl_R2")


def _neg_background(neg_counts: pd.DataFrame, n_sd: float = 2.0) -> pd.Series:
    """Compute negative control background threshold (mean + n_sd * SD)."""
    bg: dict[str, float] = {}
    for sid in neg_counts.columns:
        vals = neg_counts[sid].values.astype(float)
        bg[sid] = vals.mean() + n_sd * vals.std()
    return pd.Series(bg, name="NegBackground")


def qc(
    experiment: NanostringExperiment,
    *,
    fov_ratio_threshold: float = 0.75,
    pos_r2_threshold: float = 0.95,
    neg_sd: float = 2.0,
) -> pd.DataFrame:
    """Run QC checks and store results on the experiment.

    Parameters
    ----------
    experiment : NanostringExperiment
    fov_ratio_threshold : float
        Minimum acceptable FOV ratio.
    pos_r2_threshold : float
        Minimum R-squared for positive control linearity.
    neg_sd : float
        Number of standard deviations above mean for negative background.

    Returns
    -------
    pd.DataFrame
        QC results indexed by sample.
    """
    qc_df = pd.DataFrame(index=experiment.samples)
    qc_df.index.name = "sample"

    # FOV ratio
    if not experiment.lane_info.empty:
        fov = _fov_ratio(experiment.lane_info)
        qc_df["FovRatio"] = fov.reindex(experiment.samples)
        qc_df["FovPass"] = qc_df["FovRatio"] > fov_ratio_threshold
    else:
        qc_df["FovRatio"] = np.nan
        qc_df["FovPass"] = True

    # Binding density
    if "BindingDensity" in experiment.lane_info.columns:
        qc_df["BindingDensity"] = experiment.lane_info["BindingDensity"].reindex(
            experiment.samples
        )

    # Cartridge
    if "CartridgeID" in experiment.lane_info.columns:
        qc_df["CartridgeID"] = experiment.lane_info["CartridgeID"].reindex(
            experiment.samples
        )

    # Positive control R2
    r2 = _pos_ctrl_r2(experiment.pos_counts)
    qc_df["PosCtrl_R2"] = r2.reindex(experiment.samples)
    qc_df["PosCtrlPass"] = qc_df["PosCtrl_R2"] > pos_r2_threshold

    # Negative background
    bg = _neg_background(experiment.neg_counts, n_sd=neg_sd)
    qc_df["NegBackground"] = bg.reindex(experiment.samples)

    experiment.qc_results = qc_df
    return qc_df
