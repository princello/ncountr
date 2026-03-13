"""Normalization methods for Nanostring nCounter data."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from ncountr.experiment import NanostringExperiment


def _geomean_scale(
    counts: pd.DataFrame, samples: list[str]
) -> dict[str, float]:
    """Compute geometric-mean scaling factors across samples.

    Each sample's geometric mean is computed, and the scaling factor brings
    each sample to the grand geometric mean.
    """
    geomeans: dict[str, float] = {}
    for sid in samples:
        vals = counts[sid].values.astype(float)
        vals = vals[vals > 0]
        if len(vals) == 0:
            geomeans[sid] = np.nan
        else:
            geomeans[sid] = sp_stats.gmean(vals)

    valid = [v for v in geomeans.values() if np.isfinite(v)]
    if not valid:
        return {sid: 1.0 for sid in samples}

    grand = sp_stats.gmean(valid)
    return {
        sid: grand / geomeans[sid] if np.isfinite(geomeans[sid]) else 1.0
        for sid in samples
    }


def _apply_scale(
    df: pd.DataFrame, scale: dict[str, float], samples: list[str]
) -> pd.DataFrame:
    """Multiply each sample column by its scaling factor."""
    result = df.copy()
    for sid in samples:
        result[sid] = df[sid] * scale[sid]
    return result


def normalize(
    experiment: NanostringExperiment,
    *,
    method: Literal["pos_only", "pos_hk", "pos_hk_bg"] = "pos_hk",
    neg_bg: pd.Series | dict[str, float] | None = None,
) -> pd.DataFrame:
    """Normalize raw counts and store the result on the experiment.

    Parameters
    ----------
    experiment : NanostringExperiment
    method : str
        ``"pos_only"`` — positive control normalization only.
        ``"pos_hk"`` — positive control + housekeeping normalization.
        ``"pos_hk_bg"`` — positive control + housekeeping + background subtraction.
    neg_bg : pd.Series or dict, optional
        Per-sample negative background values. Required for ``"pos_hk_bg"``.
        If not provided, computed from ``experiment.neg_counts``.

    Returns
    -------
    pd.DataFrame
        Normalized count matrix (genes x samples).
    """
    samples = experiment.samples

    # Step 1: positive control normalization
    pos_scale = _geomean_scale(experiment.pos_counts, samples)
    normalized = _apply_scale(experiment.raw_counts, pos_scale, samples)

    if method in ("pos_hk", "pos_hk_bg"):
        # Also apply pos scaling to housekeeping before computing HK factors
        hk_pos_scaled = _apply_scale(experiment.hk_counts, pos_scale, samples)
        hk_scale = _geomean_scale(hk_pos_scaled, samples)
        normalized = _apply_scale(normalized, hk_scale, samples)

    if method == "pos_hk_bg":
        if neg_bg is None:
            # Compute from negative controls
            bg: dict[str, float] = {}
            for sid in samples:
                vals = experiment.neg_counts[sid].values.astype(float)
                bg[sid] = vals.mean() + 2.0 * vals.std()
            neg_bg = bg

        if isinstance(neg_bg, pd.Series):
            neg_bg = neg_bg.to_dict()

        for sid in samples:
            threshold = neg_bg.get(sid, 0.0)
            normalized[sid] = (normalized[sid] - threshold).clip(lower=0)

    experiment.normalized = normalized
    return normalized


def get_scaling_factors(
    experiment: NanostringExperiment,
) -> dict[str, dict[str, float]]:
    """Compute and return scaling factors without modifying the experiment.

    Returns
    -------
    dict
        ``{"pos": {sid: factor}, "hk": {sid: factor}}``.
    """
    samples = experiment.samples
    pos_scale = _geomean_scale(experiment.pos_counts, samples)

    hk_pos_scaled = _apply_scale(experiment.hk_counts, pos_scale, samples)
    hk_scale = _geomean_scale(hk_pos_scaled, samples)

    return {"pos": pos_scale, "hk": hk_scale}
