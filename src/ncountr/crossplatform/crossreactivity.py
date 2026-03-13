"""Cross-reactivity analysis using negative control samples."""

from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd

from ncountr.experiment import NanostringExperiment


def assess_crossreactivity(
    experiment: NanostringExperiment,
    *,
    neg_control_samples: list[str],
    reference_samples: list[str],
    counts: pd.DataFrame | None = None,
    thresholds: tuple[float, float] = (0.1, 0.5),
) -> pd.DataFrame:
    """Assess per-gene cross-reactivity using negative control samples.

    Cross-reactivity ratio = negative control signal / reference median signal.

    Parameters
    ----------
    experiment : NanostringExperiment
    neg_control_samples : list[str]
        Sample IDs that should have no genuine signal (e.g., NSG mouse tissue).
    reference_samples : list[str]
        Sample IDs with genuine signal for comparison.
    counts : pd.DataFrame, optional
        Count matrix to use.  Defaults to raw_counts.
    thresholds : tuple[float, float]
        Ratio thresholds for (reliable, moderate) classification.
        Default: <0.1 = reliable, <0.5 = moderate, >=0.5 = unreliable.

    Returns
    -------
    pd.DataFrame
        Per-gene cross-reactivity results with columns:
        neg_count, ref_median, xr_ratio, xr_class.
    """
    if counts is None:
        counts = experiment.raw_counts

    t_reliable, t_moderate = thresholds

    ref_median = counts[reference_samples].median(axis=1)

    results = pd.DataFrame(index=counts.index)
    for sid in neg_control_samples:
        neg_col = counts[sid]
        ratio = np.where(ref_median > 0, neg_col / ref_median, np.nan)
        results[f"neg_count_{sid}"] = neg_col
        results[f"xr_ratio_{sid}"] = ratio

    # Use first neg control for classification (or mean if multiple)
    if len(neg_control_samples) == 1:
        ratio_col = f"xr_ratio_{neg_control_samples[0]}"
    else:
        ratio_cols = [f"xr_ratio_{s}" for s in neg_control_samples]
        results["xr_ratio_mean"] = results[ratio_cols].mean(axis=1)
        ratio_col = "xr_ratio_mean"

    results["ref_median"] = ref_median

    def _classify(ratio: float) -> str:
        if pd.isna(ratio):
            return "no_signal"
        elif ratio < t_reliable:
            return "reliable"
        elif ratio < t_moderate:
            return "moderate"
        else:
            return "unreliable"

    results["xr_class"] = results[ratio_col].apply(_classify)

    return results
