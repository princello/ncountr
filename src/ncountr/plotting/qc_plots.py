"""QC summary plots."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt

from ncountr.experiment import NanostringExperiment


def plot_qc(
    experiment: NanostringExperiment,
    *,
    output: Union[str, Path, None] = None,
    fov_threshold: float = 0.75,
) -> plt.Figure:
    """Generate a 4-panel QC summary figure.

    Panels: A) FOV ratio, B) Positive control linearity,
    C) Negative background, D) Housekeeping gene totals.

    Parameters
    ----------
    experiment : NanostringExperiment
    output : str or Path, optional
        Save figure to this path.
    fov_threshold : float
        FOV ratio threshold line.

    Returns
    -------
    matplotlib.figure.Figure
    """
    samples = experiment.samples
    qc = experiment.qc_results
    pos_df = experiment.pos_counts
    neg_df = experiment.neg_counts
    hk_df = experiment.hk_counts

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # A: FOV ratios
    ax = axes[0, 0]
    if qc is not None and "FovRatio" in qc.columns:
        fov_ratios = qc.loc[samples, "FovRatio"].values.astype(float)
        colors = ["#2ca02c" if r > fov_threshold else "#d62728" for r in fov_ratios]
        ax.bar(range(len(samples)), fov_ratios, color=colors)
        ax.axhline(fov_threshold, color="red", linestyle="--", alpha=0.7,
                    label=f"Threshold ({fov_threshold})")
        ax.set_ylim(min(0.7, np.nanmin(fov_ratios) - 0.02), 1.02)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_ylabel("FOV Counted / FOV Count")
    ax.set_title("A. Field of View Ratio")
    ax.legend(fontsize=8)

    # B: Positive control linearity
    ax = axes[0, 1]
    conc_map: dict[str, float] = {}
    for name in pos_df.index:
        match = re.search(r"\(([\d.]+)\)", name)
        if match:
            conc_map[name] = float(match.group(1))
    names_with_conc = [n for n in pos_df.index if n in conc_map]
    concs = np.array([conc_map[n] for n in names_with_conc])
    for sid in samples:
        counts = np.array([pos_df.loc[n, sid] for n in names_with_conc])
        mask = (concs > 0) & (counts > 0)
        ax.scatter(np.log10(concs[mask]), np.log10(counts[mask]),
                   s=20, alpha=0.6, label=sid)
    ax.set_xlabel("log10(Expected Concentration)")
    ax.set_ylabel("log10(Observed Count)")
    ax.set_title("B. Positive Control Linearity")
    ax.legend(fontsize=6, ncol=3, loc="lower right")

    # C: Negative control background
    ax = axes[1, 0]
    if qc is not None and "NegBackground" in qc.columns:
        bg_vals = qc.loc[samples, "NegBackground"].values.astype(float)
    else:
        bg_vals = [neg_df[s].values.astype(float).mean() + 2 * neg_df[s].values.astype(float).std()
                   for s in samples]
    ax.bar(range(len(samples)), bg_vals, color="steelblue")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_ylabel("Background (mean + 2*SD)")
    ax.set_title("C. Negative Control Background")

    # D: Housekeeping gene totals
    ax = axes[1, 1]
    hk_totals = hk_df[samples].sum(axis=0)
    ax.bar(range(len(samples)), hk_totals.values, color="steelblue")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_ylabel("Total HK Gene Counts")
    ax.set_title("D. Housekeeping Gene Totals")

    fig.suptitle("Nanostring QC Summary", fontsize=14, fontweight="bold")
    plt.tight_layout()

    if output is not None:
        fig.savefig(output)

    return fig
