"""Pathway / gene set scoring plots."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_pathway_scores(
    scores: pd.Series,
    groups: dict[str, list[str]],
    *,
    group_colors: dict[str, str] | None = None,
    output: Union[str, Path, None] = None,
    title: str = "Pathway Score",
    ylabel: str = "Pathway score (z-scored)",
) -> plt.Figure:
    """Box + strip plot of pathway scores by group.

    Parameters
    ----------
    scores : pd.Series
        Per-sample scores (index = sample ID).
    groups : dict[str, list[str]]
        Group name → list of sample IDs.
    group_colors : dict, optional
        Group name → color.
    output : str or Path, optional
        Save figure to this path.
    title, ylabel : str
        Plot labels.

    Returns
    -------
    matplotlib.figure.Figure
    """
    default_colors = {"group_a": "#d62728", "group_b": "#1f77b4"}
    if group_colors is None:
        group_colors = {}
        for i, name in enumerate(groups):
            group_colors[name] = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"][i % 4]

    fig, ax = plt.subplots(figsize=(6, 5))

    group_names = list(groups.keys())
    group_scores = [scores.reindex(groups[g]).dropna().values for g in group_names]

    bp = ax.boxplot(group_scores, labels=group_names, widths=0.5, patch_artist=True)
    for i, name in enumerate(group_names):
        color = group_colors.get(name, "#888888")
        bp["boxes"][i].set_facecolor(color)
        bp["boxes"][i].set_alpha(0.5)
        vals = group_scores[i]
        ax.scatter([i + 1] * len(vals), vals, c="black", s=40, zorder=5)

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()

    if output is not None:
        fig.savefig(output)
    return fig
