"""Cross-platform correlation plots."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


def plot_correlation_scatter(
    x: pd.Series | np.ndarray,
    y: pd.Series | np.ndarray,
    *,
    labels: list[str] | None = None,
    xlabel: str = "Platform A",
    ylabel: str = "Platform B",
    title: str | None = None,
    highlight_idx: list[int] | None = None,
    output: Union[str, Path, None] = None,
) -> plt.Figure:
    """Scatter plot of two expression vectors with correlation annotation.

    Parameters
    ----------
    x, y : array-like
        Matched expression values.
    labels : list[str], optional
        Point labels (one per data point).
    xlabel, ylabel : str
    title : str, optional
    highlight_idx : list[int], optional
        Indices of points to highlight in red.
    output : str or Path, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x, y, s=8, alpha=0.4, c="grey")

    if highlight_idx is not None:
        ax.scatter(x[highlight_idx], y[highlight_idx], s=25, c="red", zorder=5)
        if labels is not None:
            for i in highlight_idx:
                ax.annotate(labels[i], (x[i], y[i]), fontsize=6, alpha=0.7)

    r, p = stats.spearmanr(x, y)
    if title is None:
        title = f"Spearman r = {r:.3f} (p = {p:.2e})"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()

    if output is not None:
        fig.savefig(output)
    return fig
