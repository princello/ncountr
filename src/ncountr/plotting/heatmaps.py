"""Heatmap plotting utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_heatmap(
    data: pd.DataFrame,
    *,
    zscore: bool = True,
    vmin: float = -2,
    vmax: float = 2,
    cmap: str = "RdBu_r",
    title: str = "",
    xlabel_rotation: int = 45,
    ylabel_fontsize: int = 8,
    output: Union[str, Path, None] = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Plot a heatmap of expression data (genes x samples).

    Parameters
    ----------
    data : pd.DataFrame
        Genes (rows) x samples (columns).
    zscore : bool
        Z-score each row across columns.
    vmin, vmax : float
        Color scale limits.
    cmap : str
        Matplotlib colormap.
    title : str
    output : str or Path, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    mat = data.values.astype(float)

    if zscore:
        row_mean = np.nanmean(mat, axis=1, keepdims=True)
        row_std = np.nanstd(mat, axis=1, keepdims=True)
        row_std[row_std == 0] = np.nan
        mat = (mat - row_mean) / row_std

    n_genes, n_samples = mat.shape
    if figsize is None:
        figsize = (max(6, n_samples * 0.6), max(4, n_genes * 0.3))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(data.columns, rotation=xlabel_rotation, ha="right",
                        fontsize=max(6, 10 - n_samples // 5))
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(data.index, fontsize=ylabel_fontsize)

    plt.colorbar(im, ax=ax, shrink=0.6, label="Z-score" if zscore else "Expression")
    ax.set_title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()

    if output is not None:
        fig.savefig(output)
    return fig
