"""Differential expression plots."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_volcano(
    de_results: pd.DataFrame,
    *,
    highlight_genes: list[str] | None = None,
    highlight_label: str = "Highlighted",
    padj_threshold: float = 0.05,
    label_top_n: int = 15,
    output: Union[str, Path, None] = None,
    title: str | None = None,
) -> plt.Figure:
    """Generate a volcano plot from DE results.

    Parameters
    ----------
    de_results : pd.DataFrame
        Must contain columns: ``gene``, ``log2FC``, ``pvalue``, ``padj``.
    highlight_genes : list[str], optional
        Genes to highlight in gold.
    padj_threshold : float
        Significance threshold.
    label_top_n : int
        Number of top genes to label by p-value.
    output : str or Path, optional
        Save figure to this path.
    title : str, optional
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    neg_log10p = -np.log10(de_results["pvalue"].clip(lower=1e-10))

    colors = np.where(
        de_results["padj"] < padj_threshold,
        np.where(de_results["log2FC"] > 0, "#d62728", "#1f77b4"),
        "#888888",
    )
    ax.scatter(de_results["log2FC"], neg_log10p, c=colors, s=10, alpha=0.5)

    if highlight_genes:
        mask = de_results["gene"].isin(highlight_genes)
        hl = de_results[mask]
        hl_nlp = -np.log10(hl["pvalue"].clip(lower=1e-10))
        ax.scatter(
            hl["log2FC"], hl_nlp, c="gold", s=40,
            edgecolors="black", linewidth=0.5, zorder=5, label=highlight_label,
        )
        for _, row in hl[hl["padj"] < 0.1].iterrows():
            ax.annotate(
                row["gene"],
                (row["log2FC"], -np.log10(max(row["pvalue"], 1e-10))),
                fontsize=7, alpha=0.8, color="darkgoldenrod",
            )

    for _, row in de_results.head(label_top_n).iterrows():
        ax.annotate(
            row["gene"],
            (row["log2FC"], -np.log10(max(row["pvalue"], 1e-10))),
            fontsize=7, alpha=0.8,
        )

    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", alpha=0.5,
               label="p = 0.05")
    ax.axvline(0, color="grey", linestyle="-", alpha=0.3)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    if title:
        ax.set_title(title)
    else:
        n_sig = (de_results["padj"] < padj_threshold).sum()
        ax.set_title(f"Volcano plot ({n_sig} significant at padj < {padj_threshold})")
    ax.legend(fontsize=9)

    plt.tight_layout()
    if output is not None:
        fig.savefig(output)
    return fig
