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
    highlight_color: str = "gold",
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
    label_top_n: int = 15,
    output: Union[str, Path, None] = None,
    title: str | None = None,
) -> plt.Figure:
    """Generate a volcano plot from DE results.

    Significant genes (padj < threshold) are colored red (up) or blue (down).
    An optional set of genes of interest can be highlighted with colored
    markers on top, useful for visualizing pathway genes (e.g. IFN/JAK-STAT)
    or custom gene lists.

    Parameters
    ----------
    de_results : pd.DataFrame
        Must contain columns: ``gene``, ``log2FC``, ``pvalue``, ``padj``.
    highlight_genes : list[str], optional
        Genes to highlight with colored markers. Can be any gene list of
        interest (pathway genes, custom markers, etc.).
    highlight_label : str
        Legend label for highlighted genes.
    highlight_color : str
        Color for highlighted gene markers.
    padj_threshold : float
        Significance threshold for coloring.
    log2fc_threshold : float
        Minimum absolute log2FC to color significant genes (default 0).
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

    is_sig = de_results["padj"] < padj_threshold
    is_up = de_results["log2FC"] > log2fc_threshold
    is_down = de_results["log2FC"] < -log2fc_threshold

    colors = np.where(
        is_sig & is_up, "#d62728",
        np.where(is_sig & is_down, "#1f77b4", "#888888"),
    )
    sizes = np.where(is_sig, 25, 10)
    alphas = np.where(is_sig, 0.8, 0.4)

    # Plot non-significant first, then significant on top
    non_sig = ~is_sig
    if non_sig.any():
        ax.scatter(
            de_results.loc[non_sig, "log2FC"], neg_log10p[non_sig],
            c=colors[non_sig], s=sizes[non_sig], alpha=0.4, zorder=2,
        )
    if is_sig.any():
        ax.scatter(
            de_results.loc[is_sig, "log2FC"], neg_log10p[is_sig],
            c=colors[is_sig], s=sizes[is_sig], alpha=0.8, zorder=3,
            edgecolors="black", linewidth=0.3,
        )

    # Highlight gene set of interest
    if highlight_genes:
        mask = de_results["gene"].isin(highlight_genes)
        hl = de_results[mask]
        if len(hl) > 0:
            hl_nlp = -np.log10(hl["pvalue"].clip(lower=1e-10))
            ax.scatter(
                hl["log2FC"], hl_nlp, c=highlight_color, s=45,
                edgecolors="black", linewidth=0.5, zorder=5,
                label=f"{highlight_label} ({len(hl)} genes)",
            )
            # Label highlighted genes that are near-significant or significant
            for _, row in hl[hl["pvalue"] < 0.1].iterrows():
                ax.annotate(
                    row["gene"],
                    (row["log2FC"], -np.log10(max(row["pvalue"], 1e-10))),
                    fontsize=7, alpha=0.9, color="darkgoldenrod",
                    fontweight="bold",
                )

    # Label top genes by p-value
    labeled = set()
    if highlight_genes:
        labeled = set(de_results[de_results["gene"].isin(highlight_genes) & (de_results["pvalue"] < 0.1)]["gene"])
    for _, row in de_results.head(label_top_n).iterrows():
        if row["gene"] not in labeled:
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
        n_sig = is_sig.sum()
        ax.set_title(f"Volcano plot ({n_sig} significant at padj < {padj_threshold})")
    ax.legend(fontsize=9)

    plt.tight_layout()
    if output is not None:
        fig.savefig(output)
    return fig
