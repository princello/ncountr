"""GSEA and pathway analysis plots."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_enrichment(
    ranked_genes: pd.Series,
    gene_set: list[str],
    running_scores: np.ndarray,
    *,
    es: float | None = None,
    pvalue: float | None = None,
    nes: float | None = None,
    gene_set_name: str = "",
    output: Union[str, Path, None] = None,
) -> plt.Figure:
    """Classic 3-panel GSEA enrichment plot.

    Parameters
    ----------
    ranked_genes : pd.Series
        Gene names as index, ranking metric as values (sorted descending).
    gene_set : list[str]
        Genes in the set.
    running_scores : np.ndarray
        Running enrichment score from ``_running_enrichment_score``.
    es : float, optional
        Enrichment score (annotated on plot).
    pvalue : float, optional
        Nominal p-value (annotated on plot).
    nes : float, optional
        Normalized enrichment score (annotated on plot).
    gene_set_name : str
        Title for the plot.
    output : str or Path, optional
        If given, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    genes = ranked_genes.index.tolist()
    metric_vals = ranked_genes.values
    n = len(genes)
    hit_mask = np.array([g in gene_set for g in genes], dtype=bool)
    hit_indices = np.where(hit_mask)[0]

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(8, 5),
        gridspec_kw={"height_ratios": [3, 0.5, 1.5], "hspace": 0.05},
        sharex=True,
    )

    # Panel 1: Running enrichment score
    ax1.plot(range(n), running_scores, color="#1a7a2e", linewidth=1.5)
    ax1.axhline(0, color="gray", linewidth=0.5, linestyle="--")

    # Mark peak
    if es is not None and es >= 0:
        peak = np.argmax(running_scores)
    elif es is not None:
        peak = np.argmin(running_scores)
    else:
        peak = np.argmax(np.abs(running_scores))
    ax1.axvline(peak, color="red", linewidth=0.5, linestyle=":", alpha=0.5)

    # Annotation
    info = []
    if es is not None:
        info.append(f"ES = {es:.3f}")
    if nes is not None:
        info.append(f"NES = {nes:.3f}")
    if pvalue is not None:
        info.append(f"p = {pvalue:.3f}")
    if info:
        ax1.text(
            0.98, 0.95 if (es or 0) >= 0 else 0.05,
            "\n".join(info),
            transform=ax1.transAxes, ha="right",
            va="top" if (es or 0) >= 0 else "bottom",
            fontsize=10, fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    title = gene_set_name or "GSEA Enrichment"
    ax1.set_title(title, fontsize=14, fontweight="bold")
    ax1.set_ylabel("Enrichment\nScore", fontsize=12)
    ax1.tick_params(labelbottom=False)

    # Panel 2: Barcode (hit positions)
    for idx in hit_indices:
        ax2.axvline(idx, color="black", linewidth=0.3)
    ax2.set_xlim(0, n)
    ax2.set_yticks([])
    ax2.set_ylabel("Hits", fontsize=10)
    ax2.tick_params(labelbottom=False)

    # Panel 3: Ranked metric
    ax3.fill_between(
        range(n), metric_vals, 0,
        where=metric_vals >= 0, color="#d62728", alpha=0.4,
    )
    ax3.fill_between(
        range(n), metric_vals, 0,
        where=metric_vals < 0, color="#2166ac", alpha=0.4,
    )
    ax3.axhline(0, color="gray", linewidth=0.5)
    ax3.set_xlabel("Gene Rank", fontsize=12)
    ax3.set_ylabel("Ranked\nMetric", fontsize=12)
    ax3.set_xlim(0, n)

    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_gsea_dotplot(
    gsea_results: pd.DataFrame,
    *,
    top_n: int = 20,
    significance_threshold: float = 0.25,
    output: Union[str, Path, None] = None,
    title: str = "GSEA Results",
) -> plt.Figure:
    """Dot plot of GSEA results.

    Parameters
    ----------
    gsea_results : pd.DataFrame
        Output from :func:`ncountr.gsea` or :func:`ncountr.competitive_test`.
        Must contain columns: gene_set, pvalue, and either nes or es.
    top_n : int
        Show top N gene sets by p-value.
    significance_threshold : float
        Gene sets with padj below this are highlighted.
    output : str or Path, optional
        If given, save figure.
    title : str
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    df = gsea_results.copy()

    # Use NES if available, else ES or mean_diff
    if "nes" in df.columns:
        x_col, x_label = "nes", "Normalized Enrichment Score"
    elif "es" in df.columns:
        x_col, x_label = "es", "Enrichment Score"
    elif "mean_diff" in df.columns:
        x_col, x_label = "mean_diff", "Mean Score Difference"
    elif "stat" in df.columns:
        # For competitive test: sign stat by direction
        if "direction" in df.columns:
            df["signed_stat"] = df["stat"] * df["direction"].map({"up": 1, "down": -1, "mixed": 0}).fillna(0)
            x_col, x_label = "signed_stat", "Signed Test Statistic"
        else:
            x_col, x_label = "stat", "Test Statistic"
    else:
        raise ValueError("GSEA results must contain nes, es, mean_diff, or stat column")

    df = df.sort_values("pvalue").head(top_n)
    df = df.sort_values(x_col)  # vertical order by score

    # Clean up names for display
    df["display_name"] = df["gene_set"].str.replace("HALLMARK_", "", regex=False)
    df["display_name"] = df["display_name"].str.replace("_", " ")

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.35)))

    neg_log10p = -np.log10(df["pvalue"].clip(lower=1e-10))
    sizes = neg_log10p * 80 + 20

    colors = np.where(df[x_col] > 0, "#d62728", "#2166ac")
    padj_col = "padj" if "padj" in df.columns else "pvalue"
    edge_colors = np.where(
        df[padj_col] < significance_threshold, "gold", "black"
    )
    edge_widths = np.where(
        df[padj_col] < significance_threshold, 2.0, 0.5
    )

    for i, (_, row) in enumerate(df.iterrows()):
        ax.scatter(
            row[x_col], i,
            s=sizes.iloc[i],
            c=colors[i],
            edgecolors=edge_colors[i],
            linewidths=edge_widths[i],
            alpha=0.7,
            zorder=3,
        )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["display_name"], fontsize=10)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="x", alpha=0.2)

    # Size legend
    from matplotlib.lines import Line2D
    for pval, label in [(0.3, "p=0.3"), (0.05, "p=0.05"), (0.005, "p=0.005")]:
        s = -np.log10(pval) * 80 + 20
        ax.scatter([], [], s=s, c="gray", alpha=0.5, label=label, edgecolors="black", linewidths=0.5)
    ax.scatter([], [], s=80, c="gray", edgecolors="gold", linewidths=2, label=f"padj<{significance_threshold}")
    ax.legend(loc="lower right", fontsize=10, title="Size = -log10(p)")

    # Overlap count annotation
    if "n_overlap" in df.columns:
        for i, (_, row) in enumerate(df.iterrows()):
            ax.text(
                ax.get_xlim()[1], i, f" n={row['n_overlap']}",
                va="center", fontsize=8, color="gray",
            )

    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_volcano_effect(
    de_results: pd.DataFrame,
    *,
    padj_threshold: float = 0.05,
    label_top_n: int = 15,
    highlight_genes: list[str] | None = None,
    highlight_label: str = "Highlighted",
    highlight_color: str = "gold",
    output: Union[str, Path, None] = None,
    title: str | None = None,
) -> plt.Figure:
    """Volcano plot with Cohen's d on the x-axis.

    Emphasizes effect size rather than statistical significance,
    which is more appropriate for small-n nCounter studies.

    Parameters
    ----------
    de_results : pd.DataFrame
        Output from ``de(effect_size=True)``.  Must contain columns:
        gene, cohens_d, pvalue, padj.
    padj_threshold : float
        Significance threshold for coloring.
    label_top_n : int
        Label top N genes by p-value.
    highlight_genes : list[str], optional
        Genes to highlight in a different color.
    highlight_label : str
        Label for highlighted genes.
    highlight_color : str
        Color for highlighted genes.
    output : str or Path, optional
    title : str, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    df = de_results.copy()
    if "cohens_d" not in df.columns:
        raise ValueError("DE results must contain 'cohens_d' column. "
                         "Run de() with effect_size=True.")

    gene_col = "gene" if "gene" in df.columns else df.index.name or "gene"
    if gene_col not in df.columns:
        df[gene_col] = df.index

    df["neg_log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8, 6))

    # Base points
    sig = df["padj"] < padj_threshold
    up = sig & (df["cohens_d"] > 0)
    down = sig & (df["cohens_d"] < 0)
    ns = ~sig

    ax.scatter(df.loc[ns, "cohens_d"], df.loc[ns, "neg_log10p"],
               c="lightgray", s=15, alpha=0.6, zorder=1)
    ax.scatter(df.loc[up, "cohens_d"], df.loc[up, "neg_log10p"],
               c="red", s=25, alpha=0.7, zorder=2, label="Up (FDR sig)")
    ax.scatter(df.loc[down, "cohens_d"], df.loc[down, "neg_log10p"],
               c="blue", s=25, alpha=0.7, zorder=2, label="Down (FDR sig)")

    # Highlight genes
    if highlight_genes:
        hl_mask = df[gene_col].isin(highlight_genes)
        ax.scatter(
            df.loc[hl_mask, "cohens_d"], df.loc[hl_mask, "neg_log10p"],
            c=highlight_color, s=40, zorder=4, edgecolor="black",
            linewidth=0.5, label=f"{highlight_label} ({hl_mask.sum()} genes)",
        )

    # Label top genes
    top = df.nlargest(label_top_n, "neg_log10p")
    for _, row in top.iterrows():
        ax.annotate(
            row[gene_col], (row["cohens_d"], row["neg_log10p"]),
            fontsize=8, alpha=0.8,
            xytext=(3, 3), textcoords="offset points",
        )

    # Reference lines
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.5,
               label="p = 0.05")
    ax.axvline(0, color="gray", linewidth=0.5)
    # Cohen's d thresholds
    for d_thresh in [-0.8, -0.5, 0.5, 0.8]:
        ax.axvline(d_thresh, color="lightgray", linewidth=0.5, linestyle=":")

    ax.set_xlabel("Cohen's d (effect size)", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(
        title or f"Effect Size Volcano ({sig.sum()} FDR sig)",
        fontsize=14, fontweight="bold",
    )
    ax.legend(loc="upper left", fontsize=10)

    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return fig
