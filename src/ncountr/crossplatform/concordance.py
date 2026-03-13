"""DE direction concordance between platforms."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def de_concordance(
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
    *,
    gene_col_a: str = "gene",
    gene_col_b: str = "gene",
    lfc_col_a: str = "log2FC",
    lfc_col_b: str = "log2FC",
    padj_col_a: str = "padj",
    padj_col_b: str = "padj",
    gene_mapping: dict[str, str] | None = None,
    gene_flags: dict[str, bool] | None = None,
    flag_col_name: str = "is_flagged",
) -> pd.DataFrame:
    """Compare DE results between two platforms.

    Parameters
    ----------
    de_a, de_b : pd.DataFrame
        DE result tables from each platform.
    gene_col_a, gene_col_b : str
        Column holding gene names.
    lfc_col_a, lfc_col_b : str
        Column holding log2 fold-changes.
    padj_col_a, padj_col_b : str
        Column holding adjusted p-values.
    gene_mapping : dict, optional
        Mapping from de_a gene names to de_b gene names.
    gene_flags : dict[str, bool], optional
        Per-gene boolean flag (e.g., IFN gene membership).
    flag_col_name : str
        Name of the flag column in output.

    Returns
    -------
    pd.DataFrame
        Per-gene concordance with columns: gene, lfc_a, lfc_b, padj_a,
        padj_b, same_direction, <flag_col_name>.
    """
    a = de_a.set_index(gene_col_a) if gene_col_a in de_a.columns else de_a
    b = de_b.set_index(gene_col_b) if gene_col_b in de_b.columns else de_b

    if gene_mapping:
        mapped_genes = [g for g in a.index if g in gene_mapping and gene_mapping[g] in b.index]
    else:
        mapped_genes = sorted(set(a.index) & set(b.index))

    rows = []
    for g in mapped_genes:
        g_b = gene_mapping[g] if gene_mapping else g
        lfc_a = a.loc[g, lfc_col_a]
        lfc_b = b.loc[g_b, lfc_col_b]
        padj_a = a.loc[g, padj_col_a]
        padj_b = b.loc[g_b, padj_col_b]

        if lfc_a != 0 and lfc_b != 0:
            same_dir = (lfc_a > 0) == (lfc_b > 0)
        else:
            same_dir = np.nan

        row = {
            "gene": g,
            "lfc_a": lfc_a,
            "lfc_b": lfc_b,
            "padj_a": padj_a,
            "padj_b": padj_b,
            "same_direction": same_dir,
        }
        if gene_flags is not None:
            row[flag_col_name] = gene_flags.get(g, False)
        rows.append(row)

    return pd.DataFrame(rows)


def concordance_summary(conc_df: pd.DataFrame) -> dict:
    """Compute concordance statistics.

    Parameters
    ----------
    conc_df : pd.DataFrame
        Output of :func:`de_concordance`.

    Returns
    -------
    dict
        Keys: overall_rate, n_concordant, n_total, lfc_spearman_r, lfc_spearman_p.
    """
    valid = conc_df.dropna(subset=["same_direction"]).copy()
    valid["same_direction"] = valid["same_direction"].astype(float)

    n_concordant = int(valid["same_direction"].sum())
    n_total = len(valid)
    rate = n_concordant / n_total if n_total > 0 else np.nan

    r, p = stats.spearmanr(valid["lfc_a"], valid["lfc_b"]) if n_total >= 3 else (np.nan, np.nan)

    return {
        "overall_rate": rate,
        "n_concordant": n_concordant,
        "n_total": n_total,
        "lfc_spearman_r": r,
        "lfc_spearman_p": p,
    }
