"""Cross-platform correlation analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def per_sample_correlation(
    nanostring: pd.DataFrame,
    external: pd.DataFrame,
    *,
    shared_genes: list[str] | None = None,
    method: str = "spearman",
) -> pd.DataFrame:
    """Compute per-sample correlation between two expression matrices.

    Parameters
    ----------
    nanostring : pd.DataFrame
        Nanostring expression (genes x samples).
    external : pd.DataFrame
        External expression (genes x samples).
    shared_genes : list[str], optional
        Genes to include.  If None, uses intersection of indices.
    method : str
        ``"spearman"`` or ``"pearson"``.

    Returns
    -------
    pd.DataFrame
        One row per sample with columns: sample, r, pvalue.
    """
    if shared_genes is None:
        shared_genes = sorted(set(nanostring.index) & set(external.index))

    shared_samples = sorted(set(nanostring.columns) & set(external.columns))

    corr_func = stats.spearmanr if method == "spearman" else stats.pearsonr

    rows = []
    for sid in shared_samples:
        n_vals = nanostring.loc[shared_genes, sid].values.astype(float)
        e_vals = external.loc[shared_genes, sid].values.astype(float)
        mask = np.isfinite(n_vals) & np.isfinite(e_vals)
        if mask.sum() < 3:
            continue
        r, p = corr_func(n_vals[mask], e_vals[mask])
        rows.append({"sample": sid, "r": r, "pvalue": p})

    return pd.DataFrame(rows)


def per_gene_correlation(
    nanostring: pd.DataFrame,
    external: pd.DataFrame,
    *,
    shared_samples: list[str] | None = None,
    method: str = "spearman",
    min_samples: int = 4,
) -> pd.DataFrame:
    """Compute per-gene correlation across shared samples.

    Parameters
    ----------
    nanostring, external : pd.DataFrame
        Expression matrices (genes x samples).
    shared_samples : list[str], optional
        Samples to include.  If None, uses intersection.
    method : str
        ``"spearman"`` or ``"pearson"``.
    min_samples : int
        Minimum samples with variation required.

    Returns
    -------
    pd.DataFrame
        One row per gene with columns: gene, r, pvalue.
    """
    if shared_samples is None:
        shared_samples = sorted(set(nanostring.columns) & set(external.columns))

    shared_genes = sorted(set(nanostring.index) & set(external.index))
    corr_func = stats.spearmanr if method == "spearman" else stats.pearsonr

    rows = []
    for gene in shared_genes:
        n_vals = nanostring.loc[gene, shared_samples].values.astype(float)
        e_vals = external.loc[gene, shared_samples].values.astype(float)
        mask = np.isfinite(n_vals) & np.isfinite(e_vals)
        if mask.sum() < min_samples:
            continue
        if np.std(n_vals[mask]) == 0 or np.std(e_vals[mask]) == 0:
            continue
        r, p = corr_func(n_vals[mask], e_vals[mask])
        rows.append({"gene": gene, "r": r, "pvalue": p})

    return pd.DataFrame(rows)
