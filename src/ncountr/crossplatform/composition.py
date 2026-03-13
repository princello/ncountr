"""Cell composition proxy from marker gene expression."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def marker_composition_proxy(
    nanostring_counts: pd.DataFrame,
    cell_proportions: pd.DataFrame,
    markers: dict[str, list[str]],
    *,
    samples: list[str] | None = None,
) -> pd.DataFrame:
    """Correlate Nanostring marker gene expression with cell proportions.

    Parameters
    ----------
    nanostring_counts : pd.DataFrame
        Nanostring expression (genes x samples).
    cell_proportions : pd.DataFrame
        Cell type proportions (samples x cell types).
    markers : dict[str, list[str]]
        Cell type name → list of marker gene names.
    samples : list[str], optional
        Subset of samples.  If None, uses intersection.

    Returns
    -------
    pd.DataFrame
        Per-cell-type correlation: cell_type, markers_used, matched_column,
        spearman_r, spearman_p.
    """
    if samples is None:
        samples = sorted(set(nanostring_counts.columns) & set(cell_proportions.index))

    rows = []
    for ct_name, marker_genes in markers.items():
        available = [m for m in marker_genes if m in nanostring_counts.index]
        if not available:
            continue

        nano_signal = nanostring_counts.loc[available, samples].mean(axis=0).values

        # Find matching column in cell_proportions
        ct_col = None
        for col in cell_proportions.columns:
            if ct_name.lower().replace("_", " ") in col.lower():
                ct_col = col
                break

        if ct_col is None:
            continue

        sc_vals = cell_proportions.loc[samples, ct_col].values.astype(float)

        if np.std(nano_signal) > 0 and np.std(sc_vals) > 0:
            r, p = stats.spearmanr(nano_signal, sc_vals)
        else:
            r, p = np.nan, np.nan

        rows.append({
            "cell_type": ct_name,
            "markers_used": ", ".join(available),
            "matched_column": ct_col,
            "spearman_r": r,
            "spearman_p": p,
        })

    return pd.DataFrame(rows)
