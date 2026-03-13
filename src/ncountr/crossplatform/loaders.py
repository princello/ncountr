"""Load external expression data for cross-platform comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Union

import numpy as np
import pandas as pd


def load_expression(
    path: Union[str, Path],
    *,
    format: Literal["h5ad", "csv", "tsv"] | None = None,
    pseudobulk_group_by: str | None = None,
    gene_column: str | None = None,
) -> pd.DataFrame:
    """Load an external expression matrix.

    Parameters
    ----------
    path : str or Path
        Path to expression data.
    format : str, optional
        File format.  Inferred from extension if not provided.
    pseudobulk_group_by : str, optional
        For h5ad files, aggregate raw counts into pseudobulk by this obs column.
        The result is log2(CPM + 1).
    gene_column : str, optional
        For CSV/TSV, column to use as gene names (default: first column or index).

    Returns
    -------
    pd.DataFrame
        Genes (rows) x samples (columns), log2-transformed.
    """
    path = Path(path)
    if format is None:
        suffix = path.suffix.lower()
        format_map = {".h5ad": "h5ad", ".csv": "csv", ".tsv": "tsv", ".txt": "tsv"}
        format = format_map.get(suffix)
        if format is None:
            raise ValueError(f"Cannot infer format from extension {suffix!r}")

    if format == "h5ad":
        return _load_h5ad(path, pseudobulk_group_by=pseudobulk_group_by)
    elif format in ("csv", "tsv"):
        sep = "," if format == "csv" else "\t"
        df = pd.read_csv(path, sep=sep, index_col=0)
        return df
    else:
        raise ValueError(f"Unsupported format: {format!r}")


def _load_h5ad(
    path: Path, *, pseudobulk_group_by: str | None = None
) -> pd.DataFrame:
    """Load an h5ad file and optionally compute pseudobulk."""
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError(
            "scanpy is required for h5ad loading. "
            "Install with: pip install ncountr[crossplatform]"
        )

    adata = sc.read_h5ad(path)

    if pseudobulk_group_by is None:
        # Return the full expression matrix
        if adata.raw is not None:
            X = adata.raw.X
            genes = adata.raw.var_names
        else:
            X = adata.X
            genes = adata.var_names
        if hasattr(X, "toarray"):
            X = X.toarray()
        return pd.DataFrame(X.T, index=genes, columns=adata.obs_names)

    # Pseudobulk aggregation
    groups = adata.obs[pseudobulk_group_by].unique()
    pseudobulk: dict[str, np.ndarray] = {}

    for group in groups:
        cells = adata[adata.obs[pseudobulk_group_by] == group]
        if cells.n_obs == 0:
            continue
        raw_X = cells.raw.X if cells.raw is not None else cells.X
        if hasattr(raw_X, "toarray"):
            raw_X = raw_X.toarray()
        total = np.asarray(raw_X.sum(axis=0)).flatten()
        pseudobulk[str(group)] = total

    genes = adata.raw.var_names if adata.raw is not None else adata.var_names
    pb_df = pd.DataFrame(pseudobulk, index=genes)

    # Normalize to CPM then log2
    pb_cpm = pb_df.div(pb_df.sum(axis=0), axis=1) * 1e6
    pb_log2 = np.log2(pb_cpm + 1)

    # Map to gene symbols if available
    var_df = adata.raw.var if adata.raw is not None else adata.var
    if "symbol" in var_df.columns:
        symbol_map = var_df["symbol"].to_dict()
        pb_log2.index = [symbol_map.get(eid, eid) for eid in pb_log2.index]
        pb_log2 = pb_log2[pb_log2.index.notna()]
        pb_log2 = pb_log2.groupby(pb_log2.index).max()

    return pb_log2
