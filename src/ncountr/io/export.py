"""Export utilities for writing results to disk."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Union

import pandas as pd

from ncountr.experiment import NanostringExperiment


def to_anndata(experiment: NanostringExperiment):
    """Convert a NanostringExperiment to an AnnData object.

    Uses normalized counts as the main matrix (X) when available,
    otherwise falls back to raw counts.  Raw counts are always stored
    in ``adata.layers["raw"]``.

    Parameters
    ----------
    experiment : NanostringExperiment
        A parsed (and optionally normalized) experiment.

    Returns
    -------
    anndata.AnnData
        Samples in ``.obs``, genes in ``.var``.

    Raises
    ------
    ImportError
        If ``anndata`` is not installed.
    """
    try:
        import anndata
    except ImportError:
        raise ImportError(
            "anndata is required for to_anndata(). "
            "Install it with: pip install ncountr[crossplatform]"
        )

    import numpy as np

    # Combine endogenous + housekeeping genes for completeness
    raw_full = pd.concat([experiment.raw_counts, experiment.hk_counts])

    # Determine the main count matrix (samples x genes)
    if experiment.normalized is not None:
        # Normalized may only include endogenous genes — align to raw_full index
        X_df = experiment.normalized.reindex(raw_full.index, fill_value=0)
    else:
        X_df = raw_full

    # AnnData expects obs (samples) x var (genes)
    X = X_df.T.values.astype(np.float32)

    # Build obs (sample metadata)
    obs = experiment.sample_meta.copy() if not experiment.sample_meta.empty else pd.DataFrame(index=X_df.columns)
    obs.index = obs.index.astype(str)
    # Merge QC results if available
    if experiment.qc_results is not None:
        for col in experiment.qc_results.columns:
            obs[f"qc_{col}"] = experiment.qc_results[col].reindex(obs.index)
    # Merge lane info
    if not experiment.lane_info.empty:
        for col in experiment.lane_info.columns:
            obs[f"lane_{col}"] = experiment.lane_info[col].reindex(obs.index)

    # Build var (gene metadata)
    gene_names = list(raw_full.index)
    hk_set = set(experiment.hk_counts.index)
    var = pd.DataFrame(
        {"housekeeping": [g in hk_set for g in gene_names]},
        index=gene_names,
    )
    var.index.name = "gene"

    adata = anndata.AnnData(
        X=X,
        obs=obs,
        var=var,
    )
    adata.layers["raw"] = raw_full.T.reindex(index=obs.index, columns=gene_names).values.astype(np.float32)

    # Store positive/negative controls in uns
    if not experiment.pos_counts.empty:
        adata.uns["pos_counts"] = experiment.pos_counts.to_dict()
    if not experiment.neg_counts.empty:
        adata.uns["neg_counts"] = experiment.neg_counts.to_dict()

    return adata


def export_counts(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> dict[str, Path]:
    """Write raw and normalized count matrices to CSV.

    Returns a dict mapping description to output path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    raw_path = output_dir / f"{prefix}_raw_counts.csv"
    combined = pd.concat([experiment.raw_counts, experiment.hk_counts])
    combined.to_csv(raw_path)
    paths["raw_counts"] = raw_path

    if experiment.normalized is not None:
        norm_path = output_dir / f"{prefix}_normalized_counts.csv"
        experiment.normalized.to_csv(norm_path)
        paths["normalized_counts"] = norm_path

    return paths


def export_qc(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> Path | None:
    """Write QC results to CSV."""
    if experiment.qc_results is None:
        return None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{prefix}_qc_results.csv"
    experiment.qc_results.to_csv(path)
    return path


def export_de(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> Path | None:
    """Write DE results to CSV."""
    if experiment.de_results is None:
        return None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{prefix}_de_results.csv"
    experiment.de_results.to_csv(path, index=False)
    return path
