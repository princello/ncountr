"""Gene set / pathway scoring."""

from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd

from ncountr.experiment import NanostringExperiment
from ncountr.datasets import get_gene_set


def score_gene_set(
    experiment: NanostringExperiment,
    *,
    gene_set: Union[str, list[str]],
    counts: pd.DataFrame | None = None,
    samples: list[str] | None = None,
    method: str = "zscore_mean",
) -> pd.Series:
    """Score samples for a gene set.

    Parameters
    ----------
    experiment : NanostringExperiment
    gene_set : str or list[str]
        A built-in gene set name (e.g. ``"IFN_JAKSTAT"``) or an explicit
        list of gene names.
    counts : pd.DataFrame, optional
        Count matrix (genes x samples).  Defaults to normalized or raw counts.
    samples : list[str], optional
        Subset of samples to score.  Defaults to all.
    method : str
        Scoring method.  Currently ``"zscore_mean"`` (z-score each gene
        across samples, then take the mean z-score per sample).

    Returns
    -------
    pd.Series
        Score per sample.
    """
    if isinstance(gene_set, str):
        genes = get_gene_set(gene_set)
    else:
        genes = list(gene_set)

    if counts is None:
        counts = experiment.normalized if experiment.normalized is not None else experiment.raw_counts

    if samples is None:
        samples = experiment.samples

    # Filter to genes present in the data
    available = [g for g in genes if g in counts.index]
    if not available:
        raise ValueError(f"No genes from the gene set are present in the data")

    mat = np.log2(counts.loc[available, samples].astype(float) + 1)

    if method == "zscore_mean":
        row_mean = mat.mean(axis=1)
        row_std = mat.std(axis=1).replace(0, np.nan)
        z = mat.subtract(row_mean, axis=0).div(row_std, axis=0)
        scores = z.mean(axis=0)
    else:
        raise ValueError(f"Unknown scoring method: {method!r}")

    scores.name = "pathway_score"
    return scores
