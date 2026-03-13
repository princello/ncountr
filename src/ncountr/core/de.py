"""Differential expression analysis for Nanostring nCounter data."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ncountr.experiment import NanostringExperiment


def de(
    experiment: NanostringExperiment,
    *,
    group_a: list[str],
    group_b: list[str],
    counts: pd.DataFrame | None = None,
    test: Literal["mannwhitneyu", "ttest"] = "mannwhitneyu",
    correction: str = "fdr_bh",
    store: bool = True,
) -> pd.DataFrame:
    """Run differential expression between two sample groups.

    Parameters
    ----------
    experiment : NanostringExperiment
    group_a, group_b : list[str]
        Sample IDs for the two groups.  Log2FC is computed as
        ``log2(mean_a + 1) - log2(mean_b + 1)``, i.e. positive values mean
        higher in group A.
    counts : pd.DataFrame, optional
        Count matrix to use.  Defaults to ``experiment.normalized`` if
        available, otherwise ``experiment.raw_counts``.
    test : str
        Statistical test: ``"mannwhitneyu"`` or ``"ttest"``.
    correction : str
        Multiple testing correction method (passed to
        ``statsmodels.stats.multitest.multipletests``).
    store : bool
        If True, store results on ``experiment.de_results``.

    Returns
    -------
    pd.DataFrame
        Columns: gene, log2FC, mean_a, mean_b, pvalue, padj.
    """
    if counts is None:
        counts = experiment.normalized if experiment.normalized is not None else experiment.raw_counts

    # Validate sample IDs
    missing_a = set(group_a) - set(counts.columns)
    missing_b = set(group_b) - set(counts.columns)
    if missing_a:
        raise ValueError(f"group_a samples not in counts: {missing_a}")
    if missing_b:
        raise ValueError(f"group_b samples not in counts: {missing_b}")

    results: list[dict] = []
    for gene in counts.index:
        vals_a = counts.loc[gene, group_a].values.astype(float)
        vals_b = counts.loc[gene, group_b].values.astype(float)

        mean_a = vals_a.mean()
        mean_b = vals_b.mean()
        log2fc = np.log2((mean_a + 1) / (mean_b + 1))

        if np.std(vals_a) > 0 or np.std(vals_b) > 0:
            try:
                if test == "mannwhitneyu":
                    _stat, p_val = stats.mannwhitneyu(
                        vals_a, vals_b, alternative="two-sided"
                    )
                else:
                    _stat, p_val = stats.ttest_ind(vals_a, vals_b)
            except ValueError:
                p_val = 1.0
        else:
            p_val = 1.0

        results.append(
            {
                "gene": gene,
                "log2FC": log2fc,
                "mean_a": mean_a,
                "mean_b": mean_b,
                "pvalue": p_val,
            }
        )

    de_df = pd.DataFrame(results)
    _, de_df["padj"], _, _ = multipletests(de_df["pvalue"], method=correction)
    de_df = de_df.sort_values("pvalue").reset_index(drop=True)

    if store:
        experiment.de_results = de_df
    return de_df
