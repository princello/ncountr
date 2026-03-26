"""Differential expression analysis for Nanostring nCounter data."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ncountr.experiment import NanostringExperiment


def effect_sizes(
    counts: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    *,
    n_bootstrap: int = 2000,
    ci_level: float = 0.95,
    seed: int = 42,
) -> pd.DataFrame:
    """Compute effect sizes for each gene.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene-by-sample count matrix.
    group_a, group_b : list[str]
        Sample IDs.
    n_bootstrap : int
        Bootstrap iterations for Cohen's d confidence intervals.
    ci_level : float
        Confidence interval level (0.95 = 95%).
    seed : int
        Random seed.

    Returns
    -------
    pd.DataFrame
        Columns: gene, cohens_d, cohens_d_ci_lo, cohens_d_ci_hi,
        rank_biserial.
    """
    rng = np.random.default_rng(seed)
    alpha = (1 - ci_level) / 2
    results = []

    for gene in counts.index:
        va = counts.loc[gene, group_a].values.astype(float)
        vb = counts.loc[gene, group_b].values.astype(float)
        na, nb = len(va), len(vb)

        # Cohen's d (pooled SD)
        pooled_std = np.sqrt(
            ((na - 1) * np.var(va, ddof=1) + (nb - 1) * np.var(vb, ddof=1))
            / max(na + nb - 2, 1)
        )
        d = (np.mean(va) - np.mean(vb)) / pooled_std if pooled_std > 0 else 0.0

        # Bootstrap CI for Cohen's d
        boot_d = np.empty(n_bootstrap)
        for i in range(n_bootstrap):
            ba = rng.choice(va, size=na, replace=True)
            bb = rng.choice(vb, size=nb, replace=True)
            ps = np.sqrt(
                ((na - 1) * np.var(ba, ddof=1) + (nb - 1) * np.var(bb, ddof=1))
                / max(na + nb - 2, 1)
            )
            boot_d[i] = (np.mean(ba) - np.mean(bb)) / ps if ps > 0 else 0.0

        ci_lo = np.percentile(boot_d, 100 * alpha)
        ci_hi = np.percentile(boot_d, 100 * (1 - alpha))

        # Rank-biserial correlation (effect size for Mann-Whitney U)
        try:
            u_stat, _ = stats.mannwhitneyu(va, vb, alternative="two-sided")
            rank_biserial = 1 - 2 * u_stat / (na * nb)
        except ValueError:
            rank_biserial = 0.0

        results.append({
            "gene": gene,
            "cohens_d": d,
            "cohens_d_ci_lo": ci_lo,
            "cohens_d_ci_hi": ci_hi,
            "rank_biserial": rank_biserial,
        })

    return pd.DataFrame(results)


def de(
    experiment: NanostringExperiment,
    *,
    group_a: list[str],
    group_b: list[str],
    counts: pd.DataFrame | None = None,
    test: Literal["mannwhitneyu", "ttest"] = "mannwhitneyu",
    correction: str = "fdr_bh",
    effect_size: bool = False,
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
    effect_size : bool
        If True, compute and include Cohen's d (with bootstrap 95% CI)
        and rank-biserial correlation.  Adds columns: cohens_d,
        cohens_d_ci_lo, cohens_d_ci_hi, rank_biserial.
    store : bool
        If True, store results on ``experiment.de_results``.

    Returns
    -------
    pd.DataFrame
        Columns: gene, log2FC, mean_a, mean_b, pvalue, padj
        (plus effect size columns if ``effect_size=True``).
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

    if effect_size:
        es_df = effect_sizes(counts, group_a, group_b)
        de_df = de_df.merge(es_df, on="gene", how="left")

    if store:
        experiment.de_results = de_df
    return de_df
