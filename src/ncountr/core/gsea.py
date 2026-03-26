"""Gene Set Enrichment Analysis for Nanostring nCounter data.

Implements permutation-based GSEA (Subramanian et al. 2005), competitive
gene set testing (CAMERA-like), and self-contained pathway tests.
Optimized for the small sample sizes typical of nCounter experiments
(n=3-12 per group).
"""

from __future__ import annotations

import warnings
from itertools import combinations
from typing import Literal, Union

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ncountr.experiment import NanostringExperiment


# ---------------------------------------------------------------------------
# Gene ranking
# ---------------------------------------------------------------------------

def rank_genes(
    counts: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    *,
    metric: Literal["signal_to_noise", "log2fc_stat", "log2fc"] = "signal_to_noise",
) -> pd.Series:
    """Rank genes by a differential expression metric.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene-by-sample count matrix.
    group_a, group_b : list[str]
        Sample IDs for the two groups.
    metric : str
        ``"signal_to_noise"`` — (mean_a - mean_b) / (std_a + std_b).
        ``"log2fc_stat"`` — log2FC * -log10(p) from Mann-Whitney U.
        ``"log2fc"`` — simple log2 fold change.

    Returns
    -------
    pd.Series
        Metric values indexed by gene name, sorted descending.
    """
    vals_a = counts[group_a].astype(float)
    vals_b = counts[group_b].astype(float)

    mean_a = vals_a.mean(axis=1)
    mean_b = vals_b.mean(axis=1)

    if metric == "signal_to_noise":
        std_a = vals_a.std(axis=1, ddof=1).clip(lower=0.5)
        std_b = vals_b.std(axis=1, ddof=1).clip(lower=0.5)
        scores = (mean_a - mean_b) / (std_a + std_b)

    elif metric == "log2fc":
        scores = np.log2((mean_a + 1) / (mean_b + 1))

    elif metric == "log2fc_stat":
        log2fc = np.log2((mean_a + 1) / (mean_b + 1))
        pvals = pd.Series(1.0, index=counts.index)
        for gene in counts.index:
            va = vals_a.loc[gene].values
            vb = vals_b.loc[gene].values
            if np.std(va) > 0 or np.std(vb) > 0:
                try:
                    _, p = stats.mannwhitneyu(va, vb, alternative="two-sided")
                    pvals[gene] = max(p, 1e-300)
                except ValueError:
                    pass
        scores = log2fc * (-np.log10(pvals))

    else:
        raise ValueError(f"Unknown metric: {metric!r}")

    scores.name = "rank_metric"
    return scores.sort_values(ascending=False)


# ---------------------------------------------------------------------------
# Core GSEA algorithm
# ---------------------------------------------------------------------------

def _running_enrichment_score(
    ranked_metric: pd.Series,
    gene_set: list[str],
    *,
    weighted_score_type: float = 1.0,
) -> tuple[float, np.ndarray, list[str]]:
    """Compute the running enrichment score for a gene set.

    Parameters
    ----------
    ranked_metric : pd.Series
        Gene names as index, metric values as data (sorted descending).
    gene_set : list[str]
        Genes in the set.
    weighted_score_type : float
        Weight exponent *p* (1.0 = standard GSEA).

    Returns
    -------
    es : float
        Enrichment score (max deviation from zero).
    running : np.ndarray
        Running enrichment score at each position.
    leading_edge : list[str]
        Leading edge genes (those contributing to the peak).
    """
    genes = ranked_metric.index.tolist()
    metric_vals = np.abs(ranked_metric.values)
    n = len(genes)
    gene_set_mask = np.array([g in gene_set for g in genes], dtype=bool)
    n_hit = gene_set_mask.sum()

    if n_hit == 0:
        return 0.0, np.zeros(n), []

    n_miss = n - n_hit

    # Weighted hit score
    hit_vals = np.where(gene_set_mask, metric_vals ** weighted_score_type, 0.0)
    hit_sum = hit_vals.sum()
    if hit_sum == 0:
        hit_sum = 1.0  # avoid division by zero

    # Running sum
    running = np.zeros(n)
    for i in range(n):
        if gene_set_mask[i]:
            running[i] = hit_vals[i] / hit_sum
        else:
            running[i] = -1.0 / n_miss

    running = np.cumsum(running)

    # Enrichment score: max deviation from zero
    max_pos = running.max()
    min_neg = running.min()
    if abs(max_pos) >= abs(min_neg):
        es = max_pos
        peak_idx = np.argmax(running)
    else:
        es = min_neg
        peak_idx = np.argmin(running)

    # Leading edge: genes up to and including the peak
    if es >= 0:
        leading_edge = [genes[i] for i in range(peak_idx + 1) if gene_set_mask[i]]
    else:
        leading_edge = [genes[i] for i in range(peak_idx, n) if gene_set_mask[i]]

    return es, running, leading_edge


def _permutation_null(
    counts: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    gene_set: list[str],
    *,
    metric: str = "signal_to_noise",
    n_perm: int = 1000,
    weighted_score_type: float = 1.0,
    seed: int = 42,
) -> np.ndarray:
    """Build null distribution by permuting sample labels.

    For small total n (≤10), enumerates all unique permutations exactly.
    """
    all_samples = list(group_a) + list(group_b)
    n_a = len(group_a)
    n_total = len(all_samples)
    rng = np.random.default_rng(seed)

    # Exact enumeration for small n
    from math import comb
    n_unique = comb(n_total, n_a)
    exact = n_unique <= max(n_perm, 200)

    if exact:
        perms = list(combinations(range(n_total), n_a))
    else:
        perms = []
        for _ in range(n_perm):
            idx = rng.permutation(n_total)
            perms.append(tuple(idx[:n_a]))

    null_es = np.empty(len(perms))
    for i, perm_a_idx in enumerate(perms):
        perm_b_idx = [j for j in range(n_total) if j not in perm_a_idx]
        perm_a = [all_samples[j] for j in perm_a_idx]
        perm_b = [all_samples[j] for j in perm_b_idx]
        ranked = rank_genes(counts, perm_a, perm_b, metric=metric)
        es, _, _ = _running_enrichment_score(
            ranked, gene_set, weighted_score_type=weighted_score_type
        )
        null_es[i] = es

    return null_es


def _compute_pvalue_and_nes(
    es: float, null_es: np.ndarray
) -> tuple[float, float]:
    """Compute nominal p-value and NES from enrichment score and null."""
    if es >= 0:
        pos_null = null_es[null_es >= 0]
        if len(pos_null) < 1:
            return 1.0, es
        pval = (pos_null >= es).sum() / len(pos_null)
        mean_pos = pos_null.mean() if pos_null.mean() != 0 else 1.0
        nes = es / abs(mean_pos)
    else:
        neg_null = null_es[null_es < 0]
        if len(neg_null) < 1:
            return 1.0, es
        pval = (neg_null <= es).sum() / len(neg_null)
        mean_neg = neg_null.mean() if neg_null.mean() != 0 else -1.0
        nes = es / abs(mean_neg)

    return max(pval, 1.0 / (len(null_es) + 1)), nes


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def gsea(
    experiment: NanostringExperiment,
    *,
    gene_sets: dict[str, list[str]] | None = None,
    group_a: list[str],
    group_b: list[str],
    counts: pd.DataFrame | None = None,
    metric: Literal["signal_to_noise", "log2fc_stat", "log2fc"] = "signal_to_noise",
    n_perm: int = 1000,
    weighted_score_type: float = 1.0,
    min_set_size: int = 5,
    max_set_size: int = 500,
    correction: str = "fdr_bh",
    seed: int = 42,
    store: bool = True,
) -> pd.DataFrame:
    """Run permutation-based GSEA (Subramanian et al. 2005).

    Parameters
    ----------
    experiment : NanostringExperiment
    gene_sets : dict[str, list[str]], optional
        Gene set name to gene list mapping.  If *None*, uses all built-in
        gene sets (custom + Hallmark).
    group_a, group_b : list[str]
        Sample IDs.  Positive ES means enriched in *group_a*.
    counts : pd.DataFrame, optional
        Count matrix; defaults to normalized or raw counts.
    metric : str
        Gene ranking metric.
    n_perm : int
        Number of permutations (ignored when exact enumeration is used).
    weighted_score_type : float
        GSEA weight parameter *p* (1.0 = standard).
    min_set_size, max_set_size : int
        Filter gene sets by overlap with measured genes.
    correction : str
        FDR correction method.
    seed : int
        Random seed.
    store : bool
        If True, stores result on ``experiment.gsea_results``.

    Returns
    -------
    pd.DataFrame
        Columns: gene_set, es, nes, pvalue, padj, n_genes, n_overlap,
        leading_edge.
    """
    if counts is None:
        counts = (
            experiment.normalized
            if experiment.normalized is not None
            else experiment.raw_counts
        )

    if gene_sets is None:
        from ncountr.genesets import get_all_gene_sets
        gene_sets = get_all_gene_sets()

    from ncountr.genesets import filter_gene_sets
    measured = list(counts.index)
    filtered = filter_gene_sets(
        gene_sets, measured, min_overlap=min_set_size, max_overlap=max_set_size
    )

    if not filtered:
        warnings.warn("No gene sets passed the overlap filter.")
        return pd.DataFrame(
            columns=["gene_set", "es", "nes", "pvalue", "padj",
                     "n_genes", "n_overlap", "leading_edge"]
        )

    # Warn about exact permutation limits
    from math import comb
    n_total = len(group_a) + len(group_b)
    n_unique = comb(n_total, len(group_a))
    if n_unique <= 20:
        warnings.warn(
            f"Only {n_unique} unique sample permutations exist "
            f"(n={len(group_a)}+{len(group_b)}). Using exact enumeration. "
            f"Minimum achievable p-value is {1/n_unique:.3f}."
        )

    # Rank genes once for observed ES
    ranked = rank_genes(counts, group_a, group_b, metric=metric)

    results = []
    for gs_name, gs_genes in filtered.items():
        # Observed ES
        es, running, le = _running_enrichment_score(
            ranked, gs_genes, weighted_score_type=weighted_score_type
        )

        # Null distribution
        null_es = _permutation_null(
            counts, group_a, group_b, gs_genes,
            metric=metric, n_perm=n_perm,
            weighted_score_type=weighted_score_type, seed=seed,
        )

        pval, nes = _compute_pvalue_and_nes(es, null_es)

        results.append({
            "gene_set": gs_name,
            "es": es,
            "nes": nes,
            "pvalue": pval,
            "n_genes": len(gene_sets[gs_name]),
            "n_overlap": len(gs_genes),
            "leading_edge": "|".join(le),
            "_running": running,
            "_ranked": ranked,
        })

    df = pd.DataFrame(results)

    # FDR correction
    if len(df) > 0:
        _, df["padj"], _, _ = multipletests(df["pvalue"], method=correction)
    else:
        df["padj"] = []

    # Clean up internal columns before returning
    internal = df[["_running", "_ranked"]].copy()
    df = df.drop(columns=["_running", "_ranked"])
    df = df.sort_values("pvalue").reset_index(drop=True)

    # Stash internals for plotting
    df.attrs["_internals"] = {
        row["gene_set"]: {"running": internal.iloc[i]["_running"],
                          "ranked": internal.iloc[i]["_ranked"]}
        for i, (_, row) in enumerate(
            pd.DataFrame(results).iterrows()
        )
    }

    if store:
        experiment.gsea_results = df

    return df


def gsea_preranked(
    ranked_genes: pd.Series,
    gene_sets: dict[str, list[str]],
    *,
    n_perm: int = 1000,
    weighted_score_type: float = 1.0,
    min_set_size: int = 5,
    max_set_size: int = 500,
    correction: str = "fdr_bh",
    seed: int = 42,
) -> pd.DataFrame:
    """Run GSEA on a pre-ranked gene list.

    Parameters
    ----------
    ranked_genes : pd.Series
        Gene names as index, ranking metric as values (sorted descending).
    gene_sets : dict[str, list[str]]
        Gene set name to gene list mapping.

    Returns
    -------
    pd.DataFrame
        Same columns as :func:`gsea`.
    """
    from ncountr.genesets import filter_gene_sets

    measured = list(ranked_genes.index)
    filtered = filter_gene_sets(
        gene_sets, measured, min_overlap=min_set_size, max_overlap=max_set_size
    )

    if not ranked_genes.index.is_monotonic_decreasing:
        ranked_genes = ranked_genes.sort_values(ascending=False)

    rng = np.random.default_rng(seed)
    results = []
    for gs_name, gs_genes in filtered.items():
        es, running, le = _running_enrichment_score(
            ranked_genes, gs_genes, weighted_score_type=weighted_score_type
        )

        # Permutation: shuffle gene labels
        null_es = np.empty(n_perm)
        for i in range(n_perm):
            perm_idx = rng.permutation(len(ranked_genes))
            shuffled = pd.Series(
                ranked_genes.values, index=ranked_genes.index[perm_idx]
            )
            null_es[i], _, _ = _running_enrichment_score(
                shuffled, gs_genes, weighted_score_type=weighted_score_type
            )

        pval, nes = _compute_pvalue_and_nes(es, null_es)

        results.append({
            "gene_set": gs_name,
            "es": es,
            "nes": nes,
            "pvalue": pval,
            "n_genes": len(gene_sets[gs_name]),
            "n_overlap": len(gs_genes),
            "leading_edge": "|".join(le),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, df["padj"], _, _ = multipletests(df["pvalue"], method=correction)
    else:
        df["padj"] = []
    return df.sort_values("pvalue").reset_index(drop=True)


def competitive_test(
    experiment: NanostringExperiment,
    *,
    gene_sets: dict[str, list[str]] | None = None,
    group_a: list[str],
    group_b: list[str],
    counts: pd.DataFrame | None = None,
    min_set_size: int = 5,
    correction: str = "fdr_bh",
) -> pd.DataFrame:
    """Competitive gene set test (CAMERA-like).

    Tests whether genes in a set are more differentially expressed than
    genes outside the set, using Wilcoxon rank-sum on per-gene test
    statistics.

    Parameters
    ----------
    experiment : NanostringExperiment
    gene_sets : dict, optional
        If None, uses all built-in gene sets.
    group_a, group_b : list[str]
        Sample IDs.
    counts : pd.DataFrame, optional
    min_set_size : int
    correction : str

    Returns
    -------
    pd.DataFrame
        Columns: gene_set, direction, stat, pvalue, padj, n_overlap.
    """
    if counts is None:
        counts = (
            experiment.normalized
            if experiment.normalized is not None
            else experiment.raw_counts
        )

    if gene_sets is None:
        from ncountr.genesets import get_all_gene_sets
        gene_sets = get_all_gene_sets()

    # Compute per-gene t-statistics
    vals_a = counts[group_a].astype(float)
    vals_b = counts[group_b].astype(float)
    mean_a = vals_a.mean(axis=1)
    mean_b = vals_b.mean(axis=1)
    std_a = vals_a.std(axis=1, ddof=1).clip(lower=0.5)
    std_b = vals_b.std(axis=1, ddof=1).clip(lower=0.5)
    n_a, n_b = len(group_a), len(group_b)
    pooled_se = np.sqrt(std_a**2 / n_a + std_b**2 / n_b)
    t_stats = (mean_a - mean_b) / pooled_se

    measured = set(counts.index)
    results = []
    for gs_name, gs_genes in gene_sets.items():
        overlap = [g for g in gs_genes if g in measured]
        if len(overlap) < min_set_size:
            continue

        in_set = t_stats.loc[overlap].values
        out_set = t_stats.loc[~t_stats.index.isin(overlap)].values

        try:
            stat, pval = stats.mannwhitneyu(
                np.abs(in_set), np.abs(out_set), alternative="greater"
            )
        except ValueError:
            stat, pval = 0.0, 1.0

        direction = "up" if np.mean(in_set) > 0 else "down"

        results.append({
            "gene_set": gs_name,
            "direction": direction,
            "stat": stat,
            "pvalue": pval,
            "n_overlap": len(overlap),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, df["padj"], _, _ = multipletests(df["pvalue"], method=correction)
    else:
        df["padj"] = []
    return df.sort_values("pvalue").reset_index(drop=True)


def self_contained_test(
    experiment: NanostringExperiment,
    *,
    gene_sets: dict[str, list[str]] | None = None,
    group_a: list[str],
    group_b: list[str],
    counts: pd.DataFrame | None = None,
    n_perm: int = 10000,
    min_set_size: int = 5,
    correction: str = "fdr_bh",
    seed: int = 42,
) -> pd.DataFrame:
    """Self-contained gene set test via permutation.

    Tests whether the mean z-score pathway score differs between groups
    by permuting sample labels.

    Parameters
    ----------
    experiment : NanostringExperiment
    gene_sets : dict, optional
    group_a, group_b : list[str]
    counts : pd.DataFrame, optional
    n_perm : int
    min_set_size : int
    correction : str
    seed : int

    Returns
    -------
    pd.DataFrame
        Columns: gene_set, mean_diff, pvalue, padj, n_overlap.
    """
    if counts is None:
        counts = (
            experiment.normalized
            if experiment.normalized is not None
            else experiment.raw_counts
        )

    if gene_sets is None:
        from ncountr.genesets import get_all_gene_sets
        gene_sets = get_all_gene_sets()

    all_samples = list(group_a) + list(group_b)
    n_a = len(group_a)
    n_total = len(all_samples)
    rng = np.random.default_rng(seed)
    measured = set(counts.index)

    # Exact enumeration for small n
    from math import comb
    n_unique = comb(n_total, n_a)
    exact = n_unique <= max(n_perm, 200)

    if exact:
        perms = list(combinations(range(n_total), n_a))
    else:
        perms = [tuple(rng.permutation(n_total)[:n_a]) for _ in range(n_perm)]

    results = []
    for gs_name, gs_genes in gene_sets.items():
        overlap = [g for g in gs_genes if g in measured]
        if len(overlap) < min_set_size:
            continue

        # Z-score the gene set genes across all samples
        mat = np.log2(counts.loc[overlap, all_samples].astype(float) + 1)
        row_mean = mat.mean(axis=1)
        row_std = mat.std(axis=1).replace(0, np.nan)
        z = mat.subtract(row_mean, axis=0).div(row_std, axis=0).fillna(0)
        sample_scores = z.mean(axis=0)

        # Observed difference
        obs_a = sample_scores[group_a].mean()
        obs_b = sample_scores[group_b].mean()
        obs_diff = obs_a - obs_b

        # Permutation null
        null_diffs = np.empty(len(perms))
        score_vals = sample_scores.values
        for i, perm_a_idx in enumerate(perms):
            perm_b_idx = [j for j in range(n_total) if j not in perm_a_idx]
            null_diffs[i] = (
                score_vals[list(perm_a_idx)].mean()
                - score_vals[list(perm_b_idx)].mean()
            )

        # Two-sided p-value
        pval = (np.abs(null_diffs) >= abs(obs_diff)).sum() / len(null_diffs)
        pval = max(pval, 1.0 / (len(perms) + 1))

        results.append({
            "gene_set": gs_name,
            "mean_diff": obs_diff,
            "pvalue": pval,
            "n_overlap": len(overlap),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, df["padj"], _, _ = multipletests(df["pvalue"], method=correction)
    else:
        df["padj"] = []
    return df.sort_values("pvalue").reset_index(drop=True)
