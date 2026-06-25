"""Tests for ncountr.core.gsea — GSEA and gene set testing."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ncountr.core.gsea import (
    rank_genes,
    gsea,
    gsea_preranked,
    competitive_test,
    self_contained_test,
    _running_enrichment_score,
    _compute_pvalue_and_nes,
)
from ncountr.datasets import get_gene_set


IFN_SET = ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1", "CXCL10"]


# ---------------------------------------------------------------------------
# rank_genes
# ---------------------------------------------------------------------------


class TestRankGenes:
    def test_returns_series_sorted_descending(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        assert isinstance(ranked, pd.Series)
        assert ranked.is_monotonic_decreasing

    def test_ifn_genes_rank_at_top(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb, metric="signal_to_noise")
        # The IFN genes are up in group A, so they should occupy the top ranks
        top = set(ranked.head(len(IFN_SET)).index)
        assert top == set(IFN_SET)

    def test_log2fc_metric(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb, metric="log2fc")
        # IFN genes ~8x up → positive log2FC near 3
        assert ranked["MX1"] > 1.0

    def test_log2fc_stat_metric(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb, metric="log2fc_stat")
        assert isinstance(ranked, pd.Series)
        assert len(ranked) == exp.n_genes

    def test_unknown_metric_raises(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.raises(ValueError, match="Unknown metric"):
            rank_genes(exp.raw_counts, ga, gb, metric="nonsense")


# ---------------------------------------------------------------------------
# _running_enrichment_score
# ---------------------------------------------------------------------------


class TestRunningEnrichmentScore:
    def test_es_within_bounds(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        es, running, le = _running_enrichment_score(ranked, IFN_SET)
        assert -1.0 <= es <= 1.0
        assert len(running) == exp.n_genes

    def test_positive_es_for_top_ranked_set(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        es, _, le = _running_enrichment_score(ranked, IFN_SET)
        assert es > 0
        # Leading edge is a subset of the gene set
        assert set(le).issubset(set(IFN_SET))

    def test_empty_overlap_returns_zero(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        es, running, le = _running_enrichment_score(ranked, ["NOT_A_GENE"])
        assert es == 0.0
        assert le == []


# ---------------------------------------------------------------------------
# _compute_pvalue_and_nes
# ---------------------------------------------------------------------------


class TestComputePvalueAndNes:
    def test_positive_es(self):
        null = np.array([0.1, 0.2, -0.1, -0.3, 0.05])
        pval, nes = _compute_pvalue_and_nes(0.5, null)
        assert 0.0 < pval <= 1.0
        assert nes > 0

    def test_negative_es(self):
        null = np.array([0.1, 0.2, -0.1, -0.3, 0.05])
        pval, nes = _compute_pvalue_and_nes(-0.5, null)
        assert 0.0 < pval <= 1.0
        assert nes < 0


# ---------------------------------------------------------------------------
# gsea
# ---------------------------------------------------------------------------


class TestGsea:
    def test_returns_expected_columns(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning):
            res = gsea(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                       group_a=ga, group_b=gb)
        expected = {"gene_set", "es", "nes", "pvalue", "padj",
                    "n_genes", "n_overlap", "leading_edge"}
        assert expected.issubset(set(res.columns))

    def test_positive_es_for_ifn(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning):
            res = gsea(exp, gene_sets={"IFN_JAKSTAT": get_gene_set("IFN_JAKSTAT")},
                       group_a=ga, group_b=gb)
        row = res[res["gene_set"] == "IFN_JAKSTAT"].iloc[0]
        assert row["es"] > 0
        assert row["n_overlap"] >= 5

    def test_sorted_by_pvalue(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning):
            res = gsea(
                exp,
                gene_sets={"IFN_JAKSTAT": IFN_SET, "OTHER": ["TNF", "IL6", "CD3D", "CD4", "GAPDH_X"]},
                group_a=ga, group_b=gb,
            )
        pvals = res["pvalue"].values
        assert all(pvals[i] <= pvals[i + 1] for i in range(len(pvals) - 1))

    def test_stores_on_experiment(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        assert exp.gsea_results is None
        with pytest.warns(UserWarning):
            res = gsea(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                       group_a=ga, group_b=gb)
        assert exp.gsea_results is res

    def test_store_false(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning):
            gsea(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                 group_a=ga, group_b=gb, store=False)
        assert exp.gsea_results is None

    def test_no_overlap_returns_empty(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning, match="No gene sets passed"):
            res = gsea(exp, gene_sets={"X": ["FOO", "BAR", "BAZ"]},
                       group_a=ga, group_b=gb)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# gsea_preranked
# ---------------------------------------------------------------------------


class TestGseaPreranked:
    def test_runs_on_ranked_series(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        res = gsea_preranked(ranked, {"IFN_JAKSTAT": IFN_SET}, n_perm=100)
        assert len(res) == 1
        assert "es" in res.columns
        assert "padj" in res.columns

    def test_sorts_unsorted_input(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        shuffled = ranked.sample(frac=1.0, random_state=0)
        res = gsea_preranked(shuffled, {"IFN_JAKSTAT": IFN_SET}, n_perm=100)
        assert res.iloc[0]["es"] > 0


# ---------------------------------------------------------------------------
# competitive_test
# ---------------------------------------------------------------------------


class TestCompetitiveTest:
    def test_returns_expected_columns(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = competitive_test(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                               group_a=ga, group_b=gb)
        expected = {"gene_set", "direction", "stat", "pvalue", "padj", "n_overlap"}
        assert expected.issubset(set(res.columns))

    def test_ifn_direction_up(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = competitive_test(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                               group_a=ga, group_b=gb)
        row = res[res["gene_set"] == "IFN_JAKSTAT"].iloc[0]
        assert row["direction"] == "up"

    def test_small_set_skipped(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = competitive_test(exp, gene_sets={"tiny": ["MX1", "IFIT1"]},
                               group_a=ga, group_b=gb, min_set_size=5)
        assert len(res) == 0


# ---------------------------------------------------------------------------
# self_contained_test
# ---------------------------------------------------------------------------


class TestSelfContainedTest:
    def test_returns_expected_columns(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = self_contained_test(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                                  group_a=ga, group_b=gb)
        expected = {"gene_set", "mean_diff", "pvalue", "padj", "n_overlap"}
        assert expected.issubset(set(res.columns))

    def test_ifn_positive_mean_diff(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = self_contained_test(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                                  group_a=ga, group_b=gb)
        row = res[res["gene_set"] == "IFN_JAKSTAT"].iloc[0]
        # IFN up in group A → group A score higher → mean_diff > 0
        assert row["mean_diff"] > 0

    def test_pvalue_bounded(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        res = self_contained_test(exp, gene_sets={"IFN_JAKSTAT": IFN_SET},
                                  group_a=ga, group_b=gb)
        assert (res["pvalue"] > 0).all()
        assert (res["pvalue"] <= 1.0).all()
