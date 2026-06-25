"""Tests for ncountr.core.de.effect_sizes and de(effect_size=True)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ncountr.core.de import de, effect_sizes
from ncountr.experiment import NanostringExperiment


def _experiment_with_signal(seed: int = 3):
    rng = np.random.default_rng(seed)
    genes = [f"Gene{i}" for i in range(5)]
    ga = ["A1", "A2", "A3", "A4"]
    gb = ["B1", "B2", "B3", "B4"]
    samples = ga + gb
    data = rng.poisson(100, size=(len(genes), len(samples))).astype(float)
    # Gene0 strongly higher in group A
    for j in range(len(ga)):
        data[0, j] = rng.poisson(800)
    raw = pd.DataFrame(data, index=genes, columns=samples)
    pos = pd.DataFrame({s: [25000, 6250, 390] for s in samples},
                       index=["POS_A(128)", "POS_B(32)", "POS_C(2)"])
    neg = pd.DataFrame({s: [8, 10] for s in samples}, index=["NEG_A(0)", "NEG_B(0)"])
    hk = pd.DataFrame({s: [5000] for s in samples}, index=["ACTB"])
    exp = NanostringExperiment(raw_counts=raw, pos_counts=pos, neg_counts=neg, hk_counts=hk)
    return exp, ga, gb


class TestEffectSizes:
    def test_expected_columns(self):
        exp, ga, gb = _experiment_with_signal()
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=200)
        expected = {"gene", "cohens_d", "cohens_d_ci_lo", "cohens_d_ci_hi", "rank_biserial"}
        assert set(es.columns) == expected

    def test_one_row_per_gene(self):
        exp, ga, gb = _experiment_with_signal()
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=200)
        assert len(es) == exp.n_genes

    def test_diff_gene_positive_cohens_d(self):
        exp, ga, gb = _experiment_with_signal()
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=300)
        row = es[es["gene"] == "Gene0"].iloc[0]
        assert row["cohens_d"] > 0

    def test_ci_brackets_ordering(self):
        exp, ga, gb = _experiment_with_signal()
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=300)
        assert (es["cohens_d_ci_lo"] <= es["cohens_d_ci_hi"]).all()

    def test_rank_biserial_in_range(self):
        exp, ga, gb = _experiment_with_signal()
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=200)
        assert (es["rank_biserial"].abs() <= 1.0 + 1e-9).all()

    def test_constant_gene_zero_d(self):
        exp, ga, gb = _experiment_with_signal()
        exp.raw_counts.loc["Gene1", :] = 50.0
        es = effect_sizes(exp.raw_counts, ga, gb, n_bootstrap=100)
        row = es[es["gene"] == "Gene1"].iloc[0]
        assert row["cohens_d"] == 0.0


class TestDeWithEffectSize:
    def test_effect_size_columns_added(self):
        exp, ga, gb = _experiment_with_signal()
        res = de(exp, group_a=ga, group_b=gb, effect_size=True, store=False)
        for col in ["cohens_d", "cohens_d_ci_lo", "cohens_d_ci_hi", "rank_biserial"]:
            assert col in res.columns

    def test_no_effect_size_columns_by_default(self):
        exp, ga, gb = _experiment_with_signal()
        res = de(exp, group_a=ga, group_b=gb, store=False)
        assert "cohens_d" not in res.columns

    def test_effect_size_merge_one_row_per_gene(self):
        exp, ga, gb = _experiment_with_signal()
        res = de(exp, group_a=ga, group_b=gb, effect_size=True, store=False)
        assert len(res) == exp.n_genes
