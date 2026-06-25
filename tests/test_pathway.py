"""Tests for ncountr.core.pathway — gene set scoring."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ncountr.core.pathway import score_gene_set


IFN_SET = ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1", "CXCL10"]


class TestScoreGeneSet:
    def test_returns_series_one_per_sample(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET)
        assert isinstance(scores, pd.Series)
        assert set(scores.index) == set(exp.samples)

    def test_score_name(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET)
        assert scores.name == "pathway_score"

    def test_group_a_scores_higher(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET)
        # IFN genes elevated in group A → higher z-score-mean for group A
        assert scores[ga].mean() > scores[gb].mean()

    def test_builtin_gene_set_name(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        # "IFN_JAKSTAT" resolves via datasets.get_gene_set
        scores = score_gene_set(exp, gene_set="IFN_JAKSTAT")
        assert len(scores) == exp.n_samples

    def test_ssgsea_method(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET, method="ssgsea")
        assert isinstance(scores, pd.Series)
        assert len(scores) == exp.n_samples

    def test_samples_subset(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET, samples=ga)
        assert set(scores.index) == set(ga)

    def test_no_genes_present_raises(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.raises(ValueError, match="No genes from the gene set"):
            score_gene_set(exp, gene_set=["NOPE1", "NOPE2"])

    def test_unknown_method_raises(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        with pytest.raises(ValueError, match="Unknown scoring method"):
            score_gene_set(exp, gene_set=IFN_SET, method="bogus")

    def test_custom_counts(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET, counts=exp.raw_counts)
        assert len(scores) == exp.n_samples
