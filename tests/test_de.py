"""Tests for ncountr.core.de — differential expression analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ncountr.core.de import de
from ncountr.core.normalize import normalize
from ncountr.experiment import NanostringExperiment


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_de_experiment(
    n_genes: int = 10,
    n_per_group: int = 3,
    diff_gene_idx: int = 0,
    diff_magnitude: float = 10.0,
    seed: int = 42,
) -> tuple[NanostringExperiment, list[str], list[str]]:
    """Create an experiment with known differential expression.

    Gene at `diff_gene_idx` will have mean counts of
    ``100 * diff_magnitude`` in group A and ``100`` in group B.
    """
    rng = np.random.default_rng(seed)

    group_a = [f"A{i}" for i in range(1, n_per_group + 1)]
    group_b = [f"B{i}" for i in range(1, n_per_group + 1)]
    samples = group_a + group_b
    genes = [f"Gene{i}" for i in range(n_genes)]

    # Base counts ~100 for all genes/samples
    data = rng.poisson(100, size=(n_genes, len(samples))).astype(float)

    # Make the diff gene strongly different in group A
    for j in range(n_per_group):
        data[diff_gene_idx, j] = rng.poisson(100 * diff_magnitude)

    raw = pd.DataFrame(data, index=genes, columns=samples)

    # Minimal controls (not used in DE but required for experiment)
    pos = pd.DataFrame(
        {s: [25000, 6250, 390] for s in samples},
        index=["POS_A(128)", "POS_B(32)", "POS_C(2)"],
    )
    neg = pd.DataFrame(
        {s: [8, 10] for s in samples},
        index=["NEG_A(0)", "NEG_B(0)"],
    )
    hk = pd.DataFrame(
        {s: [5000] for s in samples},
        index=["ACTB"],
    )

    exp = NanostringExperiment(
        raw_counts=raw,
        pos_counts=pos,
        neg_counts=neg,
        hk_counts=hk,
    )

    return exp, group_a, group_b


# ---------------------------------------------------------------------------
# Basic DE functionality
# ---------------------------------------------------------------------------


class TestDeBasic:
    """Basic tests for the de() function."""

    def test_returns_dataframe(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        expected_cols = {"gene", "log2FC", "mean_a", "mean_b", "pvalue", "padj"}
        assert set(result.columns) == expected_cols

    def test_one_row_per_gene(self):
        exp, ga, gb = _make_de_experiment(n_genes=5)
        result = de(exp, group_a=ga, group_b=gb)
        assert len(result) == 5

    def test_sorted_by_pvalue(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        pvals = result["pvalue"].values
        assert all(pvals[i] <= pvals[i + 1] for i in range(len(pvals) - 1))

    def test_stored_on_experiment(self):
        exp, ga, gb = _make_de_experiment()
        assert exp.de_results is None
        result = de(exp, group_a=ga, group_b=gb)
        assert exp.de_results is result

    def test_store_false(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb, store=False)
        assert exp.de_results is None
        assert result is not None


# ---------------------------------------------------------------------------
# Known differential expression
# ---------------------------------------------------------------------------


class TestDeKnownDiff:
    """Tests with known differentially expressed genes."""

    def test_diff_gene_has_lowest_pvalue(self):
        """Gene0 (strongly up in group A) should be the top hit."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=10.0)
        result = de(exp, group_a=ga, group_b=gb)
        top_gene = result.iloc[0]["gene"]
        assert top_gene == "Gene0"

    def test_diff_gene_positive_log2fc(self):
        """Gene0 is higher in group A, so log2FC should be positive."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=10.0)
        result = de(exp, group_a=ga, group_b=gb)
        gene0_row = result[result["gene"] == "Gene0"].iloc[0]
        assert gene0_row["log2FC"] > 0

    def test_diff_gene_significant_pvalue(self):
        """With 10x difference, p-value should be small."""
        exp, ga, gb = _make_de_experiment(
            diff_magnitude=10.0, n_per_group=5, seed=123
        )
        result = de(exp, group_a=ga, group_b=gb)
        gene0_row = result[result["gene"] == "Gene0"].iloc[0]
        assert gene0_row["pvalue"] < 0.05

    def test_nondiff_genes_high_pvalue(self):
        """Non-differential genes should generally have higher p-values."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=10.0, n_per_group=5)
        result = de(exp, group_a=ga, group_b=gb)
        # Gene0 should have much lower pvalue than the median of others
        gene0_p = result[result["gene"] == "Gene0"]["pvalue"].values[0]
        others_median = result[result["gene"] != "Gene0"]["pvalue"].median()
        assert gene0_p < others_median

    def test_mean_a_higher_for_diff_gene(self):
        exp, ga, gb = _make_de_experiment(diff_magnitude=5.0)
        result = de(exp, group_a=ga, group_b=gb)
        gene0_row = result[result["gene"] == "Gene0"].iloc[0]
        assert gene0_row["mean_a"] > gene0_row["mean_b"]


# ---------------------------------------------------------------------------
# Statistical test options
# ---------------------------------------------------------------------------


class TestDeTestMethods:
    """Tests for mannwhitneyu and ttest options."""

    def test_mannwhitneyu_default(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb, test="mannwhitneyu")
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_ttest(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb, test="ttest")
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_methods_give_different_pvalues(self):
        """Mann-Whitney and t-test should generally give different p-values."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=5.0, n_per_group=5)
        result_mw = de(exp, group_a=ga, group_b=gb, test="mannwhitneyu", store=False)
        result_tt = de(exp, group_a=ga, group_b=gb, test="ttest", store=False)

        # p-values should not be identical
        mw_pvals = result_mw.set_index("gene")["pvalue"]
        tt_pvals = result_tt.set_index("gene")["pvalue"]
        assert not np.allclose(mw_pvals.values, tt_pvals.values)

    def test_constant_gene_gets_pvalue_one(self):
        """A gene with zero variance in both groups should get p=1."""
        exp, ga, gb = _make_de_experiment(n_genes=2, seed=99)
        # Manually set Gene1 to constant across all samples
        for s in ga + gb:
            exp.raw_counts.loc["Gene1", s] = 50.0
        result = de(exp, group_a=ga, group_b=gb)
        gene1_row = result[result["gene"] == "Gene1"].iloc[0]
        assert gene1_row["pvalue"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# FDR correction
# ---------------------------------------------------------------------------


class TestDeFdrCorrection:
    """Tests for multiple testing correction."""

    def test_padj_present(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        assert "padj" in result.columns

    def test_padj_geq_pvalue(self):
        """Adjusted p-values should be >= raw p-values."""
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        assert (result["padj"] >= result["pvalue"] - 1e-15).all()

    def test_padj_leq_one(self):
        exp, ga, gb = _make_de_experiment()
        result = de(exp, group_a=ga, group_b=gb)
        assert (result["padj"] <= 1.0 + 1e-15).all()

    def test_custom_correction(self):
        """Different correction methods should give different padj."""
        exp, ga, gb = _make_de_experiment(n_genes=20, n_per_group=5)
        result_bh = de(
            exp, group_a=ga, group_b=gb, correction="fdr_bh", store=False
        )
        result_bonf = de(
            exp, group_a=ga, group_b=gb, correction="bonferroni", store=False
        )
        # Bonferroni is more conservative
        bh_padj = result_bh.set_index("gene")["padj"]
        bonf_padj = result_bonf.set_index("gene")["padj"]
        # Bonferroni should be >= BH for at least some genes
        assert (bonf_padj.reindex(bh_padj.index) >= bh_padj - 1e-10).any()


# ---------------------------------------------------------------------------
# Group validation
# ---------------------------------------------------------------------------


class TestDeGroupValidation:
    """Tests for sample group validation."""

    def test_missing_group_a_sample_raises(self):
        exp, ga, gb = _make_de_experiment()
        with pytest.raises(ValueError, match="group_a samples not in counts"):
            de(exp, group_a=["NONEXISTENT"], group_b=gb)

    def test_missing_group_b_sample_raises(self):
        exp, ga, gb = _make_de_experiment()
        with pytest.raises(ValueError, match="group_b samples not in counts"):
            de(exp, group_a=ga, group_b=["NONEXISTENT"])

    def test_partial_missing_raises(self):
        exp, ga, gb = _make_de_experiment()
        with pytest.raises(ValueError, match="group_a samples not in counts"):
            de(exp, group_a=ga + ["MISSING"], group_b=gb)


# ---------------------------------------------------------------------------
# Count matrix selection
# ---------------------------------------------------------------------------


class TestDeCountSource:
    """Tests for which count matrix gets used."""

    def test_uses_raw_when_no_normalized(self):
        """When normalized is None, DE uses raw_counts."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=10.0)
        assert exp.normalized is None
        result = de(exp, group_a=ga, group_b=gb, store=False)
        assert len(result) > 0

    def test_uses_normalized_when_available(self):
        """When normalized is set, DE uses it instead of raw."""
        exp, ga, gb = _make_de_experiment(diff_magnitude=10.0)
        normalize(exp, method="pos_only")
        assert exp.normalized is not None

        result = de(exp, group_a=ga, group_b=gb, store=False)
        # Verify it ran successfully
        assert len(result) == exp.n_genes

    def test_custom_counts_parameter(self):
        """Can provide a custom count matrix."""
        exp, ga, gb = _make_de_experiment(n_genes=3)
        # Create a custom count matrix with one gene strongly different
        custom = pd.DataFrame(
            {
                "A1": [1000, 100, 100],
                "A2": [1100, 110, 95],
                "A3": [900, 90, 105],
                "B1": [100, 100, 100],
                "B2": [110, 110, 95],
                "B3": [90, 90, 105],
            },
            index=["Gene0", "Gene1", "Gene2"],
        )
        result = de(exp, group_a=ga, group_b=gb, counts=custom, store=False)
        top_gene = result.iloc[0]["gene"]
        assert top_gene == "Gene0"


# ---------------------------------------------------------------------------
# DE with RCC experiment
# ---------------------------------------------------------------------------


class TestDeFromRcc:
    """Tests using the RCC-based experiment fixture."""

    def test_de_on_rcc_experiment(self, experiment):
        """DE runs on the RCC-parsed experiment."""
        samples = experiment.samples
        # Use first 2 as group A, last as group B
        result = de(
            experiment,
            group_a=samples[:2],
            group_b=samples[2:],
        )
        assert len(result) == experiment.n_genes

    def test_de_after_normalization(self, experiment):
        """DE works after normalization."""
        normalize(experiment)
        samples = experiment.samples
        result = de(
            experiment,
            group_a=samples[:2],
            group_b=samples[2:],
        )
        assert len(result) == experiment.n_genes

    def test_de_detects_known_diff(self, experiment_with_diff):
        """MX1 should have the largest log2FC when 10x elevated in sample 3."""
        normalize(experiment_with_diff)
        samples = experiment_with_diff.samples
        # Sample 3 has elevated MX1 → group_a = [3], group_b = [1, 2]
        result = de(
            experiment_with_diff,
            group_a=[samples[2]],  # sample "3"
            group_b=samples[:2],   # samples "1", "2"
        )
        # MX1 should have the highest absolute log2FC
        max_lfc_gene = result.loc[result["log2FC"].abs().idxmax(), "gene"]
        assert max_lfc_gene == "MX1"
