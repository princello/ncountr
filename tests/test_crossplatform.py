"""Tests for ncountr.crossplatform — cross-platform validation utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ncountr.crossplatform.correlation import (
    per_sample_correlation,
    per_gene_correlation,
)
from ncountr.crossplatform.concordance import de_concordance, concordance_summary
from ncountr.crossplatform.composition import marker_composition_proxy
from ncountr.crossplatform.crossreactivity import assess_crossreactivity
from ncountr.crossplatform.loaders import load_expression


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _matched_matrices(seed: int = 1):
    """Two genes-x-samples matrices that are strongly correlated."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(8)]
    samples = [f"S{i}" for i in range(6)]
    base = rng.uniform(1, 100, size=(len(genes), len(samples)))
    nano = pd.DataFrame(base, index=genes, columns=samples)
    external = pd.DataFrame(base * 1.5 + rng.normal(0, 1, base.shape),
                            index=genes, columns=samples)
    return nano, external


# ---------------------------------------------------------------------------
# correlation
# ---------------------------------------------------------------------------


class TestPerSampleCorrelation:
    def test_columns_and_rows(self):
        nano, ext = _matched_matrices()
        res = per_sample_correlation(nano, ext)
        assert set(res.columns) == {"sample", "r", "pvalue"}
        assert len(res) == 6

    def test_high_correlation(self):
        nano, ext = _matched_matrices()
        res = per_sample_correlation(nano, ext)
        assert res["r"].mean() > 0.5

    def test_pearson_method(self):
        nano, ext = _matched_matrices()
        res = per_sample_correlation(nano, ext, method="pearson")
        assert len(res) == 6


class TestPerGeneCorrelation:
    def test_columns(self):
        nano, ext = _matched_matrices()
        res = per_gene_correlation(nano, ext)
        assert set(res.columns) == {"gene", "r", "pvalue"}

    def test_skips_low_sample_count(self):
        nano, ext = _matched_matrices()
        # Only 2 shared samples -> below min_samples (4) -> empty
        res = per_gene_correlation(nano[["S0", "S1"]], ext[["S0", "S1"]])
        assert len(res) == 0


# ---------------------------------------------------------------------------
# concordance
# ---------------------------------------------------------------------------


class TestConcordance:
    def _de_tables(self):
        de_a = pd.DataFrame({
            "gene": ["G1", "G2", "G3", "G4"],
            "log2FC": [2.0, -1.5, 0.8, -0.2],
            "padj": [0.01, 0.02, 0.2, 0.5],
        })
        de_b = pd.DataFrame({
            "gene": ["G1", "G2", "G3", "G4"],
            "log2FC": [1.8, -1.0, -0.5, 0.3],   # G3 disagrees in direction
            "padj": [0.03, 0.04, 0.3, 0.4],
        })
        return de_a, de_b

    def test_concordance_columns(self):
        de_a, de_b = self._de_tables()
        conc = de_concordance(de_a, de_b)
        assert "same_direction" in conc.columns
        assert len(conc) == 4

    def test_direction_agreement(self):
        de_a, de_b = self._de_tables()
        conc = de_concordance(de_a, de_b).set_index("gene")
        assert conc.loc["G1", "same_direction"] == True
        assert conc.loc["G3", "same_direction"] == False

    def test_gene_flags(self):
        de_a, de_b = self._de_tables()
        conc = de_concordance(de_a, de_b, gene_flags={"G1": True}, flag_col_name="is_ifn")
        assert "is_ifn" in conc.columns
        assert conc.set_index("gene").loc["G1", "is_ifn"] == True

    def test_summary(self):
        de_a, de_b = self._de_tables()
        conc = de_concordance(de_a, de_b)
        summary = concordance_summary(conc)
        assert set(summary) == {
            "overall_rate", "n_concordant", "n_total",
            "lfc_spearman_r", "lfc_spearman_p",
        }
        assert summary["n_total"] == 4
        # G1 (+/+) and G2 (-/-) agree; G3 (+/-) and G4 (-/+) disagree
        assert summary["n_concordant"] == 2
        assert summary["overall_rate"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------


class TestMarkerComposition:
    def test_returns_correlations(self):
        samples = [f"S{i}" for i in range(6)]
        rng = np.random.default_rng(0)
        signal = rng.uniform(10, 100, size=6)
        nano = pd.DataFrame(
            {s: [signal[i], signal[i] * 1.1] for i, s in enumerate(samples)},
            index=["CD3D", "CD3E"],
        )
        proportions = pd.DataFrame(
            {"T cells": signal / signal.sum()},
            index=samples,
        )
        res = marker_composition_proxy(
            nano, proportions, {"T_cells": ["CD3D", "CD3E"]}
        )
        assert "spearman_r" in res.columns
        assert res.iloc[0]["cell_type"] == "T_cells"
        # nano signal tracks proportions perfectly -> r == 1
        assert res.iloc[0]["spearman_r"] == pytest.approx(1.0)

    def test_no_matching_column_skipped(self):
        samples = ["S0", "S1", "S2", "S3"]
        nano = pd.DataFrame({s: [10, 20] for s in samples}, index=["CD3D", "CD3E"])
        proportions = pd.DataFrame({"Bcells": [0.1, 0.2, 0.3, 0.4]}, index=samples)
        res = marker_composition_proxy(nano, proportions, {"T_cells": ["CD3D"]})
        assert len(res) == 0


# ---------------------------------------------------------------------------
# cross-reactivity
# ---------------------------------------------------------------------------


class TestCrossReactivity:
    def test_classifies_genes(self, simple_experiment):
        res = assess_crossreactivity(
            simple_experiment,
            neg_control_samples=["S3"],
            reference_samples=["S1", "S2"],
        )
        assert "xr_class" in res.columns
        assert "ref_median" in res.columns
        assert len(res) == simple_experiment.n_genes
        assert set(res["xr_class"]).issubset(
            {"reliable", "moderate", "unreliable", "no_signal"}
        )

    def test_multiple_neg_controls_mean(self, simple_experiment):
        res = assess_crossreactivity(
            simple_experiment,
            neg_control_samples=["S2", "S3"],
            reference_samples=["S1"],
        )
        assert "xr_ratio_mean" in res.columns


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------


class TestLoadExpression:
    def test_csv(self, tmp_path: Path):
        df = pd.DataFrame({"S1": [1, 2], "S2": [3, 4]}, index=["GeneA", "GeneB"])
        p = tmp_path / "expr.csv"
        df.to_csv(p)
        loaded = load_expression(p)
        assert list(loaded.index) == ["GeneA", "GeneB"]
        assert list(loaded.columns) == ["S1", "S2"]

    def test_tsv(self, tmp_path: Path):
        df = pd.DataFrame({"S1": [1, 2]}, index=["GeneA", "GeneB"])
        p = tmp_path / "expr.tsv"
        df.to_csv(p, sep="\t")
        loaded = load_expression(p)
        assert loaded.shape == (2, 1)

    def test_unknown_extension_raises(self, tmp_path: Path):
        p = tmp_path / "expr.bin"
        p.write_text("x")
        with pytest.raises(ValueError, match="Cannot infer format"):
            load_expression(p)

    def test_h5ad_full_matrix(self, tmp_path: Path):
        anndata = pytest.importorskip("anndata")
        pytest.importorskip("scanpy")
        rng = np.random.default_rng(0)
        X = rng.poisson(5, size=(6, 4)).astype(float)
        adata = anndata.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"cell{i}" for i in range(6)]),
            var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC", "GeneD"]),
        )
        p = tmp_path / "data.h5ad"
        adata.write_h5ad(p)
        loaded = load_expression(p, format="h5ad")
        # genes x cells
        assert loaded.shape == (4, 6)

    def test_h5ad_pseudobulk(self, tmp_path: Path):
        anndata = pytest.importorskip("anndata")
        pytest.importorskip("scanpy")
        rng = np.random.default_rng(1)
        X = rng.poisson(5, size=(8, 4)).astype(float)
        obs = pd.DataFrame(
            {"sample": ["P1", "P1", "P1", "P1", "P2", "P2", "P2", "P2"]},
            index=[f"cell{i}" for i in range(8)],
        )
        adata = anndata.AnnData(
            X=X, obs=obs,
            var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC", "GeneD"]),
        )
        p = tmp_path / "data.h5ad"
        adata.write_h5ad(p)
        loaded = load_expression(p, format="h5ad", pseudobulk_group_by="sample")
        assert set(loaded.columns) == {"P1", "P2"}
