"""Smoke tests for ncountr.plotting — figures render and save without error.

These are intentionally lightweight: they assert that each plotting function
returns a matplotlib Figure and writes a file, rather than checking pixels.
Run headless via the Agg backend.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from ncountr.core.de import de  # noqa: E402
from ncountr.core.gsea import rank_genes, gsea, _running_enrichment_score  # noqa: E402
from ncountr.core.pathway import score_gene_set  # noqa: E402
from ncountr.core.qc import qc  # noqa: E402
from ncountr.plotting import (  # noqa: E402
    set_style,
    plot_qc,
    plot_volcano,
    plot_pathway_scores,
    plot_heatmap,
    plot_enrichment,
    plot_gsea_dotplot,
    plot_volcano_effect,
)
from ncountr.plotting.correlation_plots import plot_correlation_scatter  # noqa: E402


IFN_SET = ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1", "CXCL10"]


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# style
# ---------------------------------------------------------------------------


class TestStyle:
    def test_set_style(self):
        set_style()
        assert plt.rcParams["axes.titlesize"] == 14

    def test_set_style_with_dpi(self):
        set_style(dpi=222)
        assert plt.rcParams["figure.dpi"] == 222


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------


class TestPlotQc:
    def test_returns_figure_and_saves(self, experiment, tmp_path: Path):
        qc(experiment)
        out = tmp_path / "qc.png"
        fig = plot_qc(experiment, output=out)
        assert isinstance(fig, plt.Figure)
        assert out.exists()


# ---------------------------------------------------------------------------
# Volcano
# ---------------------------------------------------------------------------


class TestPlotVolcano:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        de_res = de(exp, group_a=ga, group_b=gb, store=False)
        out = tmp_path / "volcano.png"
        fig = plot_volcano(de_res, output=out)
        assert isinstance(fig, plt.Figure)
        assert out.exists()

    def test_with_highlight(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        de_res = de(exp, group_a=ga, group_b=gb, store=False)
        fig = plot_volcano(de_res, highlight_genes=IFN_SET,
                           output=tmp_path / "v2.png")
        assert isinstance(fig, plt.Figure)


# ---------------------------------------------------------------------------
# Pathway scores
# ---------------------------------------------------------------------------


class TestPlotPathwayScores:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        scores = score_gene_set(exp, gene_set=IFN_SET)
        fig = plot_pathway_scores(scores, {"A": ga, "B": gb},
                                  output=tmp_path / "scores.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "scores.png").exists()


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------


class TestPlotHeatmap:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        fig = plot_heatmap(exp.raw_counts.loc[IFN_SET], output=tmp_path / "hm.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "hm.png").exists()

    def test_no_zscore(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        fig = plot_heatmap(exp.raw_counts.loc[IFN_SET], zscore=False)
        assert isinstance(fig, plt.Figure)


# ---------------------------------------------------------------------------
# GSEA plots
# ---------------------------------------------------------------------------


class TestPlotEnrichment:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        ranked = rank_genes(exp.raw_counts, ga, gb)
        es, running, _ = _running_enrichment_score(ranked, IFN_SET)
        fig = plot_enrichment(ranked, IFN_SET, running, es=es,
                              gene_set_name="IFN", output=tmp_path / "enr.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "enr.png").exists()


class TestPlotGseaDotplot:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        with pytest.warns(UserWarning):
            res = gsea(
                exp,
                gene_sets={"IFN_JAKSTAT": IFN_SET,
                           "OTHER": ["TNF", "IL6", "CD3D", "CD4", "GAPDH_X"]},
                group_a=ga, group_b=gb,
            )
        fig = plot_gsea_dotplot(res, output=tmp_path / "dot.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "dot.png").exists()


class TestPlotVolcanoEffect:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        de_res = de(exp, group_a=ga, group_b=gb, effect_size=True, store=False)
        fig = plot_volcano_effect(de_res, output=tmp_path / "ve.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "ve.png").exists()

    def test_missing_cohens_d_raises(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        de_res = de(exp, group_a=ga, group_b=gb, store=False)
        with pytest.raises(ValueError, match="cohens_d"):
            plot_volcano_effect(de_res)


# ---------------------------------------------------------------------------
# Correlation scatter
# ---------------------------------------------------------------------------


class TestPlotCorrelationScatter:
    def test_basic(self, ifn_experiment, tmp_path: Path):
        exp, ga, gb = ifn_experiment
        x = exp.raw_counts[ga[0]].values
        y = exp.raw_counts[gb[0]].values
        fig = plot_correlation_scatter(x, y, output=tmp_path / "corr.png")
        assert isinstance(fig, plt.Figure)
        assert (tmp_path / "corr.png").exists()

    def test_with_highlight_and_labels(self, ifn_experiment):
        exp, ga, gb = ifn_experiment
        x = exp.raw_counts[ga[0]].values
        y = exp.raw_counts[gb[0]].values
        labels = list(exp.raw_counts.index)
        fig = plot_correlation_scatter(x, y, labels=labels, highlight_idx=[0, 1])
        assert isinstance(fig, plt.Figure)
