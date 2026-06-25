"""Tests for ncountr.experiment.NanostringExperiment."""

from __future__ import annotations

from ncountr.core.normalize import normalize
from ncountr.core.qc import qc


class TestExperimentProperties:
    def test_samples(self, simple_experiment):
        assert simple_experiment.samples == ["S1", "S2", "S3"]

    def test_genes(self, simple_experiment):
        assert simple_experiment.genes == ["GeneA", "GeneB", "GeneC", "GeneD"]

    def test_n_samples(self, simple_experiment):
        assert simple_experiment.n_samples == 3

    def test_n_genes(self, simple_experiment):
        assert simple_experiment.n_genes == 4

    def test_result_slots_start_none(self, simple_experiment):
        assert simple_experiment.normalized is None
        assert simple_experiment.qc_results is None
        assert simple_experiment.de_results is None
        assert simple_experiment.gsea_results is None


class TestExperimentRepr:
    def test_repr_basic(self, simple_experiment):
        r = repr(simple_experiment)
        assert "NanostringExperiment" in r
        assert "4 genes x 3 samples" in r
        assert "normalized=no" in r
        assert "qc=no" in r

    def test_repr_reflects_state(self, simple_experiment):
        qc(simple_experiment)
        normalize(simple_experiment)
        r = repr(simple_experiment)
        assert "normalized=yes" in r
        assert "qc=yes" in r
