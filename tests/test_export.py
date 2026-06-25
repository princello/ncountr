"""Tests for ncountr.io.export — CSV export and AnnData conversion."""

from __future__ import annotations

from pathlib import Path

import pytest

from ncountr.core.de import de
from ncountr.core.normalize import normalize
from ncountr.core.qc import qc
from ncountr.io.export import to_anndata, export_counts, export_qc, export_de


# ---------------------------------------------------------------------------
# export_counts / export_qc / export_de
# ---------------------------------------------------------------------------


class TestExportCounts:
    def test_writes_raw(self, simple_experiment, tmp_path: Path):
        paths = export_counts(simple_experiment, tmp_path)
        assert paths["raw_counts"].exists()
        # raw + housekeeping rows are combined
        text = paths["raw_counts"].read_text()
        assert "GeneA" in text and "ACTB" in text

    def test_no_normalized_key_when_absent(self, simple_experiment, tmp_path: Path):
        paths = export_counts(simple_experiment, tmp_path)
        assert "normalized_counts" not in paths

    def test_writes_normalized_when_present(self, simple_experiment, tmp_path: Path):
        normalize(simple_experiment)
        paths = export_counts(simple_experiment, tmp_path)
        assert paths["normalized_counts"].exists()

    def test_custom_prefix(self, simple_experiment, tmp_path: Path):
        paths = export_counts(simple_experiment, tmp_path, prefix="myrun")
        assert paths["raw_counts"].name == "myrun_raw_counts.csv"


class TestExportQc:
    def test_none_when_no_qc(self, simple_experiment, tmp_path: Path):
        assert export_qc(simple_experiment, tmp_path) is None

    def test_writes_when_present(self, simple_experiment, tmp_path: Path):
        qc(simple_experiment)
        path = export_qc(simple_experiment, tmp_path)
        assert path is not None and path.exists()


class TestExportDe:
    def test_none_when_no_de(self, simple_experiment, tmp_path: Path):
        assert export_de(simple_experiment, tmp_path) is None

    def test_writes_when_present(self, simple_experiment, tmp_path: Path):
        de(simple_experiment, group_a=["S1", "S2"], group_b=["S3"])
        path = export_de(simple_experiment, tmp_path)
        assert path is not None and path.exists()


# ---------------------------------------------------------------------------
# to_anndata (requires anndata)
# ---------------------------------------------------------------------------


class TestToAnnData:
    def test_returns_anndata_samples_x_genes(self, simple_experiment):
        anndata = pytest.importorskip("anndata")
        adata = to_anndata(simple_experiment)
        assert isinstance(adata, anndata.AnnData)
        # 3 samples x (4 endogenous + 2 housekeeping) genes
        assert adata.shape == (3, 6)

    def test_raw_layer_present(self, simple_experiment):
        pytest.importorskip("anndata")
        adata = to_anndata(simple_experiment)
        assert "raw" in adata.layers
        assert adata.layers["raw"].shape == (3, 6)

    def test_housekeeping_flag_in_var(self, simple_experiment):
        pytest.importorskip("anndata")
        adata = to_anndata(simple_experiment)
        assert "housekeeping" in adata.var.columns
        assert adata.var["housekeeping"].sum() == 2  # ACTB, GAPDH

    def test_obs_carries_metadata_and_qc(self, simple_experiment):
        pytest.importorskip("anndata")
        qc(simple_experiment)
        adata = to_anndata(simple_experiment)
        assert "group" in adata.obs.columns
        assert any(c.startswith("qc_") for c in adata.obs.columns)
        assert any(c.startswith("lane_") for c in adata.obs.columns)

    def test_controls_in_uns(self, simple_experiment):
        pytest.importorskip("anndata")
        adata = to_anndata(simple_experiment)
        assert "pos_counts" in adata.uns
        assert "neg_counts" in adata.uns

    def test_uses_normalized_when_available(self, simple_experiment):
        pytest.importorskip("anndata")
        normalize(simple_experiment)
        adata = to_anndata(simple_experiment)
        # Still samples x full gene set
        assert adata.shape == (3, 6)
