"""Tests for ncountr.cli — the Click command-line interface."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

import ncountr
from ncountr.io import geo as geo_mod
from ncountr.cli import cli


@pytest.fixture()
def runner() -> CliRunner:
    return CliRunner()


def _write_counts_csv(path: Path) -> None:
    """A genes-x-samples count matrix with one clearly DE gene."""
    rng = np.random.default_rng(0)
    genes = [f"Gene{i}" for i in range(6)]
    samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
    data = rng.poisson(100, size=(len(genes), len(samples))).astype(int)
    data[0, :3] = rng.poisson(900, size=3)  # Gene0 up in S1-S3
    pd.DataFrame(data, index=genes, columns=samples).to_csv(path)


# ---------------------------------------------------------------------------
# init / --version
# ---------------------------------------------------------------------------


class TestSimpleCommands:
    def test_version(self, runner):
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert ncountr.__version__ in result.output

    def test_init_emits_valid_yaml(self, runner):
        result = runner.invoke(cli, ["init"])
        assert result.exit_code == 0
        parsed = yaml.safe_load(result.output)
        assert "input" in parsed and "normalization" in parsed


# ---------------------------------------------------------------------------
# parse
# ---------------------------------------------------------------------------


class TestParse:
    def test_parse_writes_counts(self, runner, rcc_dir: Path, tmp_path: Path):
        out = tmp_path / "raw_counts.csv"
        result = runner.invoke(
            cli, ["parse", "-d", str(rcc_dir), "-o", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        df = pd.read_csv(out, index_col=0)
        assert df.shape[1] == 3  # 3 samples


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------


class TestQcCommand:
    def test_qc_runs_on_csv(self, runner, tmp_path: Path):
        counts = tmp_path / "counts.csv"
        _write_counts_csv(counts)
        result = runner.invoke(cli, ["qc", "-c", str(counts)])
        assert result.exit_code == 0, result.output
        assert "genes x" in result.output


# ---------------------------------------------------------------------------
# de
# ---------------------------------------------------------------------------


class TestDeCommand:
    def test_de_writes_results(self, runner, tmp_path: Path):
        counts = tmp_path / "counts.csv"
        _write_counts_csv(counts)
        out = tmp_path / "de.csv"
        result = runner.invoke(
            cli,
            ["de", "-c", str(counts),
             "-g", "treated:S1,S2,S3", "-g", "control:S4,S5,S6",
             "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        de_df = pd.read_csv(out)
        assert "log2FC" in de_df.columns

    def test_de_requires_two_groups(self, runner, tmp_path: Path):
        counts = tmp_path / "counts.csv"
        _write_counts_csv(counts)
        result = runner.invoke(
            cli, ["de", "-c", str(counts), "-g", "only:S1,S2"]
        )
        assert result.exit_code != 0
        assert "exactly 2 groups" in result.output


# ---------------------------------------------------------------------------
# fetch-geo (mocked)
# ---------------------------------------------------------------------------


class TestFetchGeoCommand:
    def test_fetch_geo_reports_count(self, runner, tmp_path: Path, monkeypatch):
        fake_dir = tmp_path / "GSE123456"
        fake_dir.mkdir()
        (fake_dir / "a.RCC").write_text("x")
        (fake_dir / "b.RCC").write_text("y")

        def fake_fetch(accession, output_dir=".", **kwargs):
            return fake_dir

        monkeypatch.setattr(geo_mod, "fetch_geo", fake_fetch)
        result = runner.invoke(
            cli, ["fetch-geo", "GSE123456", "-o", str(tmp_path)]
        )
        assert result.exit_code == 0, result.output
        assert "Downloaded 2 RCC files" in result.output


# ---------------------------------------------------------------------------
# run (full pipeline)
# ---------------------------------------------------------------------------


class TestRunCommand:
    def test_full_pipeline(self, runner, rcc_dir: Path, tmp_path: Path):
        outdir = tmp_path / "results"
        config = {
            "input": {
                "rcc_dirs": [str(rcc_dir)],
                "file_pattern": "*.RCC",
                "sample_id_pattern": r"(\d+)",
            },
            "output": {"directory": str(outdir), "figure_format": "png"},
            "samples": {
                "metadata": {
                    "1": {"group": "treated"},
                    "2": {"group": "treated"},
                    "3": {"group": "control"},
                },
                "group_column": "group",
                "comparison": ["treated", "control"],
            },
            "normalization": {"method": "pos_hk"},
            "de": {"test": "mannwhitneyu", "correction": "fdr_bh"},
            "gene_sets": {"IFN_JAKSTAT": "builtin"},
        }
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text(yaml.safe_dump(config))

        result = runner.invoke(cli, ["run", str(cfg_path)])
        assert result.exit_code == 0, result.output
        assert (outdir / "qc_summary.png").exists()
        assert (outdir / "volcano.png").exists()
        assert (outdir / "nanostring_raw_counts.csv").exists()
        assert (outdir / "nanostring_de_results.csv").exists()

    def test_run_errors_without_rcc_dirs(self, runner, tmp_path: Path):
        cfg_path = tmp_path / "empty.yaml"
        cfg_path.write_text("normalization:\n  method: pos_hk\n")
        result = runner.invoke(cli, ["run", str(cfg_path)])
        assert result.exit_code != 0
        assert "no rcc_dirs" in result.output
