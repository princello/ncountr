"""Tests for ncountr.io.rcc — RCC file parsing."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from ncountr.io.rcc import parse_rcc, read_rcc
from ncountr.experiment import NanostringExperiment


# ---------------------------------------------------------------------------
# parse_rcc
# ---------------------------------------------------------------------------


class TestParseRcc:
    """Tests for the low-level parse_rcc function."""

    def test_returns_dict_with_expected_keys(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        assert set(result.keys()) == {"sample", "lane", "counts"}

    def test_sample_attributes_parsed(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        assert "ID" in result["sample"]
        assert result["sample"]["ID"] == "Sample_1"

    def test_lane_attributes_parsed(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        assert result["lane"]["FovCount"] == "280"
        assert result["lane"]["FovCounted"] == "275"
        assert result["lane"]["BindingDensity"] == "0.85"
        assert result["lane"]["CartridgeID"] == "TestCart"

    def test_counts_contains_positive_controls(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        pos_keys = [k for k in result["counts"] if k[0] == "Positive"]
        assert len(pos_keys) == 6  # POS_A through POS_F

    def test_counts_contains_negative_controls(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        neg_keys = [k for k in result["counts"] if k[0] == "Negative"]
        assert len(neg_keys) == 8  # NEG_A through NEG_H

    def test_counts_contains_housekeeping(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        hk_keys = [k for k in result["counts"] if k[0] == "Housekeeping"]
        assert len(hk_keys) == 3  # ACTB, GAPDH, RPL19

    def test_counts_contains_endogenous(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        endo_keys = [k for k in result["counts"] if k[0] == "Endogenous"]
        assert len(endo_keys) == 10

    def test_count_values_are_ints(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        for val in result["counts"].values():
            assert isinstance(val, int)

    def test_counts_keyed_by_class_name_tuple(self, rcc_dir: Path):
        rcc_file = sorted(rcc_dir.glob("*.RCC"))[0]
        result = parse_rcc(rcc_file)
        assert ("Endogenous", "MX1") in result["counts"]
        assert ("Positive", "POS_A(128)") in result["counts"]
        assert ("Negative", "NEG_A(0)") in result["counts"]
        assert ("Housekeeping", "ACTB") in result["counts"]


# ---------------------------------------------------------------------------
# read_rcc
# ---------------------------------------------------------------------------


class TestReadRcc:
    """Tests for read_rcc which builds a NanostringExperiment."""

    def test_returns_nanostring_experiment(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert isinstance(exp, NanostringExperiment)

    def test_correct_number_of_samples(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.n_samples == 3

    def test_correct_number_of_genes(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.n_genes == 10  # 10 endogenous genes in fixtures

    def test_sample_ids_extracted(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert set(exp.samples) == {"1", "2", "3"}

    def test_raw_counts_shape(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.raw_counts.shape == (10, 3)

    def test_pos_counts_shape(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.pos_counts.shape == (6, 3)

    def test_neg_counts_shape(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.neg_counts.shape == (8, 3)

    def test_hk_counts_shape(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.hk_counts.shape == (3, 3)

    def test_lane_info_populated(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert not exp.lane_info.empty
        assert "FovCount" in exp.lane_info.columns
        assert "BindingDensity" in exp.lane_info.columns

    def test_lane_info_types(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.lane_info["FovCount"].dtype in (int, "int64")
        assert exp.lane_info["BindingDensity"].dtype == float

    def test_normalized_initially_none(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.normalized is None

    def test_qc_results_initially_none(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.qc_results is None

    def test_de_results_initially_none(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert exp.de_results is None

    def test_sample_meta_provided(self, rcc_dir: Path):
        meta = {
            "1": {"group": "A"},
            "2": {"group": "A"},
            "3": {"group": "B"},
        }
        exp = read_rcc(rcc_dir, sample_meta=meta)
        assert "group" in exp.sample_meta.columns
        assert exp.sample_meta.loc["1", "group"] == "A"

    def test_counts_are_nonnegative(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir)
        assert (exp.raw_counts >= 0).all().all()
        assert (exp.pos_counts >= 0).all().all()
        assert (exp.neg_counts >= 0).all().all()
        assert (exp.hk_counts >= 0).all().all()


# ---------------------------------------------------------------------------
# sample_id_pattern
# ---------------------------------------------------------------------------


class TestSampleIdPattern:
    """Tests for the sample_id_pattern parameter."""

    def test_default_pattern_extracts_digits(self, rcc_dir: Path):
        # Default pattern is r"(\d+)"
        exp = read_rcc(rcc_dir)
        for sid in exp.samples:
            assert sid.isdigit()

    def test_custom_pattern(self, tmp_path: Path):
        """Write an RCC where the ID is 'ABC_xyz_42' and match 'xyz_42'."""
        from tests.conftest import _rcc_content

        content = _rcc_content("ABC_xyz_42", seed=1)
        (tmp_path / "test.RCC").write_text(content)

        exp = read_rcc(tmp_path, sample_id_pattern=r"(xyz_\d+)")
        assert exp.samples == ["xyz_42"]

    def test_no_match_skips_file(self, tmp_path: Path):
        """If the pattern does not match any sample ID, the file is skipped."""
        from tests.conftest import _rcc_content

        content = _rcc_content("NoNumbers", seed=1)
        (tmp_path / "test.RCC").write_text(content)

        with pytest.raises(ValueError, match="No samples could be parsed"):
            read_rcc(tmp_path, sample_id_pattern=r"(\d+)")


# ---------------------------------------------------------------------------
# file_pattern
# ---------------------------------------------------------------------------


class TestFilePattern:
    """Tests for the file_pattern parameter."""

    def test_default_pattern_matches_rcc(self, rcc_dir: Path):
        exp = read_rcc(rcc_dir, file_pattern="*.RCC")
        assert exp.n_samples == 3

    def test_custom_pattern_filters(self, rcc_dir: Path):
        # Only match sample1.RCC
        exp = read_rcc(rcc_dir, file_pattern="sample1.RCC")
        assert exp.n_samples == 1

    def test_no_files_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="No files matching"):
            read_rcc(tmp_path, file_pattern="*.RCC")

    def test_multiple_dirs(self, tmp_path: Path):
        """read_rcc can accept a list of directories."""
        from tests.conftest import _rcc_content

        dir_a = tmp_path / "a"
        dir_b = tmp_path / "b"
        dir_a.mkdir()
        dir_b.mkdir()
        (dir_a / "s1.RCC").write_text(_rcc_content("Sample_10", seed=1))
        (dir_b / "s2.RCC").write_text(_rcc_content("Sample_20", seed=2))

        exp = read_rcc([dir_a, dir_b])
        assert exp.n_samples == 2
        assert set(exp.samples) == {"10", "20"}
