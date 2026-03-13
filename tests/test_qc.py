"""Tests for ncountr.core.qc — quality control checks."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ncountr.core.qc import qc, _fov_ratio, _pos_ctrl_r2, _neg_background
from ncountr.experiment import NanostringExperiment


# ---------------------------------------------------------------------------
# FOV ratio
# ---------------------------------------------------------------------------


class TestFovRatio:
    """Tests for FOV ratio computation."""

    def test_perfect_fov(self):
        lane = pd.DataFrame(
            {"FovCount": [280], "FovCounted": [280]},
            index=pd.Index(["S1"], name="sample"),
        )
        result = _fov_ratio(lane)
        assert result["S1"] == pytest.approx(1.0)

    def test_partial_fov(self):
        lane = pd.DataFrame(
            {"FovCount": [280, 280], "FovCounted": [210, 140]},
            index=pd.Index(["S1", "S2"], name="sample"),
        )
        result = _fov_ratio(lane)
        assert result["S1"] == pytest.approx(0.75)
        assert result["S2"] == pytest.approx(0.5)

    def test_zero_fov_count_returns_nan(self):
        lane = pd.DataFrame(
            {"FovCount": [0], "FovCounted": [0]},
            index=pd.Index(["S1"], name="sample"),
        )
        result = _fov_ratio(lane)
        assert np.isnan(result["S1"])

    def test_fov_ratio_in_qc_results(self, simple_experiment):
        qc_df = qc(simple_experiment)
        assert "FovRatio" in qc_df.columns
        assert qc_df.loc["S1", "FovRatio"] == pytest.approx(275 / 280)

    def test_fov_pass_flag(self, simple_experiment):
        qc_df = qc(simple_experiment, fov_ratio_threshold=0.75)
        # S1: 275/280=0.982 PASS, S2: 270/280=0.964 PASS, S3: 200/280=0.714 FAIL
        assert qc_df.loc["S1", "FovPass"] is True or qc_df.loc["S1", "FovPass"] == True
        assert qc_df.loc["S2", "FovPass"] is True or qc_df.loc["S2", "FovPass"] == True
        assert qc_df.loc["S3", "FovPass"] is False or qc_df.loc["S3", "FovPass"] == False


# ---------------------------------------------------------------------------
# Positive control R-squared
# ---------------------------------------------------------------------------


class TestPosCtrlR2:
    """Tests for positive control linearity (R-squared)."""

    def test_perfect_linearity(self):
        """If counts are perfectly proportional to concentration, R2 ~ 1."""
        pos = pd.DataFrame(
            {"S1": [12800, 3200, 800, 200, 50, 12]},
            index=[
                "POS_A(128)",
                "POS_B(32)",
                "POS_C(8)",
                "POS_D(2)",
                "POS_E(0.5)",
                "POS_F(0.125)",
            ],
        )
        result = _pos_ctrl_r2(pos)
        assert result["S1"] == pytest.approx(1.0, abs=0.001)

    def test_good_r2_above_threshold(self, simple_experiment):
        """Synthetic positive controls should have high R2."""
        result = _pos_ctrl_r2(simple_experiment.pos_counts)
        for sid in simple_experiment.samples:
            assert result[sid] > 0.95

    def test_r2_stored_in_qc(self, simple_experiment):
        qc_df = qc(simple_experiment)
        assert "PosCtrl_R2" in qc_df.columns
        for sid in simple_experiment.samples:
            assert 0.0 <= qc_df.loc[sid, "PosCtrl_R2"] <= 1.0

    def test_pos_ctrl_pass_flag(self, simple_experiment):
        qc_df = qc(simple_experiment, pos_r2_threshold=0.95)
        assert "PosCtrlPass" in qc_df.columns
        # Our synthetic positive controls have good linearity
        for sid in ["S1", "S2", "S3"]:
            r2 = qc_df.loc[sid, "PosCtrl_R2"]
            expected_pass = r2 > 0.95
            assert qc_df.loc[sid, "PosCtrlPass"] == expected_pass

    def test_too_few_controls_returns_nan(self):
        """With fewer than 3 valid positive controls, R2 should be NaN."""
        pos = pd.DataFrame(
            {"S1": [100, 50]},
            index=["POS_A(128)", "POS_B(32)"],
        )
        result = _pos_ctrl_r2(pos)
        assert np.isnan(result["S1"])


# ---------------------------------------------------------------------------
# Negative background
# ---------------------------------------------------------------------------


class TestNegBackground:
    """Tests for negative control background threshold."""

    def test_known_values(self):
        neg = pd.DataFrame(
            {"S1": [10, 10, 10, 10]},
            index=["NEG_A(0)", "NEG_B(0)", "NEG_C(0)", "NEG_D(0)"],
        )
        # All values the same => mean=10, std=0 => threshold=10
        result = _neg_background(neg, n_sd=2.0)
        assert result["S1"] == pytest.approx(10.0)

    def test_with_variation(self):
        neg = pd.DataFrame(
            {"S1": [8, 10, 6, 12]},
            index=["NEG_A(0)", "NEG_B(0)", "NEG_C(0)", "NEG_D(0)"],
        )
        vals = np.array([8, 10, 6, 12], dtype=float)
        expected = vals.mean() + 2.0 * vals.std()
        result = _neg_background(neg, n_sd=2.0)
        assert result["S1"] == pytest.approx(expected)

    def test_custom_n_sd(self):
        neg = pd.DataFrame(
            {"S1": [8, 10, 6, 12]},
            index=["NEG_A(0)", "NEG_B(0)", "NEG_C(0)", "NEG_D(0)"],
        )
        vals = np.array([8, 10, 6, 12], dtype=float)
        expected = vals.mean() + 3.0 * vals.std()
        result = _neg_background(neg, n_sd=3.0)
        assert result["S1"] == pytest.approx(expected)

    def test_background_in_qc(self, simple_experiment):
        qc_df = qc(simple_experiment)
        assert "NegBackground" in qc_df.columns
        for sid in simple_experiment.samples:
            assert qc_df.loc[sid, "NegBackground"] > 0

    def test_background_across_samples(self, simple_experiment):
        """Different samples should get different background values."""
        qc_df = qc(simple_experiment)
        bg = qc_df["NegBackground"]
        # Not all identical (samples have different neg counts)
        assert bg.nunique() > 1


# ---------------------------------------------------------------------------
# Overall QC function
# ---------------------------------------------------------------------------


class TestQcFunction:
    """Tests for the top-level qc() function."""

    def test_returns_dataframe(self, simple_experiment):
        result = qc(simple_experiment)
        assert isinstance(result, pd.DataFrame)

    def test_index_matches_samples(self, simple_experiment):
        result = qc(simple_experiment)
        assert list(result.index) == simple_experiment.samples

    def test_stores_on_experiment(self, simple_experiment):
        assert simple_experiment.qc_results is None
        result = qc(simple_experiment)
        assert simple_experiment.qc_results is result

    def test_binding_density_present(self, simple_experiment):
        result = qc(simple_experiment)
        assert "BindingDensity" in result.columns
        assert result.loc["S1", "BindingDensity"] == pytest.approx(0.85)

    def test_cartridge_id_present(self, simple_experiment):
        result = qc(simple_experiment)
        assert "CartridgeID" in result.columns
        assert result.loc["S1", "CartridgeID"] == "Cart1"

    def test_custom_thresholds(self, simple_experiment):
        """Adjusting thresholds changes pass/fail flags."""
        # Very permissive FOV threshold
        qc_df = qc(simple_experiment, fov_ratio_threshold=0.5)
        # S3 has 200/280 = 0.714, which passes 0.5 threshold
        assert qc_df.loc["S3", "FovPass"] == True

    def test_empty_lane_info(self):
        """QC handles missing lane_info gracefully."""
        raw = pd.DataFrame({"S1": [100]}, index=["Gene1"])
        pos = pd.DataFrame(
            {"S1": [25000, 6250, 390]},
            index=["POS_A(128)", "POS_B(32)", "POS_C(2)"],
        )
        neg = pd.DataFrame({"S1": [8, 10]}, index=["NEG_A(0)", "NEG_B(0)"])
        hk = pd.DataFrame({"S1": [5000]}, index=["ACTB"])

        exp = NanostringExperiment(
            raw_counts=raw,
            pos_counts=pos,
            neg_counts=neg,
            hk_counts=hk,
        )
        qc_df = qc(exp)
        assert np.isnan(qc_df.loc["S1", "FovRatio"])
        assert qc_df.loc["S1", "FovPass"] == True  # default when no lane info

    def test_rcc_experiment_qc(self, experiment):
        """QC works with experiment built from actual RCC files."""
        qc_df = qc(experiment)
        assert qc_df.shape[0] == 3
        assert "FovRatio" in qc_df.columns
        assert "PosCtrl_R2" in qc_df.columns
        assert "NegBackground" in qc_df.columns
