"""Shared fixtures for ncountr test suite."""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ncountr.experiment import NanostringExperiment


# ---------------------------------------------------------------------------
# RCC file generation helpers
# ---------------------------------------------------------------------------

_POS_CONTROLS = [
    ("POS_A(128)", 128),
    ("POS_B(32)", 32),
    ("POS_C(8)", 8),
    ("POS_D(2)", 2),
    ("POS_E(0.5)", 0.5),
    ("POS_F(0.125)", 0.125),
]

_NEG_CONTROLS = [
    "NEG_A(0)",
    "NEG_B(0)",
    "NEG_C(0)",
    "NEG_D(0)",
    "NEG_E(0)",
    "NEG_F(0)",
    "NEG_G(0)",
    "NEG_H(0)",
]

_HOUSEKEEPING = ["ACTB", "GAPDH", "RPL19"]

_ENDOGENOUS = [
    "MX1",
    "IFIT1",
    "ISG15",
    "OAS1",
    "STAT1",
    "TNF",
    "IL6",
    "CXCL10",
    "CD3D",
    "CD4",
]


def _rcc_content(
    sample_id: str,
    *,
    fov_count: int = 280,
    fov_counted: int = 275,
    binding_density: float = 0.85,
    cartridge_id: str = "TestCart",
    pos_base: int = 25000,
    neg_base: int = 8,
    hk_base: int = 5000,
    endo_counts: dict[str, int] | None = None,
    scale: float = 1.0,
    seed: int | None = None,
) -> str:
    """Generate the text content of a synthetic RCC file.

    Parameters
    ----------
    sample_id : str
        Value for the ``ID`` field in ``<Sample_Attributes>``.
    pos_base : int
        Base count for the highest positive control (POS_A).  Lower
        controls are scaled proportionally to their nominal concentration.
    neg_base : int
        Median negative control count.
    hk_base : int
        Housekeeping gene count (before jitter).
    endo_counts : dict, optional
        Override specific endogenous gene counts.
    scale : float
        Global multiplicative scale applied to all counts (simulates
        technical variation).
    seed : int, optional
        Random seed for reproducible jitter.
    """
    rng = np.random.default_rng(seed)

    lines: list[str] = []

    # Header
    lines.append("<Header>")
    lines.append("FileVersion,1.7")
    lines.append("SoftwareVersion,4.0.0.3")
    lines.append("</Header>")

    # Sample attributes
    lines.append("<Sample_Attributes>")
    lines.append(f"ID,{sample_id}")
    lines.append("</Sample_Attributes>")

    # Lane attributes
    lines.append("<Lane_Attributes>")
    lines.append(f"FovCount,{fov_count}")
    lines.append(f"FovCounted,{fov_counted}")
    lines.append(f"BindingDensity,{binding_density}")
    lines.append(f"CartridgeID,{cartridge_id}")
    lines.append("</Lane_Attributes>")

    # Code summary
    lines.append("<Code_Summary>")
    lines.append("CodeClass,Name,Accession,Count")

    # Positive controls — counts proportional to concentration
    max_conc = _POS_CONTROLS[0][1]
    for name, conc in _POS_CONTROLS:
        expected = int(pos_base * (conc / max_conc) * scale)
        jitter = int(rng.normal(0, max(expected * 0.02, 1)))
        count = max(expected + jitter, 1)
        lines.append(f"Positive,{name},ERCC_ctrl,{count}")

    # Negative controls
    for name in _NEG_CONTROLS:
        count = max(int(neg_base * scale + rng.normal(0, 3)), 0)
        lines.append(f"Negative,{name},ERCC_ctrl,{count}")

    # Housekeeping
    for name in _HOUSEKEEPING:
        count = max(int(hk_base * scale + rng.normal(0, hk_base * 0.05)), 1)
        lines.append(f"Housekeeping,{name},NM_000000.0,{count}")

    # Endogenous
    default_endo = {g: 150 for g in _ENDOGENOUS}
    if endo_counts:
        default_endo.update(endo_counts)
    for name in _ENDOGENOUS:
        base = default_endo[name]
        count = max(int(base * scale + rng.normal(0, max(base * 0.05, 1))), 0)
        lines.append(f"Endogenous,{name},NM_000000.0,{count}")

    lines.append("</Code_Summary>")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def rcc_dir(tmp_path: Path) -> Path:
    """Create a temporary directory with 3 synthetic RCC files.

    Samples: Sample_1, Sample_2, Sample_3 (IDs: 1, 2, 3).
    All have similar counts so normalization factors should be close to 1.
    """
    for i in range(1, 4):
        content = _rcc_content(
            f"Sample_{i}",
            seed=42 + i,
            scale=0.95 + 0.05 * i,
        )
        (tmp_path / f"sample{i}.RCC").write_text(content)
    return tmp_path


@pytest.fixture()
def rcc_dir_with_diff(tmp_path: Path) -> Path:
    """RCC dir where sample 3 has MX1 strongly up-regulated (10x).

    Useful for DE testing.
    """
    for i in range(1, 4):
        endo = None
        if i == 3:
            endo = {"MX1": 1500}
        content = _rcc_content(
            f"Sample_{i}",
            seed=42 + i,
            endo_counts=endo,
        )
        (tmp_path / f"sample{i}.RCC").write_text(content)
    return tmp_path


@pytest.fixture()
def experiment(rcc_dir: Path) -> NanostringExperiment:
    """A NanostringExperiment built from the 3-sample synthetic RCC files."""
    from ncountr.io.rcc import read_rcc

    return read_rcc(rcc_dir)


@pytest.fixture()
def experiment_with_diff(rcc_dir_with_diff: Path) -> NanostringExperiment:
    """Experiment where sample 3 has strongly elevated MX1."""
    from ncountr.io.rcc import read_rcc

    return read_rcc(rcc_dir_with_diff)


@pytest.fixture()
def experiment_with_meta(rcc_dir: Path) -> NanostringExperiment:
    """Experiment with sample metadata including group assignments."""
    from ncountr.io.rcc import read_rcc

    meta = {
        "1": {"group": "treated"},
        "2": {"group": "treated"},
        "3": {"group": "control"},
    }
    return read_rcc(rcc_dir, sample_meta=meta)


@pytest.fixture()
def simple_experiment() -> NanostringExperiment:
    """Minimal hand-crafted experiment with fully deterministic values.

    3 samples (S1, S2, S3); 4 endogenous genes; 3 pos controls;
    2 neg controls; 2 housekeeping genes.
    """
    raw = pd.DataFrame(
        {
            "S1": [100, 200, 50, 300],
            "S2": [110, 190, 55, 290],
            "S3": [105, 205, 48, 310],
        },
        index=["GeneA", "GeneB", "GeneC", "GeneD"],
    )
    pos = pd.DataFrame(
        {
            "S1": [25000, 6250, 390],
            "S2": [26000, 6500, 406],
            "S3": [24000, 6000, 375],
        },
        index=["POS_A(128)", "POS_B(32)", "POS_C(2)"],
    )
    neg = pd.DataFrame(
        {
            "S1": [8, 10, 6],
            "S2": [9, 11, 7],
            "S3": [7, 9, 5],
        },
        index=["NEG_A(0)", "NEG_B(0)", "NEG_C(0)"],
    )
    hk = pd.DataFrame(
        {
            "S1": [5000, 4800],
            "S2": [5200, 5000],
            "S3": [4900, 4700],
        },
        index=["ACTB", "GAPDH"],
    )
    lane = pd.DataFrame(
        {
            "FovCount": [280, 280, 280],
            "FovCounted": [275, 270, 200],
            "BindingDensity": [0.85, 0.90, 0.80],
            "CartridgeID": ["Cart1", "Cart1", "Cart1"],
        },
        index=pd.Index(["S1", "S2", "S3"], name="sample"),
    )
    meta = pd.DataFrame(
        {"group": ["treated", "treated", "control"]},
        index=pd.Index(["S1", "S2", "S3"], name="sample"),
    )

    return NanostringExperiment(
        raw_counts=raw,
        pos_counts=pos,
        neg_counts=neg,
        hk_counts=hk,
        sample_meta=meta,
        lane_info=lane,
    )
