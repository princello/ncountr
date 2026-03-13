"""Export utilities for writing results to disk."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Union

import pandas as pd

from ncountr.experiment import NanostringExperiment


def export_counts(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> dict[str, Path]:
    """Write raw and normalized count matrices to CSV.

    Returns a dict mapping description to output path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    raw_path = output_dir / f"{prefix}_raw_counts.csv"
    combined = pd.concat([experiment.raw_counts, experiment.hk_counts])
    combined.to_csv(raw_path)
    paths["raw_counts"] = raw_path

    if experiment.normalized is not None:
        norm_path = output_dir / f"{prefix}_normalized_counts.csv"
        experiment.normalized.to_csv(norm_path)
        paths["normalized_counts"] = norm_path

    return paths


def export_qc(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> Path | None:
    """Write QC results to CSV."""
    if experiment.qc_results is None:
        return None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{prefix}_qc_results.csv"
    experiment.qc_results.to_csv(path)
    return path


def export_de(
    experiment: NanostringExperiment,
    output_dir: Union[str, Path],
    *,
    prefix: str = "nanostring",
) -> Path | None:
    """Write DE results to CSV."""
    if experiment.de_results is None:
        return None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{prefix}_de_results.csv"
    experiment.de_results.to_csv(path, index=False)
    return path
