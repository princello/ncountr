"""NanostringExperiment — central data container."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class NanostringExperiment:
    """Container for a parsed Nanostring nCounter experiment.

    Attributes
    ----------
    raw_counts : pd.DataFrame
        Endogenous gene counts, genes (rows) x samples (columns).
    pos_counts : pd.DataFrame
        Positive control counts, controls (rows) x samples (columns).
    neg_counts : pd.DataFrame
        Negative control counts, controls (rows) x samples (columns).
    hk_counts : pd.DataFrame
        Housekeeping gene counts, genes (rows) x samples (columns).
    sample_meta : pd.DataFrame
        Per-sample metadata (index = sample ID).
    lane_info : pd.DataFrame
        Per-sample lane attributes (FovCount, FovCounted, BindingDensity, etc.).
    normalized : pd.DataFrame | None
        Normalized count matrix (set after calling ``normalize``).
    qc_results : pd.DataFrame | None
        QC results per sample (set after calling ``qc``).
    de_results : pd.DataFrame | None
        DE results (set after calling ``de``).
    samples : list[str]
        Ordered sample IDs.
    """

    raw_counts: pd.DataFrame
    pos_counts: pd.DataFrame
    neg_counts: pd.DataFrame
    hk_counts: pd.DataFrame
    sample_meta: pd.DataFrame = field(default_factory=pd.DataFrame)
    lane_info: pd.DataFrame = field(default_factory=pd.DataFrame)
    normalized: Optional[pd.DataFrame] = None
    qc_results: Optional[pd.DataFrame] = None
    de_results: Optional[pd.DataFrame] = None
    gsea_results: Optional[pd.DataFrame] = None

    @property
    def samples(self) -> list[str]:
        """Return ordered sample IDs from raw_counts columns."""
        return list(self.raw_counts.columns)

    @property
    def genes(self) -> list[str]:
        """Return gene names from raw_counts index."""
        return list(self.raw_counts.index)

    @property
    def n_samples(self) -> int:
        return self.raw_counts.shape[1]

    @property
    def n_genes(self) -> int:
        return self.raw_counts.shape[0]

    def __repr__(self) -> str:
        norm_str = "yes" if self.normalized is not None else "no"
        qc_str = "yes" if self.qc_results is not None else "no"
        return (
            f"NanostringExperiment("
            f"{self.n_genes} genes x {self.n_samples} samples, "
            f"normalized={norm_str}, qc={qc_str})"
        )
