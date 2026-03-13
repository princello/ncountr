"""YAML configuration parsing for ncountr pipelines."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

import yaml


@dataclass
class NcountrConfig:
    """Parsed pipeline configuration."""

    # Input
    rcc_dirs: list[str] = field(default_factory=list)
    file_pattern: str = "*.RCC"
    sample_id_pattern: str = r"(\d+)"

    # Output
    output_dir: str = "./results"
    figure_format: str = "png"
    figure_dpi: int = 200

    # Samples
    sample_meta: dict[str, dict] = field(default_factory=dict)
    group_column: str = "group"
    comparison: list[str] = field(default_factory=list)
    exclude_samples: list[str] = field(default_factory=list)
    negative_control_samples: list[str] = field(default_factory=list)

    # QC
    fov_ratio_threshold: float = 0.75
    pos_r2_threshold: float = 0.95

    # Normalization
    normalization_method: str = "pos_hk"

    # DE
    de_test: str = "mannwhitneyu"
    de_correction: str = "fdr_bh"

    # Volcano plot
    highlight_genes: list[str] | str | None = None  # gene list, or builtin set name
    highlight_label: str = "Highlighted"
    highlight_color: str = "gold"

    # Gene sets
    gene_sets: dict[str, Any] = field(default_factory=dict)

    # Cross-platform (optional)
    cross_platform_enabled: bool = False
    cross_platform: dict[str, Any] = field(default_factory=dict)


def load_config(path: Union[str, Path]) -> NcountrConfig:
    """Load a YAML config file and return an NcountrConfig.

    Parameters
    ----------
    path : str or Path
        Path to the YAML configuration file.

    Returns
    -------
    NcountrConfig
    """
    path = Path(path)
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    cfg = NcountrConfig()

    # Input
    inp = raw.get("input", {})
    if "rcc_dirs" in inp:
        dirs = inp["rcc_dirs"]
        cfg.rcc_dirs = dirs if isinstance(dirs, list) else [dirs]
    if "file_pattern" in inp:
        cfg.file_pattern = inp["file_pattern"]
    if "sample_id_pattern" in inp:
        cfg.sample_id_pattern = inp["sample_id_pattern"]

    # Output
    out = raw.get("output", {})
    if "directory" in out:
        cfg.output_dir = out["directory"]
    if "figure_format" in out:
        cfg.figure_format = out["figure_format"]
    if "figure_dpi" in out:
        cfg.figure_dpi = int(out["figure_dpi"])

    # Samples
    samp = raw.get("samples", {})
    if "metadata" in samp:
        cfg.sample_meta = dict(samp["metadata"])
    if "group_column" in samp:
        cfg.group_column = samp["group_column"]
    if "comparison" in samp:
        cfg.comparison = list(samp["comparison"])
    if "exclude" in samp:
        cfg.exclude_samples = list(samp["exclude"])
    if "negative_control_samples" in samp:
        cfg.negative_control_samples = list(samp["negative_control_samples"])

    # QC
    qc = raw.get("qc", {})
    if "fov_ratio_threshold" in qc:
        cfg.fov_ratio_threshold = float(qc["fov_ratio_threshold"])
    if "positive_control_r2_threshold" in qc:
        cfg.pos_r2_threshold = float(qc["positive_control_r2_threshold"])

    # Normalization
    norm = raw.get("normalization", {})
    if "method" in norm:
        cfg.normalization_method = norm["method"]

    # DE
    de = raw.get("de", {})
    if "test" in de:
        cfg.de_test = de["test"]
    if "correction" in de:
        cfg.de_correction = de["correction"]

    # Volcano plot options
    volcano = raw.get("volcano", {})
    if "highlight_genes" in volcano:
        cfg.highlight_genes = volcano["highlight_genes"]
    if "highlight_label" in volcano:
        cfg.highlight_label = volcano["highlight_label"]
    if "highlight_color" in volcano:
        cfg.highlight_color = volcano["highlight_color"]

    # Gene sets
    if "gene_sets" in raw:
        cfg.gene_sets = dict(raw["gene_sets"])

    # Cross-platform
    xp = raw.get("cross_platform", {})
    if xp:
        cfg.cross_platform_enabled = xp.get("enabled", False)
        cfg.cross_platform = dict(xp)

    return cfg


def generate_template() -> str:
    """Return a YAML config template string."""
    return """\
# ncountr configuration
# See https://github.com/princello/ncountr for documentation.

input:
  rcc_dirs:
    - /path/to/RCC/files
  file_pattern: '*.RCC'
  sample_id_pattern: '(\\d+)'    # regex to extract sample ID

output:
  directory: ./results
  figure_format: png
  figure_dpi: 200

samples:
  metadata:
    # Sample ID → metadata key-value pairs
    'S1': { group: treated }
    'S2': { group: treated }
    'S3': { group: control }
    'S4': { group: control }
  group_column: group
  comparison: [treated, control]   # [group_a, group_b] for DE
  exclude: []
  negative_control_samples: []     # e.g. NSG samples for cross-reactivity

qc:
  fov_ratio_threshold: 0.75
  positive_control_r2_threshold: 0.95

normalization:
  method: pos_hk    # pos_only | pos_hk | pos_hk_bg

de:
  test: mannwhitneyu    # mannwhitneyu | ttest
  correction: fdr_bh

# Volcano plot options (optional)
# volcano:
#   highlight_genes: IFN_JAKSTAT   # builtin gene set name, or list of genes
#   highlight_label: IFN/JAK-STAT genes
#   highlight_color: gold
#   # Custom gene list example:
#   # highlight_genes:
#   #   - STAT1
#   #   - STAT2
#   #   - IRF1

gene_sets:
  IFN_JAKSTAT: builtin
  # custom_set:
  #   - GENE1
  #   - GENE2

# Cross-platform validation (optional)
# cross_platform:
#   enabled: true
#   external_data:
#     path: /path/to/data.h5ad
#     format: h5ad
#     pseudobulk_group_by: Sample
#   sample_mapping:
#     'S1': 'S1'
#     'S2': 'S2'
"""
