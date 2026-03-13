"""Tests for ncountr.config — YAML configuration handling."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from ncountr.config import load_config, generate_template, NcountrConfig


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


class TestLoadConfig:
    """Tests for loading configuration from YAML files."""

    def test_returns_ncountr_config(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("input:\n  rcc_dirs:\n    - /data/rcc\n")
        result = load_config(cfg_file)
        assert isinstance(result, NcountrConfig)

    def test_input_rcc_dirs(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "input:\n"
            "  rcc_dirs:\n"
            "    - /path/a\n"
            "    - /path/b\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.rcc_dirs == ["/path/a", "/path/b"]

    def test_input_single_rcc_dir(self, tmp_path: Path):
        """A single string for rcc_dirs gets wrapped into a list."""
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("input:\n  rcc_dirs: /single/path\n")
        cfg = load_config(cfg_file)
        assert cfg.rcc_dirs == ["/single/path"]

    def test_input_file_pattern(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "input:\n"
            "  rcc_dirs:\n"
            "    - /data\n"
            "  file_pattern: '*.rcc'\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.file_pattern == "*.rcc"

    def test_input_sample_id_pattern(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "input:\n"
            "  rcc_dirs:\n"
            "    - /data\n"
            "  sample_id_pattern: '(S\\d+)'\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.sample_id_pattern == r"(S\d+)"

    def test_output_directory(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("output:\n  directory: /output/path\n")
        cfg = load_config(cfg_file)
        assert cfg.output_dir == "/output/path"

    def test_output_figure_format(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "output:\n"
            "  figure_format: svg\n"
            "  figure_dpi: 300\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.figure_format == "svg"
        assert cfg.figure_dpi == 300

    def test_samples_metadata(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "samples:\n"
            "  metadata:\n"
            "    S1: { group: treated }\n"
            "    S2: { group: control }\n"
            "  group_column: group\n"
            "  comparison: [treated, control]\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.sample_meta == {
            "S1": {"group": "treated"},
            "S2": {"group": "control"},
        }
        assert cfg.group_column == "group"
        assert cfg.comparison == ["treated", "control"]

    def test_samples_exclude(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("samples:\n  exclude:\n    - S99\n    - S100\n")
        cfg = load_config(cfg_file)
        assert cfg.exclude_samples == ["S99", "S100"]

    def test_samples_negative_control(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "samples:\n  negative_control_samples:\n    - NSG1\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.negative_control_samples == ["NSG1"]

    def test_qc_thresholds(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "qc:\n"
            "  fov_ratio_threshold: 0.80\n"
            "  positive_control_r2_threshold: 0.99\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.fov_ratio_threshold == 0.80
        assert cfg.pos_r2_threshold == 0.99

    def test_normalization_method(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("normalization:\n  method: pos_only\n")
        cfg = load_config(cfg_file)
        assert cfg.normalization_method == "pos_only"

    def test_de_settings(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "de:\n"
            "  test: ttest\n"
            "  correction: bonferroni\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.de_test == "ttest"
        assert cfg.de_correction == "bonferroni"

    def test_gene_sets(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "gene_sets:\n"
            "  IFN_JAKSTAT: builtin\n"
            "  custom_set:\n"
            "    - GENE1\n"
            "    - GENE2\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.gene_sets["IFN_JAKSTAT"] == "builtin"
        assert cfg.gene_sets["custom_set"] == ["GENE1", "GENE2"]

    def test_cross_platform(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            "cross_platform:\n"
            "  enabled: true\n"
            "  external_data:\n"
            "    path: /data/adata.h5ad\n"
            "    format: h5ad\n"
        )
        cfg = load_config(cfg_file)
        assert cfg.cross_platform_enabled is True
        assert cfg.cross_platform["external_data"]["path"] == "/data/adata.h5ad"

    def test_empty_file_returns_defaults(self, tmp_path: Path):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("")
        cfg = load_config(cfg_file)
        assert isinstance(cfg, NcountrConfig)
        assert cfg.rcc_dirs == []
        assert cfg.normalization_method == "pos_hk"

    def test_full_config_roundtrip(self, tmp_path: Path):
        """Write a complete config and verify all fields load correctly."""
        cfg_content = """\
input:
  rcc_dirs:
    - /data/batch1
    - /data/batch2
  file_pattern: '*.RCC'
  sample_id_pattern: '(\\d+)'

output:
  directory: /results
  figure_format: pdf
  figure_dpi: 150

samples:
  metadata:
    '1': { group: infected }
    '2': { group: infected }
    '3': { group: uninfected }
  group_column: group
  comparison: [infected, uninfected]
  exclude: []
  negative_control_samples: ['697']

qc:
  fov_ratio_threshold: 0.80
  positive_control_r2_threshold: 0.98

normalization:
  method: pos_hk_bg

de:
  test: ttest
  correction: fdr_bh

gene_sets:
  IFN_JAKSTAT: builtin
"""
        cfg_file = tmp_path / "full.yaml"
        cfg_file.write_text(cfg_content)
        cfg = load_config(cfg_file)

        assert cfg.rcc_dirs == ["/data/batch1", "/data/batch2"]
        assert cfg.file_pattern == "*.RCC"
        assert cfg.output_dir == "/results"
        assert cfg.figure_format == "pdf"
        assert cfg.figure_dpi == 150
        assert cfg.group_column == "group"
        assert cfg.comparison == ["infected", "uninfected"]
        assert cfg.negative_control_samples == ["697"]
        assert cfg.fov_ratio_threshold == 0.80
        assert cfg.pos_r2_threshold == 0.98
        assert cfg.normalization_method == "pos_hk_bg"
        assert cfg.de_test == "ttest"
        assert cfg.de_correction == "fdr_bh"


# ---------------------------------------------------------------------------
# generate_template
# ---------------------------------------------------------------------------


class TestGenerateTemplate:
    """Tests for the YAML template generator."""

    def test_returns_string(self):
        result = generate_template()
        assert isinstance(result, str)

    def test_valid_yaml(self):
        """The generated template should be parseable YAML."""
        result = generate_template()
        parsed = yaml.safe_load(result)
        assert isinstance(parsed, dict)

    def test_contains_key_sections(self):
        result = generate_template()
        parsed = yaml.safe_load(result)
        assert "input" in parsed
        assert "output" in parsed
        assert "samples" in parsed
        assert "qc" in parsed
        assert "normalization" in parsed
        assert "de" in parsed

    def test_input_has_rcc_dirs(self):
        result = generate_template()
        parsed = yaml.safe_load(result)
        assert "rcc_dirs" in parsed["input"]

    def test_template_loadable_as_config(self, tmp_path: Path):
        """The template should be loadable by load_config."""
        template = generate_template()
        cfg_file = tmp_path / "template.yaml"
        cfg_file.write_text(template)
        cfg = load_config(cfg_file)
        assert isinstance(cfg, NcountrConfig)
        assert cfg.normalization_method == "pos_hk"

    def test_template_contains_comments(self):
        """The template should include helpful comments."""
        result = generate_template()
        assert "#" in result  # has comments

    def test_template_has_gene_sets(self):
        result = generate_template()
        parsed = yaml.safe_load(result)
        assert "gene_sets" in parsed


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------


class TestDefaults:
    """Tests for NcountrConfig default values."""

    def test_default_rcc_dirs(self):
        cfg = NcountrConfig()
        assert cfg.rcc_dirs == []

    def test_default_file_pattern(self):
        cfg = NcountrConfig()
        assert cfg.file_pattern == "*.RCC"

    def test_default_sample_id_pattern(self):
        cfg = NcountrConfig()
        assert cfg.sample_id_pattern == r"(\d+)"

    def test_default_output_dir(self):
        cfg = NcountrConfig()
        assert cfg.output_dir == "./results"

    def test_default_figure_format(self):
        cfg = NcountrConfig()
        assert cfg.figure_format == "png"

    def test_default_figure_dpi(self):
        cfg = NcountrConfig()
        assert cfg.figure_dpi == 200

    def test_default_fov_ratio_threshold(self):
        cfg = NcountrConfig()
        assert cfg.fov_ratio_threshold == 0.75

    def test_default_pos_r2_threshold(self):
        cfg = NcountrConfig()
        assert cfg.pos_r2_threshold == 0.95

    def test_default_normalization_method(self):
        cfg = NcountrConfig()
        assert cfg.normalization_method == "pos_hk"

    def test_default_de_test(self):
        cfg = NcountrConfig()
        assert cfg.de_test == "mannwhitneyu"

    def test_default_de_correction(self):
        cfg = NcountrConfig()
        assert cfg.de_correction == "fdr_bh"

    def test_default_group_column(self):
        cfg = NcountrConfig()
        assert cfg.group_column == "group"

    def test_default_comparison_empty(self):
        cfg = NcountrConfig()
        assert cfg.comparison == []

    def test_default_exclude_samples_empty(self):
        cfg = NcountrConfig()
        assert cfg.exclude_samples == []

    def test_default_negative_control_empty(self):
        cfg = NcountrConfig()
        assert cfg.negative_control_samples == []

    def test_default_gene_sets_empty(self):
        cfg = NcountrConfig()
        assert cfg.gene_sets == {}

    def test_default_cross_platform_disabled(self):
        cfg = NcountrConfig()
        assert cfg.cross_platform_enabled is False
        assert cfg.cross_platform == {}
