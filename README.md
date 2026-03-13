# ncountr

A Python pipeline for Nanostring nCounter data analysis.

## Features

- **RCC file parsing** — Read `.RCC` files from nCounter runs into tidy DataFrames
- **Quality control** — FOV ratio, positive control linearity, negative background, binding density checks
- **Normalization** — Positive control, housekeeping, and background-subtraction methods
- **Differential expression** — Mann-Whitney U or t-test with FDR correction
- **Gene set scoring** — Built-in IFN/JAK-STAT and cell marker gene sets, or bring your own
- **Cross-platform validation** (optional) — Compare with external expression data (h5ad, CSV)
- **Publication-ready plots** — QC panels, volcano plots, heatmaps, correlation scatter plots

## Installation

```bash
pip install ncountr

# With cross-platform validation support (requires scanpy/anndata):
pip install ncountr[crossplatform]

# Development:
pip install ncountr[dev]
```

## Quick start

```bash
# Generate a template config
ncountr init > my_config.yaml

# Edit the config, then run the full pipeline
ncountr run my_config.yaml

# Or run individual steps
ncountr parse --rcc-dir /path/to/RCC/
ncountr qc --counts raw_counts.csv
ncountr normalize --counts raw_counts.csv
ncountr de --counts normalized.csv --groups treated:S1,S2 control:S3,S4
```

## Python API

```python
import ncountr

# Parse RCC files
experiment = ncountr.read_rcc("/path/to/RCC/", sample_id_pattern=r"(\d+)")

# QC
qc_results = ncountr.qc(experiment)

# Normalize
ncountr.normalize(experiment, method="pos_hk")

# Differential expression
de_results = ncountr.de(
    experiment,
    group_a=["S1", "S2"],
    group_b=["S3", "S4"],
)

# Gene set scoring
scores = ncountr.score_gene_set(experiment, gene_set="IFN_JAKSTAT")
```

## Config format

ncountr uses YAML configuration files. See `examples/basic_config.yaml` for a full template.

## License

MIT
