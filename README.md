# ncountr

[![CI](https://github.com/princello/ncountr/actions/workflows/ci.yml/badge.svg)](https://github.com/princello/ncountr/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/ncountr)](https://pypi.org/project/ncountr/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/ncountr)](https://anaconda.org/conda-forge/ncountr)
[![Docs](https://readthedocs.org/projects/ncountr/badge/?version=latest)](https://ncountr.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A Python pipeline for analyzing Nanostring nCounter gene expression data — from raw instrument files to differential expression, pathway scoring, and publication-ready figures.

**11 R packages exist for nCounter data, but until now there was no Python option.** ncountr fills that gap.

## The problem ncountr solves

The Nanostring nCounter is a mid-throughput gene expression platform used in translational research, clinical trials, and biomarker studies. Unlike RNA-seq, it directly counts individual mRNA molecules for a pre-selected panel of genes (typically 100–800), without amplification or library preparation. This makes the data simpler to work with, but it still requires platform-specific normalization and quality control before you can trust the results.

Each nCounter run produces one `.RCC` file per sample — a proprietary text format containing raw molecule counts for your target genes, along with built-in controls:

- **Positive controls** (synthetic RNA spikes at known concentrations) to verify assay linearity
- **Negative controls** (no-template probes) to measure background noise
- **Housekeeping genes** (stably expressed reference genes) to correct for RNA input differences

Most researchers handle these steps with NanoString's nSolver software (Windows GUI) or the NanoStringNorm R package. ncountr provides a Python alternative that runs from the command line or as a library, produces all standard QC and normalization steps, and adds differential expression, gene set scoring, and optional cross-platform validation in a single reproducible config file.

## Installation

```bash
pip install ncountr

# With cross-platform validation support (adds scanpy/anndata):
pip install ncountr[crossplatform]
```

Requires Python 3.9+.

## What's new

- **`to_anndata()` export** — convert nCounter data to AnnData for scanpy/scverse integration
- **`ncountr fetch-geo GSE275334`** — download RCC files directly from NCBI GEO
- **`sample_id_from="filename"`** — extract sample IDs from filenames when internal RCC IDs are inconsistent
- **[Documentation](https://ncountr.readthedocs.io)** — Sphinx docs with API reference and CLI docs
- **5 validated vignettes** — reproduced published results from GEO datasets (see below)

## Validated on real-world data

ncountr has been tested against 5 published nCounter datasets totaling 1,458 samples:

| Dataset | Panel | Samples | Key result |
|---------|-------|---------|------------|
| [GSE275334](examples/vignettes/GSE275334_long_covid.ipynb) | Immune Exhaustion (773 genes) | 47 | Long COVID / ME/CFS 3-group design |
| [GSE140901](examples/vignettes/GSE140901_hcc_immunotherapy.ipynb) | PanCancer Immune (730 genes) | 24 | ICI responder vs non-responder |
| [GSE117751](examples/vignettes/GSE117751_autoimmune_retinopathy.ipynb) | Human Immunology (579 genes) | 42 | AIR vs RP vs Control (also in NanoTube tutorials) |
| [GSE268012](examples/vignettes/GSE268012_ifn_macrophages.ipynb) | Human Metabolism (748 genes) | 24 | IFN factorial: 232 DE genes for IFNβ |
| [GSE74821](examples/vignettes/GSE74821_breast_cancer_pam50.ipynb) | PAM50 Custom (50 genes) | 1,321 | **Stress test**: 1,321 samples parsed in 3s |

## Typical workflow

A standard nCounter analysis follows four steps. ncountr handles all of them, either through a single config-driven pipeline or as individual commands.

### 1. Parse the raw data

nCounter instruments write one `.RCC` file per sample. ncountr reads a directory of these files and organizes them into a structured experiment object:

```bash
ncountr parse --rcc-dir /path/to/RCC/ --id-pattern '(\d+)'
```

The `--id-pattern` regex extracts sample IDs from filenames (e.g., `Sample_668_Lung.RCC` becomes sample `668`).

### 2. Quality control

Before trusting any results, you need to verify that the assay worked:

```bash
ncountr qc --counts raw_counts.csv
```

This checks four things per sample and produces a 4-panel summary figure:
- **FOV ratio** — Did the scanner image enough fields of view? (threshold: >75%)
- **Positive control linearity** — Do the spike-in controls track their expected concentrations? (R² > 0.95)
- **Negative background** — How much non-specific signal is present?
- **Housekeeping stability** — Are the reference genes consistent across samples?

Samples failing these checks may need to be excluded or flagged.

### 3. Normalize

Raw nCounter counts vary between samples due to differences in RNA input, hybridization efficiency, and imaging. Normalization corrects for these technical factors in two stages:

1. **Positive control normalization** — scales each sample so the synthetic spike-ins match their expected ratios, correcting for assay-level variation
2. **Housekeeping normalization** — further adjusts for RNA input differences using stably expressed reference genes

```bash
ncountr normalize --counts raw_counts.csv --method pos_hk
```

Three methods are available: `pos_only`, `pos_hk` (recommended default), and `pos_hk_bg` (adds background subtraction based on negative controls).

### 4. Differential expression and downstream analysis

With normalized counts, you can compare groups:

```bash
ncountr de --counts normalized.csv --groups treated:S1,S2,S3 control:S4,S5,S6
```

This runs a per-gene statistical test (Mann-Whitney U by default, or t-test), applies FDR correction, and generates a volcano plot. Genes of interest — such as interferon pathway genes — can be highlighted on the plot for quick visual assessment.

## Running the full pipeline

Rather than running each step separately, you can define everything in a single YAML config:

```bash
# Generate a starter config
ncountr init > my_config.yaml

# Edit it with your sample info, then run
ncountr run my_config.yaml
```

A minimal config looks like this:

```yaml
input:
  rcc_dirs:
    - /path/to/RCC/files
  file_pattern: '*.RCC'
  sample_id_pattern: '(\d+)'

output:
  directory: ./results
  figure_format: png
  figure_dpi: 200

samples:
  metadata:
    '668': { group: infected }
    '678': { group: infected }
    '587': { group: control }
    '3601': { group: control }
  group_column: group
  comparison: [infected, control]

normalization:
  method: pos_hk

de:
  test: mannwhitneyu
  correction: fdr_bh

gene_sets:
  IFN_JAKSTAT: builtin
```

This will parse your RCC files, run QC, normalize, test for differential expression between infected and control, score each sample for IFN/JAK-STAT pathway activity, and save all results and figures to `./results/`.

## Python API

For integration into notebooks or custom scripts:

```python
import ncountr

# Load data
experiment = ncountr.read_rcc("/path/to/RCC/", sample_id_pattern=r"(\d+)")

# QC and normalize
ncountr.qc(experiment)
ncountr.normalize(experiment, method="pos_hk")

# Differential expression
de_results = ncountr.de(
    experiment,
    group_a=["668", "678", "680"],
    group_b=["587", "3601", "3603"],
)

# Gene set scoring
scores = ncountr.score_gene_set(experiment, gene_set="IFN_JAKSTAT")

# Access built-in gene sets
ncountr.list_gene_sets()          # ['IFN_JAKSTAT']
ncountr.get_gene_set("IFN_JAKSTAT")  # list of 48 genes
```

### Plotting

```python
from ncountr.plotting import plot_qc, plot_volcano, plot_pathway_scores

# 4-panel QC summary
plot_qc(experiment, output="qc_summary.png")

# Volcano with IFN pathway genes highlighted
plot_volcano(
    de_results,
    highlight_genes=ncountr.get_gene_set("IFN_JAKSTAT"),
    highlight_label="IFN/JAK-STAT genes",
    highlight_color="gold",
    output="volcano.png",
)

# Pathway scores by group
plot_pathway_scores(scores, output="ifn_scores.png")
```

## Download data from GEO

```bash
ncountr fetch-geo GSE275334 -o data/
```

Or from Python:

```python
from ncountr.io.geo import fetch_geo
rcc_dir = fetch_geo("GSE275334", output_dir="data/")
```

## Export to AnnData (scverse integration)

```python
adata = ncountr.to_anndata(experiment)
# adata.X = normalized counts (samples x genes)
# adata.layers["raw"] = raw counts
# adata.obs = sample metadata + QC results + lane info
# adata.var["housekeeping"] = housekeeping gene flag
```

This enables seamless downstream analysis with scanpy, squidpy, decoupler, and other scverse tools.

## Cross-platform validation

When you have expression data from another platform (e.g., single-cell RNA-seq) on the same samples, ncountr can assess how well the two platforms agree. This is useful for confirming that nCounter results are not artifacts of the technology.

Enable it in your config:

```yaml
cross_platform:
  enabled: true
  external_data:
    path: /path/to/adata.h5ad     # scanpy h5ad, or CSV/TSV
    format: h5ad
    pseudobulk_group_by: Sample   # aggregate single cells by sample
  sample_mapping:
    '668': '668'                  # nCounter ID → external ID
    '678': '678'
  negative_control_samples: ['697']  # non-target species controls
```

This adds:
- **Per-sample correlation** — Spearman/Pearson r between platforms for shared genes
- **DE concordance** — what fraction of differentially expressed genes change in the same direction on both platforms
- **Cell composition proxy** — how well nCounter marker gene expression tracks cell type proportions from the other platform
- **Cross-reactivity assessment** — identifies genes with non-specific binding using negative control samples

Install the optional dependencies: `pip install ncountr[crossplatform]`

## Built-in gene sets

| Name | Genes | Description |
|------|-------|-------------|
| `IFN_JAKSTAT` | 48 | Interferon and JAK-STAT signaling pathway (MX1, IFIT1-3, ISG15, STAT1/2, IRF1/7, OAS1-3, CXCL9-11, GBP1-5, JAK1/2, ...) |

Cell type markers are also available for composition estimation:

| Cell type | Markers |
|-----------|---------|
| T cells | CD3D, CD3E, CD3G |
| CD4 T | CD4 |
| CD8 T | CD8A, CD8B |
| Monocytes | CD14, FCGR3A, CD68 |
| B cells | MS4A1, CD19, CD79A |
| NK cells | GNLY, NKG7, KLRD1 |

Custom gene sets can be defined directly in the config or passed as Python lists.

## Configuration reference

Full list of config options:

```yaml
input:
  rcc_dirs: []              # directories containing .RCC files
  file_pattern: '*.RCC'     # glob pattern for RCC files
  sample_id_pattern: '(\d+)'  # regex with one capture group for sample ID

output:
  directory: ./results
  figure_format: png        # png, pdf, svg
  figure_dpi: 200

samples:
  metadata:                 # per-sample key-value pairs
    'S1': { group: treated, batch: A }
  group_column: group       # metadata field used for comparisons
  comparison: [treated, control]  # [group_a, group_b] for DE
  exclude: []               # sample IDs to skip
  negative_control_samples: []  # for cross-reactivity (e.g., NSG controls)

qc:
  fov_ratio_threshold: 0.75
  positive_control_r2_threshold: 0.95

normalization:
  method: pos_hk            # pos_only | pos_hk | pos_hk_bg

de:
  test: mannwhitneyu        # mannwhitneyu | ttest
  correction: fdr_bh        # FDR method from statsmodels

volcano:
  highlight_genes: IFN_JAKSTAT  # builtin name or list of genes
  highlight_label: 'IFN/JAK-STAT genes'
  highlight_color: gold

gene_sets:
  IFN_JAKSTAT: builtin
  my_custom_set: [GENE1, GENE2, GENE3]

cross_platform:             # optional, requires scanpy
  enabled: false
  external_data:
    path: /path/to/data
    format: h5ad            # h5ad | csv | tsv
    pseudobulk_group_by: Sample
  sample_mapping: {}
```

## License

MIT
