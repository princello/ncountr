---
title: 'ncountr: A Python pipeline for Nanostring nCounter gene expression analysis'
tags:
  - Python
  - bioinformatics
  - gene expression
  - nanostring
  - ncounter
  - normalization
  - differential expression
authors:
  - name: Zicheng Wang
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Biomedical Informatics, Columbia University Irving Medical Center, New York, NY, USA
    index: 1
date: 24 March 2026
bibliography: paper.bib
---

# Summary

The Nanostring nCounter platform is a mid-throughput gene expression system
that directly counts mRNA molecules for pre-selected gene panels (typically
50--800 genes) without amplification or library preparation. With roughly
1,100 instruments deployed worldwide and over 100 publications per year,
nCounter remains widely used in translational research, clinical trials,
and biomarker studies [@Chilimoniuk2024]. However, the existing software
ecosystem consists exclusively of Windows-only vendor tools (nSolver) and
11 R packages, with no functional Python alternative. As Python becomes the
dominant language for computational biology --- particularly through the
scanpy/scverse ecosystem --- this gap increasingly forces researchers to
maintain parallel R and Python workflows for what should be a single
analysis.

`ncountr` fills this gap. It is a Python package that provides the complete
nCounter analysis workflow: parsing raw instrument files (.RCC format),
quality control, normalization, differential expression testing, gene set
scoring, and optional cross-platform validation against RNA-seq data.

# Statement of need

A recent systematic review of nCounter data analysis tools [@Chilimoniuk2024]
evaluated 11 R packages across the five-step nCounter workflow and explicitly
identified the absence of a standardized, comprehensive pipeline as a barrier
to wider adoption. That review found no Python packages among the available
tools. `ncountr` addresses this unmet need by providing:

- **End-to-end pipeline** from raw .RCC files to differential expression
  results and publication-ready figures, driven by a single YAML
  configuration file or Python API
- **AnnData export** (`to_anndata()`) that converts nCounter experiments into
  AnnData objects, enabling integration with the scverse ecosystem (scanpy,
  squidpy, decoupler)
- **GEO integration** (`ncountr fetch-geo`) for downloading and extracting
  nCounter datasets directly from NCBI GEO
- **Cross-platform validation** tools for assessing concordance between
  nCounter and RNA-seq measurements on shared samples

# Validation

We validated `ncountr` against five published nCounter datasets spanning
different gene panels, sample sizes, and biological contexts:

| Dataset | Panel | Samples | Key result |
|---------|-------|---------|------------|
| GSE275334 | Immune Exhaustion (773 genes) | 47 | 3-group Long COVID/ME/CFS design; SOCS3, TNFAIP3 top hits |
| GSE140901 | PanCancer Immune (730 genes) | 24 | ICI responders show higher LAG3, HLA-DPA1, PDCD1LG2 |
| GSE117751 | Human Immunology (579 genes) | 42 | 2 FDR-significant genes for AIR vs control (IKBKB, MIF) |
| GSE268012 | Human Metabolism (748 genes) | 24 | IFN$\beta$ drives 232 DE genes vs 0 for IFN$\alpha$ subtypes |
| GSE74821 | PAM50 Custom (50 genes) | 1,321 | Parsed 1,321 files in 3s; 99.5% QC pass rate |

All datasets were parsed, QC-checked, normalized, and analyzed without
modification to the pipeline. One bug was discovered and fixed during
validation (inconsistent sample ID handling across instrument versions),
demonstrating the value of systematic testing against diverse real-world data.
The GSE117751 dataset is also used in NanoTube R package tutorials, enabling
direct cross-tool comparison of results.

# Implementation

`ncountr` is implemented in pure Python with minimal dependencies (NumPy,
pandas, SciPy, statsmodels, matplotlib, PyYAML, Click). It provides both a
command-line interface (`ncountr run config.yaml`) and a library API. The
package follows a modular architecture with separate modules for I/O, QC,
normalization, differential expression, gene set scoring, plotting, and
cross-platform validation.

Normalization implements the standard nCounter protocol: positive control
scaling corrects for assay-level variation, followed by housekeeping gene
normalization for RNA input differences. Differential expression uses
per-gene Mann-Whitney U or t-tests with FDR correction. Gene set scoring
applies z-score normalization followed by per-sample mean aggregation.

# Availability

`ncountr` is available on PyPI (`pip install ncountr`), with documentation
at https://ncountr.readthedocs.io and source code at
https://github.com/princello/ncountr under the MIT license.

# References
