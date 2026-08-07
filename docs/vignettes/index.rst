Vignettes
=========

Five worked examples that run the ncountr pipeline end to end on public NanoString
nCounter datasets from GEO, spanning 1,458 samples and five different assay panels.

Each notebook downloads its own RCC files with :func:`ncountr.io.geo.fetch_geo`, so they
can be run without any manual data preparation. They are stored without execution
outputs; run them locally or open them in Colab or Binder to see results.

The notebooks live under ``examples/vignettes/`` in the repository and are linked
individually below.

.. list-table::
   :header-rows: 1
   :widths: 16 26 10 24 24

   * - Dataset
     - Biological context
     - Samples
     - Panel
     - Shows
   * - `GSE275334 <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE275334_long_covid.ipynb>`_
     - Long COVID and ME/CFS vs healthy controls (PBMC)
     - 47
     - Human Immune Exhaustion v1.0 (773 genes)
     - Full pipeline, three-group design, **AnnData → scanpy handoff**
   * - `GSE140901 <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE140901_hcc_immunotherapy.ipynb>`_
     - Checkpoint-inhibitor response in hepatocellular carcinoma
     - 24
     - PanCancer Immune Profiling v1.1 (730 genes)
     - Two-group responder vs non-responder comparison
   * - `GSE117751 <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE117751_autoimmune_retinopathy.ipynb>`_
     - Autoimmune retinopathy and retinitis pigmentosa vs controls
     - 42
     - Human Immunology v2 (579 genes)
     - Three-group design; dataset also used by the R NanoTube tutorials
   * - `GSE268012 <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE268012_ifn_macrophages.ipynb>`_
     - Interferon-treated monocyte-derived macrophages
     - 24
     - Human Metabolism v1.0 (748 genes)
     - Four-arm factorial design with paired donors
   * - `GSE74821 <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE74821_breast_cancer_pam50.ipynb>`_
     - FFPE breast tumours from the CALGB 9741 trial
     - 1,321
     - RUO PAM50 (50 genes)
     - Scale test across 116 cartridges

.. _vignette-long-covid:

GSE275334 — Long COVID / ME/CFS
-------------------------------

`Open the notebook <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE275334_long_covid.ipynb>`_

47 PBMC samples on the Human Immune Exhaustion panel (773 endogenous + 12 housekeeping
genes): 18 healthy controls, 15 Long COVID, 14 ME/CFS. This is the most complete of the
five vignettes and the recommended starting point. It covers GEO download, RCC parsing,
QC (including one sample excluded for poor positive-control linearity), positive-control
plus housekeeping normalization, two differential expression contrasts, IFN/JAK-STAT gene
set scoring, and a heatmap of top genes.

**This is also the scverse interoperability example.** Section 9 carries the exported
:func:`ncountr.to_anndata` object into scanpy for PCA, principal components coloured by
the QC and lane metadata that travel across in ``.obs``, a group dendrogram,
``rank_genes_groups``, and a matrixplot — then cross-checks that re-running the same
contrast in scanpy reproduces the ``ncountr.de()`` gene ranking.

GSE140901 — HCC immunotherapy response
--------------------------------------

`Open the notebook <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE140901_hcc_immunotherapy.ipynb>`_

24 hepatocellular carcinoma tumours profiled on the PanCancer Immune Profiling panel
(730 endogenous + 40 housekeeping genes) from patients treated with anti-PD-1/PD-L1
checkpoint inhibitors, split into 13 responders and 11 non-responders. A compact
two-group worked example: parsing, QC, normalization, differential expression, pathway
scoring, and a heatmap of the top 25 genes.

GSE117751 — Autoimmune retinopathy
----------------------------------

`Open the notebook <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE117751_autoimmune_retinopathy.ipynb>`_

42 whole blood (PAXgene) samples on the Human Immunology v2 panel (579 endogenous + 15
housekeeping genes): 14 autoimmune retinopathy, 14 retinitis pigmentosa, 14 healthy
controls. Demonstrates a three-group design with all pairwise contrasts. This dataset is
also used in the tutorials for the R package NanoTube, which makes it convenient for
comparing ncountr's output against an established implementation.

GSE268012 — Interferon-treated macrophages
------------------------------------------

`Open the notebook <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE268012_ifn_macrophages.ipynb>`_

24 human monocyte-derived macrophages on the Human Metabolism panel (748 endogenous + 20
housekeeping genes), in a four-arm factorial design with paired donors: control, IFNα2a,
IFNα2b, and IFNβ. The built-in IFN/JAK-STAT gene set is directly relevant here, so this
vignette is the clearest illustration of gene set scoring against a known perturbation.

GSE74821 — PAM50 breast cancer scale test
-----------------------------------------

`Open the notebook <https://github.com/princello/ncountr/blob/main/examples/vignettes/GSE74821_breast_cancer_pam50.ipynb>`_

1,321 FFPE breast tumours from the CALGB 9741 clinical trial on the RUO PAM50 panel
(50 endogenous + 8 housekeeping genes), across 116 cartridges. This vignette exists to
stress-test the parser and QC at a sample size well beyond a typical nCounter study, and
to check that cartridge-level batch structure is handled sensibly.
