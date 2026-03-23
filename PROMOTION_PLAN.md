# ncountr Promotion Plan

## Current Position

- ~1,135 nCounter instruments worldwide, ~2,000–5,000 active analysts/year, ~100–130 publications/year
- 11 R packages exist (NanoStringNorm, NanoStringDiff, NanoTube, NACHO, nanostringr, etc.)
- Zero functional Python packages — ncountr is the only working option
- Platform stable under Bruker (acquired NanoString for $393M in 2024, actively developing new panels)
- A 2024 review paper (Chilimoniuk et al., CSBJ) explicitly identified the lack of a standardized pipeline as a barrier to wider nCounter adoption

## Track 1: Self-Validation (Build Credibility First)

Download published nCounter datasets from GEO, reproduce their results with ncountr, and turn each into a vignette. This gives three things at once: real-world testing to find bugs/gaps, reproducible case studies for the paper, and tutorial content that attracts users.

### Priority Datasets (all confirmed RCC files on GEO)

| GSE | Topic | Panel | N | Why this one |
|-----|-------|-------|---|-------------|
| GSE275334 | Long COVID / ME/CFS | Immune Exhaustion (785 genes) | 47 | Directly relevant to LCM manuscript; 3-group design (healthy vs long COVID vs ME/CFS) |
| GSE140901 | HCC immunotherapy | PanCancer Immune (770 genes) | 24 | Classic panel; responder vs non-responder to checkpoint inhibitors |
| GSE117751 | Autoimmune retinopathy | Human Immunology | 42 | 3-group design; already used in NanoTube tutorials (comparability) |
| GSE268012 | IFN in macrophages | Human Metabolism (768 genes) | 24 | Clean 4-arm factorial design; IFN biology matches our built-in gene set |
| GSE74821 | Breast cancer PAM50 | PAM50 Custom | 1,321 | Stress-test at scale; used in NACHO demos |

### Additional Datasets Worth Exploring

| GSE | Topic | Panel | N | Notes |
|-----|-------|-------|---|-------|
| GSE165745 | Melanoma anti-PD-1 | Wnt Pathways | 24 | Responder vs non-responder; ~67 citations |
| GSE229605 | Psoriasis T cell stimulation | Custom | 92 | Multi-factorial paired design |
| GSE135535 | Prostate cancer miRNA | miRNA | 320 | Two independent cohorts; tests miRNA support |
| GSE70970 | Nasopharyngeal carcinoma | miRNA v1.0 | 263 | Large; training + validation cohort |
| GSE130557 | Colorectal cancer | PanCancer Pathways (770 genes) | 2 | Tiny but quick smoke test |

### What to Do with Each Dataset

1. Download RCC files from GEO (`GSE*_RAW.tar`)
2. Run full ncountr pipeline (parse → QC → normalize → DE → gene set scoring)
3. Compare results to published figures/tables — can we reproduce their key findings?
4. Document what worked, what broke, what's missing
5. Save each as a Jupyter notebook in `examples/vignettes/`

### Expected Improvements to Discover

- Edge cases in RCC parsing (different instrument versions, miRNA panels)
- Missing normalization options (geNorm reference selection, RUV)
- Panel-specific features (miRNA ligation normalization — only nSolver handles this currently)
- Batch correction needs (GSE74821 has 1,321 samples across many cartridges)
- Better default QC thresholds calibrated against real failure modes

Each fix makes ncountr more robust. Each vignette becomes promotional content.

## Track 2: Promotion

### Phase 1: Paper (Month 1–2)

**Option A — JOSS (faster, lower barrier):**
- 2-page summary: what ncountr does, why Python, one benchmark figure
- Review takes 2–4 weeks, published open-access
- Requirements already met: working tests, docs, example data
- Gets a citable DOI quickly

**Option B — Bioinformatics application note (more impactful):**
- 2 pages + supplementary benchmarks
- Show ncountr reproducing results from 3–4 published datasets
- Compare normalization accuracy against NanoStringNorm/NanoTube on the same data
- Review takes 2–3 months

**Recommended**: JOSS first (get a citable DOI fast), then Bioinformatics benchmark paper later with full validation.

### Phase 2: Ecosystem Integration (Month 1, parallel with paper)

1. **Add `to_anndata()` export** — wrap count matrix + metadata into AnnData object. Trivial code, but unlocks scverse ecosystem listing
2. **Submit to scverse community packages** — gets listed on scverse.org/packages, visible to every scanpy user
3. **Conda-forge recipe** — `conda install ncountr` reaches users who don't use pip
4. **ReadTheDocs site** — auto-generated from docstrings + vignette notebooks
5. **Binder/Colab notebook** — one-click demo using example RCC files, zero install needed

### Phase 3: Community Seeding (Ongoing)

- Answer NanoString questions on Biostars / StackOverflow / GitHub Discussions with ncountr solutions (~20–30 existing NanoString questions on Biostars)
- Post vignettes on Bluesky/Twitter: "Reproduced [paper X] results in 15 lines of Python" with a figure. Hashtags: #Bioinformatics, #NanoString, #Python, #GeneExpression
- Submit to awesome-bioinformatics and awesome-single-cell lists
- Email authors of the 2024 review paper (Chilimoniuk et al.) — they catalogued 11 R packages and zero Python options; ncountr is exactly what they noted was missing
- Tag Bruker Spatial Biology on social media posts

### Phase 4: Comparative Benchmark (Month 3–4)

Run the same datasets through ncountr, NanoTube, NACHO, and nanostringr. Compare:

- Normalization accuracy (correlation with RNA-seq ground truth if available)
- DE sensitivity/specificity
- Runtime performance
- Lines of code needed for equivalent analysis

This becomes the Bioinformatics paper and positions ncountr as a serious alternative, not just "another tool."

### The Pitch (use everywhere)

> 11 R packages exist for nCounter data, but until now there was no Python option. ncountr fills that gap — parse RCC files, run QC, normalize, and test for differential expression in a single config file or Python API. Includes cross-platform validation against scRNA-seq.

## Execution Timeline

```
Week 1–2:  Download GSE275334 + GSE140901 + GSE117751
           Run through ncountr, fix bugs, create vignette notebooks

Week 3–4:  Add to_anndata() export, set up ReadTheDocs, conda-forge recipe
           Download GSE268012 + GSE74821, more vignettes

Month 2:   Submit JOSS paper (vignettes as supplementary)
           Submit to scverse community packages
           Post first vignette on social media

Month 3:   Start comparative benchmark against R packages
           Post more vignettes, answer community questions

Month 4:   Submit Bioinformatics application note with benchmarks
```

## Realistic Expectations

- nf-core/nanostring (Nextflow pipeline): 18 stars
- PyDESeq2 (closest analogue — Python port of R-only workflow): 742 stars after Bioinformatics paper + scverse listing
- Target: 50–200 stars is achievable for a well-positioned niche tool
- Each vignette reproducing a published dataset is both a quality check and a promotion opportunity

## Reference: Comparable Packages and Their Strategies

| Package | Stars | Strategy |
|---------|-------|----------|
| PyDESeq2 | 742 | scverse listing + Bioinformatics paper + anndata integration |
| CellBender | 384 | Broad Institute brand + fills real pain point |
| decoupler | 263 | scverse ecosystem + used by downstream tools |
| diffxpy | 206 | scverse adjacent |
| nf-core/nanostring | 18 | Nextflow ecosystem + Bioinformatics paper |

## Reference: Key Review Paper

Chilimoniuk J, Erol A, Rodiger S, Burdukiewicz M. "Challenges and opportunities in processing NanoString nCounter data." *Computational and Structural Biotechnology Journal*, 2024. PMID: 38736697.

- Evaluates 11 R packages across 5-step workflow
- Recommends NanoTube as most robust for mRNA data
- Notes zero Python options exist
- Benchmark framework could be adapted for ncountr comparison

## Reference: nCounter vs RNA-seq Concordance Studies

| Study | Design | Key Result |
|-------|--------|------------|
| Ebola NHP (BMC Genomics 2025) | 62 NHP samples, 769 genes | Spearman rho 0.78–0.88 |
| Respiratory virus organoids (Front Genet 2024) | 754 overlapping genes | Both highlight ISG15, MX1, RSAD2, OAS |
| Heart failure PBMC (J Biol Methods 2020) | 770 genes, 25 patients | Intra-sample r = 0.904 |
