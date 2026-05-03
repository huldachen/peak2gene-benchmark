# GutCRE Benchmark

[![Python](https://img.shields.io/badge/python-3.11-blue)](https://www.python.org)
[![scanpy](https://img.shields.io/badge/scanpy-≥1.10-1c6b3c)](https://scanpy.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

If you're picking between Cicero, the Activity-by-Contact model, ArchR peak2gene, or ENCODE-rE2G to link regulatory peaks to genes from single-cell data, the existing head-to-head benchmarks (Gschwind et al. 2025; Nasser et al. 2021) compare some of these methods on CRISPR-perturbation ground truth in cell lines - predominantly K562, with smaller HCT116 and eQTL panels. **This project complements those benchmarks with a different lens** - a six-method comparison (including Cicero canonical, ArchR's faster Pearson shortcut, paired-Multiome correlation, and cross-cohort anchor-transfer alongside ABC and rE2G) on a primary-tissue single-cell atlas (Hickey et al. 2023), ranked by inflammatory bowel disease (IBD) GWAS heritability partitioning (Finucane et al. 2015; de Lange et al. 2017) rather than CRISPRi gene-by-gene ground truth. The result is three concrete methodological claims that the original publications do not address: ArchR's raw-Pearson shortcut is not statistically equivalent to canonical Cicero for tau\*-based analyses; cross-cohort anchor-transfer integration silently fails at the cell-type level for most cell types; and the K562-calibrated ENCODE-rE2G threshold collapses onto a distance baseline when applied to gut tissue.

---

## What's at stake

A team picking a peak-gene linkage method today reads four or five method papers, each validated against its preferred CRISPRi or eQTL panel in its preferred cell type, and is left to guess which one transfers to their disease and tissue. Choose wrong, and the downstream fine-mapping or target-prioritisation analysis inherits a silent bias from upstream linkage. The community has converged on stratified LD-score regression (Finucane et al. 2015) with the baseline-LD v2.2 model (Gazal et al. 2017) as the standard tool for asking *"how much of a trait's heritability concentrates in this annotation?"* - a metric that doesn't care which method produced the annotation. That makes it the right lens for a comparison the field has been postponing.

## Three claims this benchmark tests

Each of the following is a claim that recurs in the peak-gene linkage literature, and each has a specific empirical answer the original publications do not provide.

### 1. ArchR's raw-Pearson shortcut is "Cicero-equivalent" - *not for tau\* analyses.*

Pliner et al. (2018) define Cicero as graphical lasso applied to overlapping 500 kb windows of the genome with a distance-based penalty: the regularisation `rho_{ij}` decreases for peaks closer in genomic distance, so nearby pairs can retain strong edges while distant pairs need stronger evidence. SCENIC+ (Bravo González-Blas et al. 2023) uses canonical Cicero internally and inherits this behaviour. ArchR's `getCoAccessibility()` (Granja et al. 2021) takes a faster path - Pearson correlation between meta-cell aggregated accessibility profiles, computed via the C++ `rowCorCpp` helper - and is widely used as a Cicero substitute when the canonical R package is impractical at scale. **The two are not statistically equivalent.** Running both as separate methods on this benchmark, the GLASSO step shifts per-SNP IBD enrichment by up to 4.3x and flips the tau\* coefficient (the unique heritability contribution beyond baseline-LD) from -3.4 in one cell - meaning the raw-Pearson annotation is *redundant with baseline-LD* - to +0.9 - meaning the GLASSO annotation is benign. Headline rank-order against the other linkage families is preserved on this dataset, but anyone running tau\*-based or baseline-LD-conditioned analyses with the ArchR shortcut should know they're trading something measurable away.

### 2. Anchor-transfer integration is reliable at the cell-type level - *only for cells with distinctive transcriptomes.*

Stuart et al. (2019) introduced anchor-transfer integration and validated it at the cluster level; Seurat, Signac, ArchR, and SCENIC+ all build on this foundation, and downstream uses rarely audit anchor quality at the per-cell, per-cell-type level. **At that level the integration silently fails for most cell types.** Mean per-cell anchor agreement on the Hickey et al. (2023) colon atlas is 0.6-9.7% across compartments. Out of 44 analysed cell types, only Enterocytes (61%), Myofibroblasts/SM-DES-high (68%), and Glia (29%) clear 25%; the other 41 cell types have under 10% correctly typed anchors. The aggregate compartment-level metrics that anchor-transfer pipelines report are not wrong - but they hide a silent failure mode for most cell types within the compartment. The diagnostic that catches this is <30 lines of code and runs in under a minute (Figure 4); we recommend it as a standard sanity check before any cell-type-specific claim from cross-cohort integration.

### 3. The ENCODE-rE2G K562-calibrated threshold transfers to other tissues - *not to gut.*

Gschwind et al. (2025) calibrate the rE2G score threshold (0.179) at 70% recall on the K562 CRISPRi screen and benchmark rE2G against ABC, distance, and several other features on >10,000 element-gene pairs from CRISPR perturbation. **Applied off-the-shelf to colon, the K562 threshold is over-permissive by enough to make rE2G numerically indistinguishable from a distance baseline.** ~46% of candidate pairs in the immune compartment pass the threshold; the resulting binary annotation overlaps Method 1 (distance window) to within 0.01x per-SNP enrichment in all nine trait x compartment cells. The continuous-score variant doesn't fix it - the score distribution among non-zero SNPs is saturated near 1.0, so binary and continuous formulations yield numerically equivalent heritability partitions. Until gut-matched CRISPRi recalibration exists, rank-based downstream use is more defensible than the K562 binary call.

Together these findings move the peak-gene linkage literature from a series of *"my method works in my hands"* reports into a ranking with three concrete, reproducible caveats - and each caveat costs less than an hour to incorporate into an existing pipeline.

## Six methods, one benchmark

| # | Method | Core idea | Input modalities |
|---|--------|-----------|-----------------|
| 1 | **Distance window** | Every peak within +/-500 kb of a gene TSS | ATAC |
| 2a | **Co-accessibility (Pearson)** | Peak-peak Pearson r in KNN meta-cells (ArchR-style) | ATAC |
| 2b | **Co-accessibility (GLASSO)** | Partial correlation with distance penalty (canonical Cicero) | ATAC |
| 3 | **Activity-by-Contact (ABC)** | Activity x distance-decayed contact | ATAC |
| 4A | **Paired-Multiome correlation** | Direct ATAC x RNA correlation (same cells) | paired ATAC+RNA |
| 4B | **Cross-cohort anchor transfer** | Integrate RNA from a reference cohort -> correlation | ATAC + independent RNA |
| 5 | **ENCODE-rE2G** | Pretrained supervised model (CRISPR-validated) | ATAC + rE2G pickle |

All six methods are implemented as self-contained Python scripts (`scripts/linkage/method1_distance.py` through `method5_re2g.py`, plus `method2b_coaccess_glasso.py`), with a shared unified output schema so their outputs are directly comparable. Each script reproduces the canonical algorithm from its primary citation - details in [`docs/methods.md`](docs/methods.md).

Methods 2a and 2b together quantify a methodological question that recurs in the multi-omics literature but is rarely answered empirically: *does the canonical Cicero GLASSO step change the heritability signal versus the raw-Pearson approximation used by ArchR?* On this benchmark the answer is "yes, in two specific and reproducible ways" - see Section 1 of [`results/DISCUSSION.md`](results/DISCUSSION.md).

## Headline findings

**1. Paired Multiome gives the cleanest signal per SNP.**
Method 4A concentrates IBD heritability 10-21x in just 0.2-0.5 % of SNPs (p ≤ 0.05 across all three gut compartments).

**2. Cross-cohort anchor-transfer integration can silently fail at cell-type resolution.**
Compartment-aggregate metrics look reasonable, but a per-cell anchor-agreement diagnostic reveals that most cell types have < 10 % correctly typed anchors - a field-level caveat for any pipeline using standard gene-activity-bridge integration.

**3. Pretrained supervised models (rE2G) don't transfer off-the-shelf.**
The K562-calibrated threshold makes rE2G indistinguishable from the distance baseline on gut data, at either binary or continuous scoring.

**4. Negative controls validate the pipeline.**
Height GWAS enrichment correctly concentrates in the stromal compartment (p ≤ 10^-9) - the biology predicted by bone/connective-tissue regulation.

See [`results/DISCUSSION.md`](results/DISCUSSION.md) for full interpretation and [`results/figures/fig3_ldsc_benchmark.png`](results/figures/fig3_ldsc_benchmark.png) for the main benchmark figure.

## Figures

Four publication-quality figures in `results/figures/` (PNG). See [`results/FIGURES_AND_TABLES.md`](results/FIGURES_AND_TABLES.md) for the full index and [`results/DISCUSSION.md`](results/DISCUSSION.md) for interpretation.

### Figure 1 - Dataset overview and QC
<p align="center">
  <img src="results/figures/fig1_data_overview.png" alt="Dataset overview and QC" width="880"/>
</p>

*Two-cohort dataset structure (RNA-only + Multiome), per-compartment cell composition, RNA UMAP coloured by compartment, and per-modality QC metrics.*

### Figure 2 - Data flow across the six methods
<p align="center">
  <img src="results/figures/fig2_data_flow.png" alt="Data flow schematic" width="880"/>
</p>

*How RNA and ATAC from the two cohorts feed each of the six peak-gene linkage methods, with a ribbon showing the concrete question each method answers.*

### Figure 3 - sLDSC benchmark (main result)
<p align="center">
  <img src="results/figures/fig3_ldsc_benchmark.png" alt="sLDSC benchmark forest plot" width="880"/>
</p>

*Two-row forest plot. Row A: fold-enrichment +/- SE across 7 method keys (1, 2a, 2b, 3, 4A, 4B, 5) x 3 compartments x 3 traits. Row B: tau\* (coefficient z-score) - the unique heritability contribution beyond the baseline-LD v2.2 model, size-independent. Dot size encodes annotation coverage (prop_snps). IBD is primary; Height and Educational Attainment are non-gut negative controls.*

### Figure 4 - Cross-cohort anchor-agreement diagnostic
<p align="center">
  <img src="results/figures/fig4_anchor_agreement.png" alt="Anchor-agreement diagnostic" width="720"/>
</p>

*Per-cell-type anchor-agreement in cross-cohort integration: most cell types have near-zero agreement (red bars), revealing a silent failure mode that aggregate metrics miss.*

## Repository structure

```
github.ver/
├── README.md                    # this file
├── environment.yaml             # conda env spec (Python + libraries) - used once at install
├── config/config.yaml           # analysis parameters + reference paths - used every run
├── docs/
│   ├── project_overview.md      # audience: HR / hiring manager
│   ├── methods.md               # educational walkthrough of each method
│   └── results_summary.md       # plain-English findings
├── scripts/
│   ├── data/                    # Dryad downloader + Hickey loader
│   ├── linkage/                 # methods 1, 2a, 2b, 3, 4A, 4B, 5 + cross-cohort diagnostic
│   ├── ldsc/                    # sLDSC annotation -> LD scores -> h^2 pipeline
│   └── figures/                 # reproducible figure + table generators
├── notebooks/
│   └── 01_hickey_qc_executed.ipynb      # QC + clustering + cell-type annotation
└── results/
    ├── DISCUSSION.md                 # narrative writeup (~1,260 words, 16 citations)
    ├── FIGURES_AND_TABLES.md         # index + quick-reference
    ├── figures/fig1-4.png            # 4 figures (PNG)
    ├── tables/table1-4.tsv           # 4 tables (TSV, machine-readable)
    └── linkage/ + ldsc/              # raw per-method outputs + sLDSC results
```

**A note on `environment.yaml` vs `config/config.yaml`:** `environment.yaml` is the tooling spec (Python version + library versions, used only by `conda env create`). `config/config.yaml` is the analysis parameter file (thresholds, windows, reference paths, loaded by scripts at runtime). Keeping them separate is standard practice - they have different lifecycles and different consumers, and merging them would conflate environment reproducibility with analysis tuning.

## Quick-start

```bash
conda env create -f environment.yaml
conda activate omics

# Download reference files (one-time, ~3 GB)
# see scripts/data/download_hickey_dryad.py for details

# Run the six linkage methods (each takes --all-compartments)
python scripts/linkage/method1_distance.py            --all-compartments --gtf ...
python scripts/linkage/method2_cicero.py              --all-compartments --gtf ...   # 2a: raw Pearson
python scripts/linkage/method2b_coaccess_glasso.py    --all-compartments --gtf ...   # 2b: GLASSO partial
python scripts/linkage/method3_abc.py                 --all-compartments --gtf ...
python scripts/linkage/method4_paired.py              --all-compartments --gtf ... \
    --rna-h5ad data/processed/hickey2023/multiome_rna_merged.h5ad
python scripts/linkage/method4_crosscohort.py         --all-compartments --gtf ... \
    --rna-h5ad data/processed/hickey2023/rna_merged.h5ad
python scripts/linkage/method5_re2g.py                --all-compartments --gtf ...

# Build sLDSC annotations, compute LD scores, run partitioned heritability
python scripts/ldsc/build_annotations.py --all
python scripts/ldsc/compute_ldscores.py --all --parallel 4
python scripts/ldsc/run_h2.py --all --parallel 6
python scripts/ldsc/aggregate_results.py

# Regenerate figures and tables from scratch
python scripts/figures/figure1_data_overview.py
python scripts/figures/figure3_ldsc_benchmark.py   # main benchmark figure
python scripts/figures/make_tables.py
```

## Runtime reference

Wall-clock times from a development run on Apple M-series (macOS, `omics` conda env, 10-core parallel where applicable). Data-loading / sLDSC steps dominate.

### Per-script runtimes

| Stage | Script | Typical runtime | Notes |
|-------|--------|----------------:|-------|
| Data | `scripts/data/download_hickey_dryad.py` | manual (user clicks) | ~10 min of clicks + downloads |
| Data | `scripts/data/load_hickey.py --mode rna` | ~2 min | 4 RNA-only samples |
| Data | `scripts/data/load_hickey.py --mode rna-multiome` | ~10-15 min | 42 Multiome samples; large MTX files |
| Data | `scripts/data/load_hickey.py --mode atac --source peaks` | ~3-5 min | 3 compartment peak matrices |
| Linkage | `scripts/linkage/method1_distance.py --all-compartments` | ~3.5 min | bedtools slop + intersect |
| Linkage | `scripts/linkage/method2_cicero.py --all-compartments` | ~15-20 min | Method 2a: LSI + KNN meta-cells + per-chrom Pearson |
| Linkage | `scripts/linkage/method2b_coaccess_glasso.py --all-compartments` | ~1.5-2 h | Method 2b: GLASSO with distance penalty; chromosome-parallel |
| Linkage | `scripts/linkage/method3_abc.py --all-compartments` | ~1-2 min | vectorised per-gene normalisation |
| Linkage | `scripts/linkage/method4_paired.py --all-compartments` | ~10-15 min | paired barcode matching + meta-cell correlation |
| Linkage | `scripts/linkage/method4_crosscohort.py --all-compartments` | ~18-25 min | anchor KNN transfer is the bottleneck |
| Linkage | `scripts/linkage/method4b_anchor_diagnostic.py --all` | ~3 min | per-compartment anchor agreement |
| Linkage | `scripts/linkage/method5_re2g.py --all-compartments` | ~1-2 min | features -> pretrained LR predict |
| sLDSC | `scripts/ldsc/build_annotations.py --all` | ~70 min | 21 annotations x 22 chroms (liftover + make_annot) |
| sLDSC | `scripts/ldsc/build_method5_continuous.py --all` | ~1 min | 3 continuous annotations (output dropped from final reports as numerically equivalent to binary) |
| sLDSC | `scripts/ldsc/compute_ldscores.py --all --parallel 4` | **~6 hours** | 21 annotations x 22 chroms; dominant step |
| sLDSC | `scripts/ldsc/run_h2.py --all --parallel 6` | ~30-40 min | 72 (trait, annotation) combinations; idempotent (skips existing .results) |
| sLDSC | `scripts/ldsc/aggregate_results.py` | ~1 s | parses .results/.log files |
| Figures | `scripts/figures/figure{1,2,3,4}_*.py` | 5-30 s each | matplotlib rendering |
| Figures | `scripts/figures/make_tables.py` | ~1 s | pandas -> TSV (Tables 1-4) |

### End-to-end total
Full pipeline from raw Dryad downloads to final figures: **~9-10 hours wall-clock**, of which ~6 hours is LD score computation.

### First-time setup (one-off, not re-run)
- Reference data downloads (1000G, baseline-LD, HapMap3, GWAS sumstats): ~10 min of clicking + ~1 GB download
- Conda env build (`environment.yaml`): ~5 min
- LDSC env build (`bioconda ldsc` via Rosetta): ~5 min
- ENCODE-rE2G clone: ~1 min

---

## What this project demonstrates

This is a methods-benchmark project involving **multi-modal single-cell data integration**, **statistical genetics (sLDSC partitioned heritability)**, **reproducible pipeline design**, and **careful interpretation of method trade-offs**. Technical skills shown include:

- Single-cell RNA-seq + ATAC-seq QC, clustering, and cell-type annotation (scanpy, snapATAC2)
- Multi-modal integration via anchor transfer, KNN in shared latent spaces, and paired-barcode correlation
- Implementation of published enhancer-prediction methods from primary literature (Cicero, ABC, ArchR peak2gene, rE2G)
- Stratified LD-score regression pipeline (ldsc, Python 2.7 via Rosetta)
- Reproducible, scriptable figures and tables (matplotlib, pandas)
- Method comparison framed around a concrete biological question (disease heritability partitioning)

See [`docs/project_overview.md`](docs/project_overview.md) for an expanded educational description.

## Citations

- Hickey J. W. et al. *Nature* **619**, 572 (2023) - primary dataset
- de Lange K. M. et al. *Nat. Genet.* **49**, 256 (2017) - IBD GWAS
- Finucane H. K. et al. *Nat. Genet.* **47**, 1228 (2015) - sLDSC
- Fulco C. P. et al. *Nat. Genet.* **51**, 1664 (2019); Nasser J. et al. *Nature* **593**, 238 (2021) - ABC
- Pliner H. A. et al. *Mol. Cell* **71**, 858 (2018) - Cicero
- Granja J. M. et al. *Nat. Genet.* **53**, 403 (2021) - ArchR
- Stuart T. et al. *Cell* **177**, 1888 (2019) - Seurat/Signac anchor transfer
- Gschwind A. R. et al. *Nature* (2025) - ENCODE-rE2G

Full reference list in [`results/DISCUSSION.md`](results/DISCUSSION.md).

