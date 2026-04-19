# Project overview

## In plain English

Every human cell carries the same DNA, but different cell types use different parts of it at different times. The regulatory elements that decide *which* gene gets turned on *where* are called **enhancers**. When disease-associated genetic variants fall in enhancers (as most of them do), understanding which gene each enhancer controls is the missing link between a GWAS hit and actionable biology.

Several computational methods have been proposed to link enhancers (open-chromatin peaks, measured by ATAC-seq) to their target genes (measured by RNA-seq) from single-cell data. They disagree on how to do it:

- some use raw **geometry** — just the genomic distance,
- some measure **correlated openness** between peaks,
- some combine **activity × contact** using approximate 3D-genome models,
- some compute **direct accessibility–expression correlation** when paired data is available,
- and some use **machine-learning models pretrained on CRISPR-validated enhancer–gene pairs**.

Which method should a disease-gene-discovery team actually use? Each publication describes their own method, but a direct head-to-head on a common dataset, evaluated against a concrete biological question, is uncommon. **This project builds that benchmark.**

## What was built

A reproducible benchmark comparing five established peak–gene linkage methods on the Hickey et al. 2023 colon single-cell atlas (~100,000 cells, paired RNA + ATAC multiome + a separate RNA-only cohort). Each method produces peak–gene links; each link set is converted into genome-wide SNP annotations; each annotation set is evaluated via stratified LD-score regression (sLDSC) for its ability to partition inflammatory-bowel-disease (IBD) heritability.

The benchmark includes:

- 5 peak–gene linkage methods (Methods 1–5, plus two variants of Method 4 testing paired vs cross-cohort integration).
- 3 target GWAS traits: IBD (primary), Height (negative control favouring stromal tissue), Educational Attainment (negative control).
- 3 cell-type compartments: epithelial, immune, stromal.
- 63 (method × compartment × trait) sLDSC runs in total.
- A novel cross-cohort integration diagnostic that reveals a previously under-recognised failure mode in standard anchor-transfer pipelines.

## What the results say

Three clear conclusions:

1. **Method 4A (paired Multiome correlation)** concentrates IBD heritability 10–21× in just 0.2–0.5 % of SNPs. When paired data is available, it dominates per-SNP enrichment.

2. **Cross-cohort integration** (Method 4B) looks reasonable at the compartment-aggregate level but fails at cell-type resolution. We introduce an *anchor-agreement diagnostic* (asking "do each ATAC cell's anchor RNA cells share its cell-type label?") that empirically exposes this failure and shows the gene-activity bridge feature captures housekeeping-gene variance rather than cell-type-distinguishing regulation. This is a concrete caveat for any pipeline using the standard Seurat / Signac / ArchR cross-cohort integration without paired data.

3. **The rE2G supervised model** trained on K562 CRISPR data transfers poorly to gut cell types. At its published threshold (calibrated for 70 % recall on K562), it produces an annotation numerically identical to the distance baseline. Better practice: continuous-score downstream use or per-cell-type threshold recalibration.

See [`results/DISCUSSION.md`](../results/DISCUSSION.md) for full interpretation.

## Technical skills demonstrated

- **Single-cell multi-omics**: QC, clustering, cell-type annotation, cross-modality barcode pairing
- **Published-method reimplementation**: each of the 5 methods reproduces its primary citation's algorithm
- **Statistical genetics**: sLDSC pipeline (ldsc, GWAS sumstats munging, baseline-LD v2.2, h² partitioning)
- **Reproducibility**: single conda environment, scripted figures and tables, unified output schema
- **Quality control**: anchor-agreement diagnostic, negative-control trait validation, per-gene FDR correction
- **Large-scale data handling**: 102 k cells × 1.12 M peaks, 45 k genes × 11 k cells, ~10 GB raw + 5 GB LDSC reference data, 2 h of LD score computation across 18 annotations × 22 chromosomes

## Technology stack

| Area | Tools |
|------|-------|
| Language | Python 3.11 |
| Single-cell | scanpy, anndata, snapATAC2 |
| Genomics | pybedtools, pyliftover, bedtools, scipy.sparse |
| Statistics | scikit-learn (LR + KNN), scipy.stats, statsmodels |
| Heritability | ldsc (Python 2.7 via Rosetta on Apple Silicon) |
| Figures | matplotlib (publication-ready PDF + PNG) |
| Reproducibility | conda/mamba environments, CLI-driven scripts |

## Scale

- ~2,500 lines of Python code across 19 scripts
- 4 publication-quality figures (multi-panel)
- 3 tables (TSV format, machine-readable)
- 63-row enrichment result table
- ~1,900-word narrative discussion with 16 citations
- Full reference data + GWAS sumstats download process documented end-to-end

## What this shows I can do

- Design a multi-method computational biology benchmark from scratch
- Read five independent method papers and implement their algorithms faithfully in a shared framework
- Bridge wet-lab biology (single-cell multi-omics) with statistical genetics (heritability partitioning)
- Identify and quantify a non-obvious failure mode in a standard integration pipeline
- Communicate technical work with figures, tables, and prose calibrated to the audience
