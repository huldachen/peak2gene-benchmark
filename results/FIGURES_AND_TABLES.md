# Results deliverables - Hickey colon peak-gene linkage benchmark

All figures, tables, and writeup sections for the paper / thesis chapter.

- **Figures** in PNG (for presentations) and PDF (for publication).
- **Tables** in TSV, Markdown, and LaTeX.
- **Discussion** in Markdown (`DISCUSSION.md`).

See `DISCUSSION.md` for the narrative writeup (~1,900 words, 16 citations).

---

## Figures

Four figures total:

| Fig | Name | Purpose |
|-----|------|---------|
| 1 | Dataset QC + composition | "What data do we have?" |
| 2 | Data flow schematic | "How do methods use the data?" |
| 3 | sLDSC benchmark forest plot | **Main finding** - method rankings |
| 4 | Anchor-agreement bar chart | **Finding E visual** - integration failure |

---

### Figure 1 - Dataset QC + composition
**File:** `results/figures/fig1_data_overview.{png,pdf}`
**Generator:** `scripts/figures/figure1_data_overview.py`

Four-panel overview of the Hickey 2023 colon single-cell dataset:

- **Panel A** - Two-cohort schematic: RNA cohort (B001/B004/B005, 3 donors, 4 colon samples, 11,604 cells) vs Multiome cohort (B006/B008-B012, 6 donors, 42 samples, 102,453 cells).
- **Panel B** - Cells per compartment x cohort. Stacked bars showing epithelial/immune/stromal distribution in each cohort.
- **Panel C** - RNA UMAP coloured by cell type (top-12 types in legend).
- **Panel D** - Four QC violins: RNA n_genes, RNA %mito, ATAC log₁₀(n_peaks), ATAC log₁₀(total_accessibility).

### Figure 2 - Data flow across the 6 methods
**File:** `results/figures/fig2_data_flow.{png,pdf}`
**Generator:** `scripts/figures/figure2_data_flow.py`

Schematic showing how RNA and ATAC data flow from the two cohorts into the six peak-gene linkage methods. Colour-coded arrows indicate integration quality (green = paired/valid, amber = compartment-only, red = integration fails, violet = external prior).

Each method card shows name, algorithm summary, and final link count. A bottom-row ribbon lists the specific **question each method answers** (e.g., Method 1: "Is the peak within the cis-regulatory window of the gene?").

### Figure 3 - sLDSC benchmark (headline result)
**File:** `results/figures/fig3_ldsc_benchmark.{png,pdf}`
**Generator:** `scripts/figures/figure3_ldsc_benchmark.py`

Two-row panel layout (Row A: enrichment +/- SE; Row B: tau\* coefficient z-score) x three traits (IBD primary, Height negative control, EA negative control). 21 rows per trait-panel (7 method keys x 3 compartments - Methods 1, 2a, 2b, 3, 4A, 4B, 5) showing enrichment +/- SE with significance stars in Row A and the size-invariant tau\* in Row B. Dot size encodes `prop_snps` (annotation coverage), so the reader sees the specificity-vs-coverage trade-off at a glance. Compartments are colour-coded (orange=epithelial, blue=immune, purple=stromal). Dashed vertical line at 1x marks the null in Row A; at 0 marks the null in Row B.

The main takeaway reads off directly:
- **Row A:** Method 4A (paired Multiome) enrichment dots sit far to the right (21x for IBD epithelial) but with wide error bars (small annotation -> higher SE); Method 2b (canonical Cicero with GLASSO) sits one tier lower (6-8x) with intermediate annotation size; Methods 1, 2a, 3 cluster near 3-7x with tight SEs; Method 5 overlaps Method 1 exactly because the K562 threshold isn't calibrated for gut tissue.
- **Row B:** Method 4A is the only method with tau\* > 2 (significant unique contribution beyond baseline-LD); Method 2a's Height-immune cell sits at tau\* = −3.44 (annotation redundant with baseline-LD), while Method 2b's same cell is tau\* = +0.85 - concrete evidence the canonical Cicero step removes baseline-LD-redundant edges that the Pearson shortcut retains.

### Figure 4 - Anchor-agreement diagnostic (Finding E visual)
**File:** `results/figures/fig4_anchor_agreement.{png,pdf}`
**Generator:** `scripts/figures/figure4_anchor_agreement.py`

Three-panel horizontal bar chart (one per compartment) showing per-cell-type anchor-agreement. Bars coloured green (>=50 %), amber (10-50 %), red (<10 %). The pattern is striking: most cell types are red bars near zero, with only a handful of green "winners" (Enterocytes 61 %, Myofibroblasts/SM DES High 68 %, etc.). Overall mean agreement per compartment is listed in each panel title (epithelial 9.7 %, immune 0.6 %, stromal 7.6 %).

This figure operationalises the Finding E caveat: *most cell types have zero agreement despite Method 4B's aggregate metrics looking reasonable.*

---

## Tables

All tables are in `results/tables/`. Each table is available in 3 formats: `.tsv`, `.md`, `.tex`.

### Table 1 - sLDSC benchmark results (headline finding)
**Files:** `table1_ldsc_benchmark.{tsv,md,tex}` and `table1_ldsc_benchmark_extended.{tsv,md,tex}`

21 rows (7 methods x 3 compartments) x 3 traits (IBD, Height, EA). Each cell: enrichment +/- SE with significance stars (\* p<0.05, \*\* p<0.01, \*\*\* p<0.001).

The extended version adds explicit p-value columns.

### Table 2 - Method design reference
**Files:** `table2_method_design.{tsv,md,tex}`

One row per method. 8 columns: Method, Inputs, Algorithm, Key parameter, Threshold, Canonical source, Our setting, Defensibility. Serves as a standalone methods-supplementary reference.

### Table 3 - Cross-cohort anchor-agreement diagnostic (Finding E quantitative)
**Files:** `table3_anchor_agreement.{tsv,md,tex}`

Per-compartment + per-cell-type anchor-agreement rates. For each compartment:
- One summary row with weighted-mean agreement across all cells
- Top 3 most abundant cell types
- Any additional "winners" (mean agreement >= 0.25)

Highlights the two non-obvious findings: (a) mostly-zero agreement across compartments, (b) only a handful of cell types with distinctive markers (Enterocytes, Myofibroblasts/SM DES High, Glia) get correct anchors.

### Table 4 - Total heritability captured (`prop_h2` supplement)
**Files:** `table4_prop_h2_supplement.{tsv,md,tex}`

21 rows (7 method keys x 3 compartments) x 3 traits (IBD, Height, EA), columns: Method, Compartment, prop_snps, IBD prop_h2, Height prop_h2, EA prop_h2. Complements Table 1's per-SNP enrichment view with the *total coverage* view: how much of each trait's heritability does each annotation actually capture.

Highlights the size-coverage trade-off:
- Method 1 (~8% of SNPs) captures 28-40 % of IBD h^2 - broad coverage, dilute per SNP.
- Method 4A (~0.25-0.54 % of SNPs) captures only 5-6 % of IBD h^2 - concentrated per SNP, low total coverage.
- Method 2b (~1.1 % of SNPs) sits between, capturing 6-9 %.
- Method 5 ≈ Method 1 within rounding (the K562 threshold passes nearly every candidate region for gut tissue).
- Method 4B stromal: prop_h2 = −0.3 % (worse than null), corroborating the Finding-E anchor collapse.

---

## Quick reference: what the benchmark shows

From Table 1 for the primary trait (IBD):

| Method              | Per-SNP enrichment (best compartment) | Caveat |
|---------------------|--------------------------------------:|--------|
| Method 4A (paired)  | **21.0x** (epithelial)                | true paired Multiome required |
| Method 4B (cross-cohort) | 21.6x (epithelial)               | Finding E: cell-type anchors mostly fail |
| Method 2 (Cicero)   | 7.1x (epithelial)                      | threshold not externally calibrated |
| Method 3 (ABC)      | 7.2x (epithelial)                      | promoter-biased (power-law Hi-C) |
| Method 1 / 5        | 4.5x (epithelial)                      | rE2G indistinguishable from distance at K562 threshold |

- **Best specificity/SNP**: Method 4A - 10-21x in just 0.2-0.5% of SNPs.
- **Best coverage + signal**: Methods 1-2 - 3-9% of SNPs, 25-40% of heritability.
- **Strongest negative control**: Height correctly enriches stromal (5-8x, p <= 10⁻⁹) across Methods 1-3.

See `lab_notebook/entries/2026-04-16.md` for the full discussion, decisions (DEC-014-020), and findings (A-E).
