# Results inventory

Quick index of everything in `results/`. See [`DISCUSSION.md`](DISCUSSION.md) for interpretation; the main [`README`](../README.md) embeds all four figures with captions.

## Figures

Each figure is available in **PNG** (for presentations) and **PDF** (for publication).

| # | File | What it shows |
|---|------|---------------|
| 1 | [`figures/fig1_data_overview`](figures/fig1_data_overview.png) | Two-cohort dataset structure, cell-count composition, RNA UMAP, per-modality QC metrics |
| 2 | [`figures/fig2_data_flow`](figures/fig2_data_flow.png) | How RNA and ATAC cohorts flow into each of the 5 linkage methods + question each answers |
| 3 | [`figures/fig3_ldsc_benchmark`](figures/fig3_ldsc_benchmark.png) | Forest plot of sLDSC heritability enrichment for IBD + 2 negative controls (main result) |
| 4 | [`figures/fig4_anchor_agreement`](figures/fig4_anchor_agreement.png) | Per-cell-type anchor-agreement showing cross-cohort integration failure mode |

Regenerate with `python scripts/figures/figure{1,2,3,4}_*.py`.

## Tables

All tables in **TSV**, **Markdown**, and **LaTeX** (`results/tables/`).

| # | File | What it contains |
|---|------|------------------|
| 1 | `table1_ldsc_benchmark` | 18 rows × 3 traits: enrichment + significance stars per (method, compartment) |
| 1b | `table1_ldsc_benchmark_extended` | Same as Table 1 with explicit p-value columns |
| 2 | `table2_method_design` | One row per method: inputs, algorithm, key parameter, threshold, source |
| 3 | `table3_anchor_agreement` | Per-compartment anchor-agreement with notable winner/loser cell types |

Regenerate with `python scripts/figures/make_tables.py`.

## Raw outputs

| Path | Content |
|------|---------|
| `ldsc/enrichment_table.tsv` | 63-row raw sLDSC output (every trait × method × compartment) |
| `linkage/method{1–5}/{compartment}.tsv` | Per-method peak–gene linkage tables |
| `linkage/method4b_anchor_diagnostic/{compartment}*.tsv` | Per-cell anchor-agreement values feeding Figure 4 |

## Quick reference (IBD primary trait)

| Method | Best-compartment enrichment | Caveat |
|--------|:---------------------------:|--------|
| Method 4A (paired Multiome) | **21.0×** (epithelial) | requires true paired data |
| Method 4B (cross-cohort) | 21.6× (epithelial) | cell-type integration fails (see Figure 4) |
| Method 2 (Cicero) | 7.1× (epithelial) | threshold not externally calibrated |
| Method 3 (ABC) | 7.2× (epithelial) | promoter-biased (power-law Hi-C) |
| Method 1 (distance) / 5 (rE2G) | 4.5× (epithelial) | rE2G's K562 threshold collapses onto distance in gut |

See [`DISCUSSION.md`](DISCUSSION.md) for full interpretation.
