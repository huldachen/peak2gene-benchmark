# Discussion

## Summary

Five peak–gene linkage methods — distance-window (Method 1), Cicero co-accessibility (Method 2; Pliner et al., 2018), Activity-by-Contact (Method 3; Fulco et al., 2019), ArchR-style paired correlation (Method 4A; Granja et al., 2021), cross-cohort anchor-transfer correlation (Method 4B; Stuart et al., 2019), and the ENCODE-rE2G supervised model (Method 5; Gschwind et al., 2025) — were benchmarked on the Hickey et al. (2023) colon atlas by their ability to partition inflammatory bowel disease (IBD; de Lange et al., 2017) heritability using stratified LD-score regression (Finucane et al., 2015) with the baseline-LD v2.2 model (Gazal et al., 2017). Height (Yengo et al., 2022) and educational attainment (Rietveld et al., 2013) were included as non-gut negative controls.

Across 63 (method × compartment × trait) runs, 40 reached p ≤ 0.05 and 29 reached p ≤ 10⁻⁴. Three findings stand out.

## 1. Paired Multiome correlation concentrates heritability per SNP

Method 4A — direct peak–gene Pearson correlation on meta-cells from paired multiome nuclei — produced the sharpest per-SNP IBD enrichment: **21.0 ± 8.7× (epithelial, p = 0.02), 10.7 ± 4.7× (immune, p = 0.04), 14.3 ± 5.3× (stromal, p = 0.01)** in just 0.2–0.5 % of HapMap3 SNPs. Each annotated SNP carries roughly 3–5× the IBD heritability of a distance-baseline SNP. The trade-off is wide standard errors from the small annotation.

When paired multiome data is available, Method 4A is the best choice for analyses where per-SNP prior on causal status matters (e.g., heritability-informed fine-mapping; Weissbrod et al., 2022). It is not an option for atlases without paired measurements — which is where the next finding matters.

## 2. Cross-cohort integration fails at cell-type resolution

The practical alternative when paired data is unavailable is cross-cohort anchor-transfer integration (Method 4B), as implemented in Seurat/Signac (Stuart et al., 2019) and ArchR's `addGeneIntegrationMatrix()` (Granja et al., 2021). At the compartment-aggregate level, Method 4B's IBD enrichment looks reasonable (13.1–21.6× for immune and epithelial).

A per-cell diagnostic tells a different story. For each ATAC cell, we asked what fraction of its top-30 RNA anchors share its cell-type label. Results:

- Mean anchor agreement across compartments: **0.6–9.7 %**.
- 41 of 44 analysed cell types show < 10 % agreement; 29 show zero.
- Only a handful of cell types with very distinctive transcriptomes are anchored correctly: mature Enterocytes (61 %), Myofibroblasts/SM DES-high (68 %), Glia (29 %).

The gene-activity bridge feature is dominated by housekeeping-gene accessibility variance rather than cell-type-distinguishing regulation, so the KNN returns housekeeping-matched rather than cell-type-matched neighbours. Meta-cell aggregation rescues compartment-level signals by averaging across the within-compartment anchor noise, but any **cell-type-specific** downstream claim from Method 4B is unsupported by the underlying integration.

Most published peak–gene benchmarks do not explicitly diagnose anchor quality at the cell-type level. Our result suggests this diagnostic should be routine, especially when reference cohorts are small (here, 11 k RNA cells vs 102 k ATAC cells) or when compartments span a continuous transcriptional gradient (as in the gut epithelium). Alternative integration strategies — scVI (Gayoso et al., 2022), GLUE (Cao & Gao, 2022), CCA-MNN (Stuart et al., 2019) — may partially recover cell-type structure but deserve the same diagnostic.

## 3. Supervised-method transfer does not carry off-the-shelf

Method 5's pretrained `atac_megamap` rE2G model was calibrated on K562 CRISPR data with a 0.179 decision threshold for 70 % recall (Gschwind et al., 2025). Applied directly to colon: 1.6 M of 3.5 M candidate pairs pass threshold (46 % pass rate), producing an annotation that covers nearly every peak. The resulting sLDSC enrichment matches Method 1 (distance baseline) to within 0.01× across all nine trait × compartment cells. A continuous-valued annotation using the max rE2G probability per SNP gives the same result — the probability distribution is bimodal near 1, so continuous ≈ binary.

The model's *ranking* of pairs still contains information (features like ABC.Score and numTSSEnhGene generalise), but the *decision threshold* is miscalibrated for non-training tissue. Rank-based downstream use (top-N % by probability) is more defensible than binary calls at the K562 threshold until cell-type-matched CRISPR data enables recalibration.

## Supporting observations

**ABC with power-law Hi-C is structurally promoter-biased.** Median enhancer–TSS distance is 15–25 kb regardless of candidate window (±500 kb or ±5 Mb). Per-gene normalisation combined with `(d + 5000)^(−0.87)` decay means distal peaks cannot clear the 0.02 threshold unless their activity is disproportionately large. Real Hi-C would partially flatten this (Nasser et al., 2021) but the bias is structural to the proxy, not an implementation artefact. Published ABC analyses reporting short median distances should read these as method properties, not biological measurements of enhancer geometry.

**Specificity–coverage trade-off.** Methods 1–3 and Method 5 annotate 3–9 % of SNPs and capture 11–40 % of IBD heritability with 2–7× per-SNP enrichment; Methods 4A and 4B annotate 0.2–0.5 % of SNPs at 10–21×. Neither end dominates — choose based on downstream use. High-specificity Method 4A is better for fine-mapping priors; higher-coverage Methods 1–2 are better for broad tissue-of-action mapping (Finucane et al., 2018).

**Negative controls validate the pipeline.** Height enrichment concentrates in stromal (5.5× vs 2.2× immune for Method 1, p ≤ 10⁻⁸) — consistent with bone/connective-tissue biology. Educational attainment shows no compartment-specific signal (all p > 0.05) as expected for a brain-driven trait. These rule out pipeline artefacts.

## Limitations

- European ancestry only (reference LD + GWAS sumstats); extension to other ancestries requires ancestry-matched panels.
- Compartment-level pooling; finer cell-type resolution is supported only for Method 4A (paired data), not for methods that depend on cross-cohort integration.
- No cell-type-matched CRISPR ground truth for colon — method rankings are based on heritability-partitioning efficiency, an indirect metric.
- Hi-C proxy in Methods 3 and 5 is power-law, not real Hi-C, at an estimated 10–15 % AUPRC cost per Nasser et al. (2021).

## Future directions

- **Cell-type-specific claims**: rerun Method 4B with direct cell-type-label anchoring or an alternative integration method, each subject to the anchor-agreement diagnostic introduced here.
- **Method 3 sensitivity**: rerun with the Nasser et al. (2021) Avg-Hi-C cooler to isolate the contribution of the power-law proxy.
- **Method 5 recalibration**: extend the Gschwind et al. (2025) HCT116 CRISPR screen to colonocytes for gut-specific threshold calibration.
- **Broader traits**: add colorectal cancer, celiac disease, and irritable bowel syndrome to test generalisability.

## The contribution

A reproducible five-method benchmark that produces a defensible ranking (Method 4A > 2, 3 > 1 ≈ 5 > 4B for per-SNP IBD heritability), plus a quantitative diagnostic that reveals a silent failure mode in standard cross-cohort integration pipelines. The diagnostic costs ~30 lines of code and ~30 minutes of runtime, and should become part of standard practice.

---

## References

- Cao, Z.-J. & Gao, G. *Nat. Biotechnol.* **40**, 1458 (2022).
- de Lange, K. M. *et al.* *Nat. Genet.* **49**, 256 (2017).
- Finucane, H. K. *et al.* *Nat. Genet.* **47**, 1228 (2015).
- Finucane, H. K. *et al.* *Nat. Genet.* **50**, 621 (2018).
- Fulco, C. P. *et al.* *Nat. Genet.* **51**, 1664 (2019).
- Gayoso, A. *et al.* *Nat. Biotechnol.* **40**, 163 (2022).
- Gazal, S. *et al.* *Nat. Genet.* **49**, 1421 (2017).
- Granja, J. M. *et al.* *Nat. Genet.* **53**, 403 (2021).
- Gschwind, A. R. *et al.* *Nature* (2025).
- Hickey, J. W. *et al.* *Nature* **619**, 572 (2023).
- Nasser, J. *et al.* *Nature* **593**, 238 (2021).
- Pliner, H. A. *et al.* *Mol. Cell* **71**, 858 (2018).
- Rietveld, C. A. *et al.* *Science* **340**, 1467 (2013).
- Stuart, T. *et al.* *Cell* **177**, 1888 (2019).
- Weissbrod, O. *et al.* *Nat. Genet.* **54**, 450 (2022).
- Yengo, L. *et al.* *Nature* **610**, 704 (2022).
