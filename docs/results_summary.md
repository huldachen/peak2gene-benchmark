# Results summary

A plain-English walkthrough of the benchmark findings. For full technical detail see [`results/DISCUSSION.md`](../results/DISCUSSION.md). For raw numbers see [`results/tables/table1_ldsc_benchmark.md`](../results/tables/table1_ldsc_benchmark.md).

---

## One-sentence headline

When paired multi-omic data is available, **direct peak–gene correlation (Method 4A) is the most informative per-SNP method for partitioning IBD heritability** — 10–21× enrichment in under 0.5 % of the genome — while simpler methods capture a larger total share of heritability at lower per-SNP density.

## The benchmark result (Figure 3)

<p align="center">
  <img src="../results/figures/fig3_ldsc_benchmark.png" alt="sLDSC benchmark forest plot" width="920"/>
</p>

Reading the figure:
- **IBD (left panel)**: every method reaches significance (p ≤ 0.05) in at least one compartment. Method 4A's dots are farthest to the right at 10–21×, with wider error bars because the annotation is small. Methods 1–3 cluster around 3–7× with tight SEs.
- **Height (middle panel, negative control)**: stromal compartment correctly dominates (5–8× for Methods 1–3, p ≤ 10⁻⁸) — consistent with bone/connective-tissue biology. Method 4A/4B lose the signal because paired-cell correlation is not the dominant mechanism for height variation.
- **EA (right panel, second negative control)**: all methods show weak, non-significant enrichments, as expected for a trait driven by brain tissue not represented in the annotations.

## Method rankings by question asked

Different downstream questions favour different methods. Rough guidance:

| Downstream use | Method to pick | Why |
|-----------------|:--------------:|-----|
| Fine-mapping prior (per-SNP specificity matters) | **Method 4A paired Multiome r** | 21× enrichment in 0.2 % of SNPs |
| Tissue-of-action map across many traits | Methods 1–2 | Broadest annotation, best recall per compartment |
| Quick baseline / no RNA available | Method 1 distance | Zero RNA cost, baseline signal |
| Proximal-enhancer discovery | Method 3 ABC | Biased to cis but captures activity-weighted links |
| Cell-type-specific claims without paired data | ⚠ Run the anchor-agreement diagnostic first | See below |

## Finding: cross-cohort integration can silently fail at cell-type resolution

**This is the methodological caveat that a disease-gene discovery pipeline needs to know about.**

Many single-cell multi-omics pipelines face a version of this problem: you have ATAC from one cohort and RNA from another cohort, and you want to combine them to study cell-type-specific regulation. The standard approach (gene-activity bridge + KNN anchor transfer, used by Seurat, Signac, ArchR, SCENIC+) looks reasonable at the aggregate level but can fail at cell-type resolution in a way that standard quality controls don't catch.

To test this, we ran an **anchor-agreement diagnostic**: for every ATAC cell, find its top 30 RNA anchors and ask, what fraction of those anchors actually share the ATAC cell's cell-type label? If the integration is working, the answer should be high (most anchors should be the matching cell type).

**Result (Figure 4):**

<p align="center">
  <img src="../results/figures/fig4_anchor_agreement.png" alt="Anchor-agreement diagnostic" width="720"/>
</p>

For most cell types in our dataset, **anchor agreement is below 10 %** — effectively random. Only a handful of cell types with very distinctive transcriptional signatures (mature Enterocytes at 61 %, Myofibroblasts/SM DES High at 68 %, Glia at 29 %) have correct anchors. Everything else — the abundant transit-amplifying cells, progenitors, and closely related cell-type variants — is anchored to wrong cells.

### What this means in practice

- **Aggregate metrics can look OK while cell-type integration is broken**: Method 4B's sLDSC enrichment for IBD was 13–21× at the compartment level, which passed standard sanity checks. But at the cell-type level, the integration was not discriminating correctly — most ATAC cells were matched to RNA cells with the wrong identity.
- **The gene-activity bridge captures housekeeping variance more than cell-type variance**: ATAC signal at housekeeping genes is similar across cell types; ATAC signal at cell-type-specific regulatory regions is very different. If the bridge feature is dominated by housekeeping-gene accessibility patterns, the KNN will find "most similar housekeeping profile" rather than "most similar cell type." In this dataset, that's what happened.
- **Paired multiome, where available, avoids this entirely**: Method 4A uses same-cell barcode matching, which is always correctly "anchored."
- **Without paired data, the diagnostic should be run before downstream interpretation**: if most cell types have <10 % anchor agreement, any cell-type-specific claim is not supported by the integration.

This diagnostic is a small addition to the standard pipeline (~30 lines of code + 30 minutes to run) but provides empirical backing for whether cross-cohort integration is actually doing what it's supposed to.

## Finding: pretrained supervised models don't transfer off-the-shelf

The ENCODE-rE2G model (Gschwind et al. 2025) is a logistic regression trained on CRISPR-validated enhancer–gene pairs from K562 cells, with a decision threshold calibrated for 70 % recall on that training set. Applied to colon cells, the model's probability distribution saturates near 1 for most candidate pairs; 46 % of all candidate pairs pass the K562-calibrated threshold. The resulting binary annotation covers essentially every peak — making it numerically identical to the distance-window baseline in our sLDSC results.

Using a continuous-valued annotation (where each SNP's value is the maximum rE2G probability of any peak containing it) gives the same result, because the rE2G probability distribution is strongly bimodal (either near 0 or near 1) — so the continuous version behaves like binary.

**What this means in practice**: pretrained supervised models should be recalibrated on target-cell-type CRISPR data before binary thresholding. Until such data is available, the safer uses of Method 5 are (a) rank-based filtering (top N % by probability) or (b) incorporating the probability as one feature in a broader meta-analysis, not as a binary classifier.

## Finding: ABC model with power-law Hi-C has structural promoter bias

The ABC model's per-gene normalisation combined with a power-law contact term produces a median enhancer-TSS distance of 15-25 kb, much shorter than other methods produce. This isn't a bug in our implementation — it's inherent to how the power-law proxy weights distal candidates. Real Hi-C data would partially correct this; it's an important caveat when interpreting published ABC results that use power-law contact.

## Final pruning ladder

Across 5 methods × 3 compartments × 3 traits = 63 sLDSC runs, 40 reach p ≤ 0.05 and 29 reach p ≤ 10⁻⁴. The specificity–coverage trade-off is clean:

| Method | Cover this much of genome | Capture this much of IBD h² | Per-SNP enrichment |
|--|:-:|:-:|:-:|
| Distance window | 7.4 % | 31 % | 4.1× |
| Cicero co-access | 5.8 % | 23 % | 3.9× |
| ABC (power-law) | 4.1 % | 11 % | 2.7× |
| Cross-cohort r (4B) | 0.5 % | 7 % | 13.1× |
| **Paired Multiome r (4A)** | **0.5 %** | **6 %** | **10.7×** |

A team picking a method should ask: "which end of this trade-off serves my next step?" There's no universally best answer, but the benchmark makes the trade-off explicit.

---

## Takeaway

The project delivers three concrete deliverables for a disease-gene-discovery team:

1. **A ranked comparison** of five established peak–gene linkage methods on a common dataset and a biologically meaningful evaluation metric (IBD heritability partitioning).
2. **A quality-control diagnostic** (anchor-agreement) that reveals a previously under-recognised failure mode in standard cross-cohort integration pipelines.
3. **A reproducible pipeline** (single conda env, scripted figures and tables) that can be applied to any other multi-omic atlas + GWAS combination.

Full narrative and citations in [`results/DISCUSSION.md`](../results/DISCUSSION.md).
