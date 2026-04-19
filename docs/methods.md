# Methods — an educational walkthrough

Each of the five peak–gene linkage methods answers a slightly different question about enhancer–gene regulation. This document describes what each method is actually computing, in plain language, with enough technical detail to follow the code.

---

## Shared setup

- **Input peaks**: ATAC-seq peak calls per compartment (epithelial, immune, stromal), aggregated across ~100 k cells.
- **Input genes**: GENCODE v44 protein-coding gene annotations (hg38).
- **Candidate window**: ±500 kb around each gene's transcription start site (TSS) — the standard cis-regulatory window in the literature.
- **Output**: a standardised TSV per (method × compartment) with columns: `peak_id, gene_id, gene_name, chrom, peak_start, peak_end, tss_pos, strand, distance_bp, compartment, method, score`.

Having a shared output schema means every downstream step (sLDSC annotation generation, figure code, result aggregation) is identical across methods — only the scoring logic differs.

---

## Method 1 — Distance window

**Question answered:** *Is this peak inside the gene's cis-regulatory window?*

**Algorithm:** For each protein-coding gene, identify all ATAC peaks whose centre falls within ±500 kb of the TSS. Every such pair is a link; every link gets score 1.0.

**When it's the right method:** Whenever the downstream analysis just needs "SNPs in the vicinity of expressed genes" — classic sLDSC cell-type-specificity analyses do exactly this. It's also the natural null comparator against which fancier methods must demonstrate added value.

**Implementation:** `bedtools slop` to extend each TSS to a ±500 kb window, then `bedtools intersect` to find overlapping peaks. Output ~3 M links for immune compartment.

**File:** [`scripts/linkage/method1_distance.py`](../scripts/linkage/method1_distance.py)

---

## Method 2 — Cicero-style co-accessibility

**Question answered:** *Do these two peaks open together across cells?*

**Algorithm:**
1. Binarise the peak × cell matrix.
2. Compute a low-dimensional embedding via latent semantic indexing (TF-IDF → truncated SVD), dropping component 1 which correlates with sequencing depth.
3. Build KNN meta-cells: sample 2,000 random seed cells; each seed's meta-cell is the column-sum of its 50 nearest neighbours in LSI space. This recovers the statistical power of aggregation while preserving cell-type structure.
4. Standardise the (meta-cell × peak) matrix column-wise.
5. For each peak, compute Pearson correlation with every other peak within ±500 kb on the same chromosome. Threshold at r ≥ 0.25.
6. Convert peak–peak co-accessibility to peak–gene links: if either peak in a pair overlaps a TSS-proximal region (±2 kb), the other peak is linked to that gene. Promoter peaks also get a self-link at score 1.0.

**When it's the right method:** ATAC-only data where RNA is unavailable. Co-accessibility is a well-validated regulatory signal at the open-chromatin level.

**Citation:** Pliner et al. *Mol. Cell* **71**, 858–871 (2018).

**Implementation note:** The canonical R Cicero applies a graphical-lasso shrinkage with a distance-dependent penalty. This Python implementation emits raw Pearson correlations on the meta-cell aggregate — equivalent to what ArchR's `getCoAccessibility()` produces and what most downstream consumers actually use. The threshold is accordingly scaled (0.25 on raw correlations vs. Cicero's 0.05 on shrunk scores).

**File:** [`scripts/linkage/method2_cicero.py`](../scripts/linkage/method2_cicero.py)

---

## Method 3 — Activity-by-Contact (ABC)

**Question answered:** *Does peak activity, weighted by 3D-genome contact with the gene's promoter, exceed a threshold?*

**Algorithm:**
$$\mathrm{ABC}(E, G) = \frac{\mathrm{Activity}(E) \times \mathrm{Contact}(E, G)}{\sum_{e \in \text{candidates}} \mathrm{Activity}(e) \times \mathrm{Contact}(e, G)}$$

where **Activity** is the pseudobulk ATAC signal at the enhancer peak (ATAC-only variant — the original formula uses geometric mean of ATAC and H3K27ac, available when H3K27ac is measured), and **Contact** is approximated by the Nasser 2021 power-law `(d + 5000)^(−0.87)` since real Hi-C is not available for the Hickey cohort.

Per-gene normalisation (dividing by the sum over all candidate enhancers within ±5 Mb) is the key innovation of ABC: it makes scores calibrated across genes and lets one universal threshold (0.02) apply.

**When it's the right method:** Whenever a physical-contact interpretation of enhancer–promoter linkage is desired and Hi-C (or a reasonable proxy) is available.

**Citations:** Fulco et al. *Nat. Genet.* **51**, 1664 (2019); Nasser et al. *Nature* **593**, 238 (2021).

**Implementation note:** The power-law contact proxy produces a systematic bias toward promoter-proximal links (median distance 15–25 kb). This is inherent to the proxy and would be partially corrected by using a real Hi-C cooler.

**File:** [`scripts/linkage/method3_abc.py`](../scripts/linkage/method3_abc.py)

---

## Method 4A — Paired-Multiome correlation

**Question answered:** *Does peak accessibility covary with gene expression in the same cells?*

**Algorithm:**
1. Match barcodes between the Multiome RNA h5ad and the per-compartment ATAC h5ad — 100% pairing because Multiome captures both modalities from the same nuclei.
2. LSI embed the binarised ATAC matrix (same as Method 2).
3. Build 2,000 KNN meta-cells — importantly, the **same** meta-cell indicator is used for both modalities, so meta-cell *i* represents the same set of cells in both ATAC and RNA.
4. Aggregate: binary ATAC accessibility → meta-cell sum; RNA counts → log1p(meta-cell sum). Standardise each column.
5. For each gene expressed in ≥ 5 meta-cells, find all peaks within ±500 kb of its TSS. Compute Pearson r between the standardised peak column and gene column across meta-cells.
6. Convert r to a t-statistic, then to a p-value with n_meta − 2 degrees of freedom. Apply Benjamini–Hochberg FDR correction per-gene.
7. Keep links with |r| ≥ 0.45 AND FDR ≤ 0.05.

**When it's the right method:** The gold standard when paired multiome data is available. Every peak–gene link is supported by observed co-variation in the same nuclei.

**Citation:** Equivalent to ArchR's `addPeak2GeneLinks()` (Granja et al. *Nat. Genet.* **53**, 403, 2021).

**File:** [`scripts/linkage/method4_paired.py`](../scripts/linkage/method4_paired.py)

---

## Method 4B — Cross-cohort anchor transfer

**Question answered:** *If paired RNA data is unavailable, can we still estimate peak–gene correlation by borrowing RNA from a different cohort?*

**Algorithm:**
1. Compute a gene-activity matrix from ATAC: for each gene, sum the ATAC signal across peaks overlapping its body + 2 kb upstream.
2. Use gene-activity as a **bridge feature** to link the two cohorts in a shared space.
3. For each ATAC cell, find its *k* = 30 nearest RNA-cohort cells (cosine distance in the normalised gene-activity space). These are the "anchor" cells.
4. Impute RNA expression for each ATAC cell as a Gaussian-weighted average over its 30 anchors' RNA profiles.
5. With the ATAC cell now carrying imputed RNA, run the same meta-cell aggregation and peak–gene correlation as Method 4A.

**When it's the right method:** When ATAC is from one cohort and RNA is from another (no paired multiome available). It's the standard cross-cohort workflow in Seurat, Signac, ArchR, and SCENIC+.

**Citation:** Stuart et al. *Cell* **177**, 1888 (2019) describes the anchor-transfer formalism.

**A key caveat surfaced by this benchmark:** The anchor-transfer integration can look reasonable at the compartment-aggregate level while failing at cell-type resolution. See [`docs/results_summary.md`](results_summary.md) and `scripts/linkage/method4b_anchor_diagnostic.py` for the anchor-agreement diagnostic that reveals this.

**File:** [`scripts/linkage/method4_crosscohort.py`](../scripts/linkage/method4_crosscohort.py)

---

## Method 5 — ENCODE-rE2G (supervised)

**Question answered:** *Would a supervised machine-learning model, trained on thousands of CRISPR-validated enhancer–gene pairs, predict this as a regulatory link?*

**Algorithm:**
1. For every candidate (peak, gene) pair within ±500 kb, compute 8 features:
   - `numTSSEnhGene`: number of TSSs between the peak and the target
   - `distanceToTSS`: |peak_midpoint − tss|
   - `normalizedATAC_prom`: quantile-normalised pseudobulk ATAC at the gene's promoter peak
   - `sumNearbyEnhancers`: sum of peak activity within ±5 kb of the candidate peak
   - `ubiquitousExpressedGene`: 1 if the gene is a housekeeping gene (per rE2G's lookup table)
   - `numCandidateEnhGene`: number of peaks between the peak and the target
   - `contactFrequency`: power-law Hi-C proxy
   - `ABC.Score`: per-gene-normalised Activity × Contact
2. Log-transform the feature matrix (`np.log(|X| + 0.01)` — the same transformation used in rE2G's training pipeline).
3. Feed to the pretrained `atac_megamap` logistic-regression model from ENCODE-rE2G.
4. Apply the authors' binary threshold of 0.179.

**When it's the right method:** Whenever one trusts supervised calibration against CRISPR truth to generalise to a target cell type — and when the target cell type isn't far from the training cell type (K562).

**Citation:** Gschwind A. R. et al. *Nature* (2025).

**Implementation note:** The pretrained model is a scikit-learn pickle (`atac_megamap/model.pkl`) from the EngreitzLab repo. Our implementation computes the 8 required features natively, applies rE2G's training-time log-transform, and calls `model.predict_proba()`. Continuous-score (Prob ∈ [0, 1]) and binary (Prob ≥ 0.179) variants are both produced.

**File:** [`scripts/linkage/method5_re2g.py`](../scripts/linkage/method5_re2g.py)

---

## Cross-cohort integration diagnostic

**Not a linkage method, but a quality-control tool**: `scripts/linkage/method4b_anchor_diagnostic.py` tests whether each ATAC cell's top-30 RNA anchors share its cell-type label. This operationalises a natural quality-control question ("is the integration putting the right cells together?") that is rarely explicitly tested in published cross-cohort pipelines.

**Output:** per-compartment + per-cell-type anchor-agreement rates. The results reveal that most cell types in our dataset have near-zero agreement, with only a handful of cell types (Enterocytes, Myofibroblasts/SM DES High, Glia) showing strong agreement. See [`docs/results_summary.md`](results_summary.md) and Figure 4.

---

## Comparison matrix

| | Method 1 | Method 2 | Method 3 | Method 4A | Method 4B | Method 5 |
|--|:---:|:---:|:---:|:---:|:---:|:---:|
| Needs ATAC | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Needs RNA | — | — | — | ✓ (paired) | ✓ (any) | — |
| Needs Hi-C | — | — | proxy OK | — | — | proxy OK |
| Supervised | — | — | — | — | — | ✓ |
| Cell-type validity | n/a | aggregate | aggregate | ✓ | ⚠ diagnostic req'd | transfer-limited |
| Typical # links (immune) | 3.5 M | 600 K | 150 K | 48 K | 230 K | 1.6 M |
| Median distance | 230 kb | 170 kb | 25 kb | 130 kb | 230 kb | 120 kb |

The ladder runs from loose (Method 1) to strict (Method 4A) in specificity, and conversely from high to low in genome coverage. See [`results/DISCUSSION.md`](../results/DISCUSSION.md) for how to pick between them for specific downstream analyses.
