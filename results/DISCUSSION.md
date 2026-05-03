# Discussion

## Summary

A team choosing a peak-gene linkage method today reads four or five method papers, each validated against its preferred CRISPRi or eQTL panel in its preferred cell type, and is left to guess which one transfers to their disease and tissue. The existing head-to-head benchmarks - Gschwind et al. (2025) for ENCODE-rE2G against ABC and distance baselines, Nasser et al. (2021) for ABC against several alternatives - test some of these methods on cell-line CRISPR perturbation data, predominantly K562 with smaller HCT116 and fine-mapped eQTL panels. **This benchmark complements those with a different lens**: a broader six-method scope (including Cicero canonical, ArchR's faster Pearson shortcut, paired-Multiome correlation, and cross-cohort anchor-transfer alongside ABC and rE2G), on a primary-tissue single-cell atlas, using disease-relevant GWAS heritability partitioning as the readout rather than CRISPRi gene-by-gene ground truth. The two ground truths probe different aspects of regulatory validity - CRISPRi tests element-level perturbation effects on a small number of genes; sLDSC tests whether the SNPs in an annotation explain a disproportionate share of trait heritability genome-wide - and a method-comparison conclusion that holds across both is more robust than one anchored on either alone.

We ran six methods on the Hickey et al. (2023) colon atlas - the distance-window baseline (Method 1), Pearson co-accessibility on KNN meta-cells equivalent to ArchR's `getCoAccessibility` (Method 2a; Granja et al., 2021), graphical-lasso partial-correlation co-accessibility (Method 2b; canonical Cicero, Pliner et al., 2018), Activity-by-Contact (Method 3; Fulco et al., 2019; Nasser et al., 2021), paired-Multiome correlation (Method 4A; Granja et al., 2021), cross-cohort anchor-transfer correlation (Method 4B; Stuart et al., 2019), and the supervised ENCODE-rE2G logistic-regression model (Method 5; Gschwind et al., 2025) - and ranked them by stratified LD-score regression (Finucane et al., 2015) partitioned heritability for inflammatory bowel disease (de Lange et al., 2017) on the baseline-LD v2.2 model (Gazal et al., 2017), with Height (Yengo et al., 2022) and Educational Attainment (Rietveld et al., 2013) as non-gut negative controls.

Across 72 (method x compartment x trait) sLDSC runs, 45 combinations reach p <= 0.05 and 31 reach p <= 10^-4. The ranking on its own is informative: Method 4A concentrates IBD heritability the most per SNP (10-21x in <=0.5% of tested SNPs), Methods 1, 2a, 3 cover a larger fraction of the genome at lower per-SNP enrichment (2-7x in 3-9%), Method 2b sits between these extremes (5-8x in ~1%), and Method 5 at the K562-calibrated threshold collapses numerically onto Method 1. **But the more useful contribution is what the comparison reveals about three claims that recur in the original publications and do not survive sLDSC scrutiny:** (i) ArchR's raw-Pearson shortcut for Cicero is *not* statistically equivalent to canonical GLASSO Cicero when conditioned on baseline-LD; (ii) anchor-transfer integration is *not* reliable at the cell-type level for most cell types in primary tissue; (iii) the K562-calibrated ENCODE-rE2G threshold *does not* transfer to gut tissue - the resulting annotation collapses onto a distance baseline. Each finding has a concrete, testable downstream consequence and is developed below.

## Paired Multiome correlation concentrates heritability per SNP

Method 4A, which computes per-gene Pearson correlations between peak
accessibility and gene expression on KNN-aggregated meta-cells drawn
from the multiome cohort (same nuclei, same barcodes), produced the
cleanest per-SNP heritability signal for IBD. Enrichment reached 21.0 +/-
8.7x in the epithelial compartment (p = 0.02), 10.7 +/- 4.7x in immune
(p = 0.04), and 14.3 +/- 5.3x in stromal (p = 0.01), despite annotating
only 0.2-0.5 % of HapMap3 SNPs. In relative terms, each annotated SNP
from Method 4A carries roughly 3-5x the heritability signal of an
annotated SNP from the distance baseline (Method 1), mirroring the
specificity gain reported in the original ArchR publication
(Granja et al., 2021) and consistent with the interpretation that
correlated peak-gene pairs enrich for genuine cis-regulatory activity.
Method 4A is therefore the strongest candidate when per-SNP enrichment
matters - for example, when constructing compact annotations for
heritability-informed fine-mapping (Weissbrod et al., 2022) - though the
small annotation size also produces correspondingly wide standard
errors. The same point is sharper in the tau\* (coefficient z-score) view
shown in Row B of Figure 3: tau\* measures the *unique* heritability signal
contributed beyond the baseline-LD v2.2 model and is invariant to
annotation size, which makes it the more rigorous comparison metric across
methods that span a 25x range in coverage. Method 4A's tau\* values exceed
their respective standard errors in all three compartments (i.e., the
annotation contributes signal not already captured by the ~97 features
of baseline-LD), while Methods 1, 2a, 3, and 5 cluster near or below zero
in most cells - confirming that the size-corrected per-SNP signal is
genuinely strongest in paired Multiome correlation rather than an
artefact of the small-annotation enrichment inflation seen in Row A.
A complementary view of total heritability *coverage* (prop_h2; Table 4)
shows the converse trade-off: Method 4A captures only 5-6 % of IBD
heritability per compartment despite its large per-SNP enrichment,
while Method 1 captures 28-40 %; the appropriate metric depends on the
downstream question.

Importantly, this benefit is only available on datasets with
true paired multiome; adapting the approach to atlases without paired
measurements requires either the cross-cohort integration we evaluate
next, or alternative bridges such as SCENIC+-style TF-target
reconstruction (Bravo González-Blas et al., 2023).

## ArchR's raw-Pearson shortcut is not statistically equivalent to canonical Cicero

Pliner et al. (2018) define Cicero as graphical lasso applied to overlapping 500 kb genome windows with a distance-based penalty: the regularisation `rho_{ij}` decreases for peaks closer in genomic distance, so nearby pairs can retain strong edges while distant pairs need stronger evidence. This is the canonical algorithm cited as "Cicero" in the multi-omics literature, and SCENIC+ (Bravo González-Blas et al., 2023) uses canonical Cicero internally for its co-accessibility step. ArchR's `getCoAccessibility()` (Granja et al., 2021) takes a faster path - Pearson correlation between meta-cell aggregated accessibility profiles, computed via the C++ `rowCorCpp` helper - and is widely used as a Cicero substitute when the canonical R package is impractical at scale. The downstream literature is comfortable treating the two paths as interchangeable: many analyses cite Cicero (Pliner et al., 2018) as the methodological reference while running ArchR's shortcut in practice.

**At the level of sLDSC heritability partitioning, they are not interchangeable.** Running both as separate methods on the same input - Method 2a (raw Pearson, ArchR-equivalent) at r >= 0.25, and Method 2b (GLASSO with alpha = 0.5 and a 250 kb distance scale) at partial-r >= 0.05 - reveals two reproducible, quantitative differences and one qualitative one.

The first is structural. Method 2b retains 7-14 % of Method 2a's edges (epithelial: 91 k of 1.20 M; immune: 85 k of 599 k; stromal: 79 k of 641 k). The GLASSO step is doing what Cicero claims it does - aggressively pruning indirect / mediated edges - and the surviving 7-14 % are the direct partial correlations. The second is a per-SNP enrichment shift. The smaller GLASSO annotation is more concentrated per SNP: IBD enrichment moves upward going from 2a to 2b across all three compartments, with an average shift of +1.1x and a maximum of +4.3x in the immune compartment (3.91 +/- 0.51x to 8.18 +/- 4.95x). The raw-Pearson approximation is *conservative* in the per-SNP enrichment sense - it under-reports how concentrated co-accessibility-driven heritability really is.

The third difference, and the one with the largest practical consequence, is in the tau\* coefficient (Row B of Figure 3). tau\* is the contribution of an annotation to heritability *beyond* the ~97 features of baseline-LD v2.2 (Gazal et al., 2017) - the standard size-invariant metric for asking "is this annotation telling us something baseline-LD doesn't already cover?" For Height in the immune compartment, Method 2a's tau\* is **-3.4** (highly significant negative): the raw-Pearson annotation is *redundant with baseline-LD*, and the model fits better when the annotation is down-weighted. For the same cell, Method 2b's tau\* is **+0.9** (benign): the GLASSO-pruned annotation no longer leaks into baseline-LD redundancy. GLASSO is removing precisely the indirect peak-peak edges that overlap with baseline-LD's existing annotation set - which means raw Pearson is silently importing those overlaps into the heritability model.

The headline ranking among method *families* (Method 4A > 2 >= 3 > 1 >= 5 > 4B for IBD per-SNP enrichment) is preserved across the 2a to 2b switch on this dataset. **But that rank stability is a property of the Hickey 2023 colon atlas combined with the de Lange 2017 IBD GWAS** - well-powered immune-relevant hits, dense immune cells, strong compartment specificity - **not of the algorithms.** On a noisier dataset (smaller cohort, weaker GWAS, less compartment specificity, or a less immune-driven trait), the 4.3x IBD-immune gap is plausibly large enough to flip a ranking; we cannot rule that out without replication on a second dataset. The honest reading is therefore narrow: on this benchmark, the choice between Methods 2a and 2b matters quantitatively for per-SNP enrichment and qualitatively for tau\* in at least one cell, but it does not change the rank order against the other linkage families. Anyone running tau\*-based or baseline-LD-conditioned analyses with ArchR's shortcut on a different dataset should run both methods and compare, rather than treating either as a universal substitute for the other.

## Cross-cohort integration fails at cell-type resolution (Finding E)

The central methodological contribution of this benchmark is a direct
measurement of the cost of substituting cross-cohort anchor-transfer
integration for true paired data. Method 4B implements a standard
approach - gene-activity as a bridge feature, cosine-distance KNN,
weighted expression transfer - that closely resembles the default
pipelines in Seurat/Signac (Stuart et al., 2019) and ArchR's
`addGeneIntegrationMatrix()` (Granja et al., 2021). At the compartment
level, Method 4B's IBD enrichment results look broadly reasonable:
13.1 +/- 5.2x for immune (p = 0.02), 21.6 +/- 10.4x for epithelial (p = 0.05).
However, when we asked whether each ATAC cell's top-30 RNA anchors
share its known cell-type label - a property the integration must
satisfy if it is to support cell-type-specific downstream claims - the
answer is that they do not (Table 3; Figure 4). Across all three
compartments, mean anchor agreement is 0.6-9.7 %: for 41 of the 44
analysed cell types, fewer than 10 % of anchors are correctly typed,
and 29 cell types have zero agreement at all. Only a handful of cell
types with very distinctive transcriptomes get correct anchors:
mature Enterocytes (61 %), Myofibroblasts/SM DES-high (68 %), Glia
(29 %). This pattern indicates that the gene-activity bridge feature
is dominated by housekeeping-gene accessibility variance rather than
cell-type-distinguishing regulation, and the KNN therefore returns
housekeeping-matched rather than cell-type-matched neighbours.

The implication is that the aggregate peak-gene correlations from
Method 4B reflect *compartment-level* rather than cell-type-level
covariation: meta-cell aggregation rescues the signal by averaging
across within-compartment noise, but any claim about a specific cell
type's regulatory architecture is not supported. Few published
peak-gene benchmarks explicitly diagnose anchor quality at the cell-
type level; our result suggests this diagnostic should become
routine, particularly when the reference cohort is small (here, 11 k
RNA cells vs 102 k ATAC cells) or cell-type diversity within a
compartment is high (as in the gut epithelium, which spans a
continuous stem -> mature gradient with 10+ transcriptional states).
Alternative integration strategies - deep joint latent spaces
(scVI; Gayoso et al., 2022), adversarial matching (GLUE; Cao & Gao,
2022), or CCA-MNN (Stuart et al., 2019) - may partially recover
cell-type structure but should be subject to the same anchor-
agreement sanity check before downstream cell-type-specific
interpretation.

## Supervised method transfer does not carry off-the-shelf

The pretrained rE2G `atac_megamap` logistic-regression model was
developed and calibrated on K562 CRISPR-validated enhancer-gene pairs
(Gschwind et al., 2025), with a threshold of 0.179 chosen to achieve
70 % recall on that ground truth. Applied directly to colon data, the
model's probability distribution saturates near 1.0 for most candidate
pairs: 1.6 M of 3.5 M immune candidates pass the K562-calibrated
threshold (46 % pass rate), yielding a binary annotation that covers
essentially every peak. The resulting sLDSC enrichment is numerically
identical to that of Method 1 (distance-window baseline) to within
0.01x across all nine trait x compartment cells (Table 1; Figure 3).
Rerunning the analysis with a continuous annotation (each SNP's
annotation value = max rE2G probability of any peak containing it)
did not change the enrichment - the probability distribution among
non-zero SNPs is strongly saturated near 1.0, so the continuous
formulation is effectively binary.

This is a concrete, quantitative illustration of a broader concern
regarding pretrained supervised models applied to new cell types:
threshold calibration and feature-scale distributions are both
cell-type-specific properties that do not automatically transfer.
In our setting, Method 5's *ranking* of pairs by probability still
contains information - the underlying logistic regression captures
features like ABC score, numTSSEnhGene, and normalizedATAC_prom that
do generalise - but the *decision threshold* is miscalibrated for
gut tissue, producing an over-permissive annotation. A gut-matched
CRISPR dataset (e.g., extending the Gschwind et al., 2025 HCT116 screen)
would enable recalibration; until then, rank-based downstream use
(e.g., top-N % by probability) is more defensible than binary calls
at the K562 threshold.

## Structural promoter bias in ABC with power-law Hi-C

The ABC model's per-gene normalisation combined with a power-law
contact term `(d + 5 000)^(−γ)` produces a median enhancer-TSS
distance of 15-25 kb in our output - consistent across compartments
and robust to the candidate window size (+/-500 kb vs +/-5 Mb). Because
the normalisation divides each peak's Activity x Contact product by
the sum across all candidates for that gene, and because the power-
law decays by a factor of ~67 across the 500 kb-TSS range, distal
peaks can only clear the 0.02 threshold when their activity is
disproportionately large. Real Hi-C typically has a longer-tailed
contact distribution and would partially flatten this bias (Nasser
et al., 2021), but the effect is structural to the proxy, not an
implementation artefact. Published ABC analyses that report
"median enhancer-TSS distance = X kb" should therefore be read as
a property of the method's contact term rather than a biological
measurement of enhancer geometry; in our data this shift was ~2x
(median 16 kb with 500 kb window, 8 kb with 5 Mb window). This
caveat has direct bearing on tissue-of-action inference from
ABC-style annotations: any conclusion that gut disease heritability
is "promoter-concentrated" requires independent validation with
Hi-C-aware methods or paired multiome correlation before being
accepted as biology rather than artefact.

## Specificity-coverage trade-off in heritability partitioning

Taken together, the eight method variants occupy distinct positions
on a specificity-coverage plane in heritability units (Table 1; Table 4;
Figure 3): Methods 1, 2a, 3, and the collapsed Method 5 annotate 3-9 % of
HapMap3 SNPs and capture 11-40 % of IBD heritability with per-SNP
enrichment of 2-7x; Method 2b sits at an intermediate point - ~1 % of
SNPs, 6-9 % of IBD heritability, 6-8x per-SNP enrichment - neither broad
nor maximally compact; Methods 4A and 4B annotate 0.2-0.5 % of SNPs with
per-SNP enrichment of 10-21x but capture only 5-7 % of heritability.
Neither end dominates; the appropriate method depends on the
downstream question. For heritability-informed fine-mapping, where
per-SNP prior on causal status is the primary input, the high-
specificity Method 4A annotation is preferable despite its small
footprint. For broad tissue-of-action mapping across many traits
(e.g., PanGWAS-style partitioning; Finucane et al., 2018), the
higher-coverage Methods 1 and 2 are likely more informative. Hybrid
strategies - using the union of multiple method annotations as a
joint sLDSC model - would quantify this trade-off directly and
deserve systematic evaluation.

## Negative controls validate the pipeline

Height's stromal-dominant enrichment pattern (5.5x in stromal vs
2.2x in immune for Method 1, p <= 10⁻⁸ for the difference) matches the
biology of height-GWAS loci acting through bone/connective-tissue
regulatory elements and provides positive evidence that the pipeline
correctly discriminates tissue of action. Years of education showed
no compartment-specific enrichment above 3x (all p > 0.05 except one
marginal case), as expected for a trait whose genetic architecture
primarily acts in brain tissue not represented in the annotations.
Together these results argue that the IBD enrichments we observe are
biology-driven rather than pipeline artefacts.

## Limitations

Six limitations bound the interpretation of this benchmark.
First, our analyses use European-ancestry reference LD and
European-ancestry GWAS summary statistics; extension to other
ancestries requires ancestry-matched reference panels. Second,
compartment-level pooling lumps developmentally related cell types
together; finer resolution would require either per-cell-type
analysis (feasible for Method 4A, limited for others by the anchor-
agreement issue) or per-state pseudo-bulking. Third, without
cell-type-matched CRISPR validation for colon cells, we cannot
anchor method rankings to an external ground truth; rankings here
are based on heritability partitioning efficiency, which is itself
a meaningful but indirect metric of regulatory validity. Fourth,
we use the power-law Hi-C proxy in both Method 3 (ABC) and Method 5
(rE2G) rather than a real Hi-C cooler, at an estimated 10-15 %
AUPRC cost for enhancer-gene prediction (Nasser et al., 2021);
this decision keeps Methods 3 and 5 comparable but weakens the
absolute claims each can support. Fifth, the GLASSO penalty in
Method 2b uses a single fixed parameter setting (alpha = 0.5,
distance scale = 250 kb) per Pliner et al. (2018); sensitivity to
these choices is not characterised here, and tighter or looser
penalties may shift the 2a-vs-2b comparison. Sixth, the agreement
between Methods 2a and 2b on rank-order against the other linkage
families is observed on a single benchmark (Hickey et al., 2023
colon atlas, three traits); the robustness of this rank-order
agreement across datasets, tissues, ancestries, and traits is not
tested. Method 2a's per-SNP enrichment is up to 4.3x lower than
Method 2b's in matched cells (IBD-immune), a margin large enough
to plausibly flip rankings on noisier datasets, and Method 2a
additionally exhibits a baseline-LD redundancy failure (Height-immune
tau\* = −3.44) that Method 2b does not. Users planning tau\*-based or
baseline-LD-conditioned analyses on other datasets should run both
Methods 2a and 2b and compare, rather than treating either as a
universal substitute for the other.

## Future directions

Three follow-ups are conditional on specific downstream needs and
have been deferred with explicit triggers. If cell-type-specific
claims become necessary, Method 4B should be re-run with either
cell-type-label anchoring or a modern integration method (scVI,
GLUE, or CCA-MNN), each subjected to the same anchor-agreement
diagnostic introduced here. If Method 3 becomes a leading or
trailing method in extended trait lists, a real Hi-C run with the
Nasser et al. (2021) Avg-Hi-C would clarify whether its promoter
bias is primarily a contact-proxy artefact. If Method 5 becomes a
primary result, cell-type-matched CRISPR data - either the
Gschwind et al. (2025) HCT116 screen extended to colonocytes, or a
de novo screen - would enable proper threshold recalibration.
Beyond these, extending the benchmark to other supervised frameworks
(SCENIC+, PolyFun, PoPS) and to additional gut-relevant traits
(colorectal cancer, celiac disease, irritable bowel syndrome) would
test the generalisability of the method rankings reported here.
A small grid scan over the GLASSO penalty parameters (alpha, distance scale)
for Method 2b would also characterise the robustness of the 2a-vs-2b
heritability comparison, and replication on at least one additional
dataset (e.g., a non-gut tissue or a non-IBD trait) would establish
whether the dataset-conditional rank stability we report here holds
beyond Hickey 2023 + de Lange IBD.

---

## References

- Bravo González-Blas, C. *et al.* SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nat. Methods* **20**, 1355-1367 (2023). doi:10.1038/s41592-023-01938-4
- Cao, Z.-J. & Gao, G. Multi-omics single-cell data integration and regulatory inference with graph-linked embedding. *Nat. Biotechnol.* **40**, 1458-1466 (2022). doi:10.1038/s41587-022-01284-4
- de Lange, K. M. *et al.* Genome-wide association study implicates immune activation of multiple integrin genes in inflammatory bowel disease. *Nat. Genet.* **49**, 256-261 (2017). doi:10.1038/ng.3760
- Finucane, H. K. *et al.* Partitioning heritability by functional annotation using genome-wide association summary statistics. *Nat. Genet.* **47**, 1228-1235 (2015). doi:10.1038/ng.3404
- Finucane, H. K. *et al.* Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. *Nat. Genet.* **50**, 621-629 (2018). doi:10.1038/s41588-018-0081-4
- Fulco, C. P. *et al.* Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. *Nat. Genet.* **51**, 1664-1669 (2019). doi:10.1038/s41588-019-0538-0
- Gayoso, A. *et al.* A Python library for probabilistic analysis of single-cell omics data. *Nat. Biotechnol.* **40**, 163-166 (2022). doi:10.1038/s41587-021-01206-w
- Gazal, S. *et al.* Linkage disequilibrium-dependent architecture of human complex traits shows action of negative selection. *Nat. Genet.* **49**, 1421-1427 (2017). doi:10.1038/ng.3954
- Granja, J. M. *et al.* ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nat. Genet.* **53**, 403-411 (2021). doi:10.1038/s41588-021-00790-6
- Gschwind, A. R. *et al.* An encyclopedia of enhancer-gene regulatory interactions in the human genome (ENCODE-rE2G). *Nature* (2025). doi:10.1038/s41586-024-08227-w
- Hickey, J. W. *et al.* Organization of the human intestine at single-cell resolution. *Nature* **619**, 572-584 (2023). doi:10.1038/s41586-023-05915-x
- Nasser, J. *et al.* Genome-wide enhancer maps link risk variants to disease genes. *Nature* **593**, 238-243 (2021). doi:10.1038/s41586-021-03446-x
- Pliner, H. A. *et al.* Cicero predicts cis-regulatory DNA interactions from single-cell chromatin accessibility data. *Mol. Cell* **71**, 858-871.e8 (2018). doi:10.1016/j.molcel.2018.06.044
- Rietveld, C. A. *et al.* GWAS of 126,559 individuals identifies genetic variants associated with educational attainment. *Science* **340**, 1467-1471 (2013). doi:10.1126/science.1235488
- Stuart, T. *et al.* Comprehensive integration of single-cell data. *Cell* **177**, 1888-1902.e21 (2019). doi:10.1016/j.cell.2019.05.031
- Weissbrod, O. *et al.* Leveraging fine-mapping and multipopulation training data to improve cross-population polygenic risk scores. *Nat. Genet.* **54**, 450-458 (2022). doi:10.1038/s41588-022-01036-9
- Yengo, L. *et al.* A saturated map of common genetic variants associated with human height. *Nature* **610**, 704-712 (2022). doi:10.1038/s41586-022-05275-y
