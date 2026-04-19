# Peak–gene linkage outputs

Each method writes one TSV per compartment (epithelial / immune / stromal) with the unified schema:

```
peak_id  gene_id  gene_name  chrom  peak_start  peak_end  tss_pos
strand   distance_bp  compartment  method  score
```

Methods 4A and 4B additionally emit `r`, `pvalue`, and `fdr`; Method 5 adds `rE2G_probability` and `rE2G_binary`.

## What's in this repository

Only the first **1,000 rows** of each full output is committed, as `{method}/{compartment}_sample.tsv`. Reviewers can inspect the schema, peek at real rows, and check column types without downloading gigabytes.

**Why not the full files?** Full linkage tables reach 764 MB (Method 1 epithelial), well above GitHub's 100 MB per-file limit, for a cumulative ~3 GB that is trivially regenerated from the scripts.

## Regenerating the full outputs

After cloning and setting up the `omics` conda environment (see root `README.md`):

```bash
python scripts/linkage/method1_distance.py       --all-compartments --gtf … --chrom-sizes … --out-dir results/linkage/method1/
python scripts/linkage/method2_cicero.py         --all-compartments --gtf …                --out-dir results/linkage/method2/
python scripts/linkage/method3_abc.py            --all-compartments --gtf …                --out-dir results/linkage/method3/
python scripts/linkage/method4_paired.py         --all-compartments --rna-h5ad … --gtf …   --out-dir results/linkage/method4_paired/
python scripts/linkage/method4_crosscohort.py    --all-compartments --rna-h5ad … --gtf …   --out-dir results/linkage/method4_crosscohort/
python scripts/linkage/method5_re2g.py           --all-compartments --gtf …                --out-dir results/linkage/method5/
```

Total wall-clock for the full regeneration of all 18 outputs is ~50–60 minutes; see per-script times in the root README.

## Full-size reference

| Method | compartment | full rows | full size |
|--------|-------------|----------:|----------:|
| method1_distance       | epithelial | 6,669,636 | 764 MB |
| method1_distance       | immune     | 3,518,511 | 390 MB |
| method1_distance       | stromal    | 4,511,169 | 504 MB |
| method2_cicero         | epithelial | 1,201,659 | 177 MB |
| method2_cicero         | immune     |   599,017 |  86 MB |
| method2_cicero         | stromal    |   640,553 |  93 MB |
| method3_abc            | epithelial |    93,989 |  12 MB |
| method3_abc            | immune     |   111,164 |  13 MB |
| method3_abc            | stromal    |    97,669 |  12 MB |
| method4_paired         | epithelial |    24,431 | 3.3 MB |
| method4_paired         | immune     |    48,332 | 6.4 MB |
| method4_paired         | stromal    |    24,353 | 3.2 MB |
| method4_crosscohort    | epithelial |    37,200 | 5.0 MB |
| method4_crosscohort    | immune     |    71,003 | 9.1 MB |
| method4_crosscohort    | stromal    |    36,971 | 4.8 MB |
| method5_re2g           | epithelial | 3,313,375 | 481 MB |
| method5_re2g           | immune     | 1,611,080 | 227 MB |
| method5_re2g           | stromal    | 2,214,745 | 315 MB |

Method 4B's anchor-agreement diagnostic (`method4b_anchor_diagnostic/`) is small enough to ship in full — those files back Figure 4 and Table 3 directly.
