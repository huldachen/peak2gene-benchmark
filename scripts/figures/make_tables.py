"""
Generate the 4 tables for the writeup in TSV, Markdown, and LaTeX formats.

Tables
------
  1. sLDSC benchmark results (main finding, one row per method x compartment)
  2. Method design reference (inputs, algorithm, key parameters, citations)
  3. Cross-cohort anchor-agreement diagnostic (Finding E quantitative)

Outputs written to results/tables/:
  table1_ldsc_benchmark.{tsv,md,tex}
  table2_method_design.{tsv,md,tex}
  table3_anchor_agreement.{tsv,md,tex}
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
TBL_DIR = ROOT / "results/tables"
TBL_DIR.mkdir(parents=True, exist_ok=True)

ENRICH_TSV = ROOT / "results/ldsc/enrichment_table.tsv"
ANCHOR_DIR = ROOT / "results/linkage/method4b_anchor_diagnostic"


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _write_all(df: pd.DataFrame, stem: str, *,
               md_kwargs=None, tex_kwargs=None, tsv_kwargs=None) -> None:
    """Write one DataFrame to TSV, Markdown, and LaTeX simultaneously."""
    md_kwargs = md_kwargs or {}
    tex_kwargs = tex_kwargs or {}
    tsv_kwargs = tsv_kwargs or {}

    tsv_path = TBL_DIR / f"{stem}.tsv"
    md_path  = TBL_DIR / f"{stem}.md"
    tex_path = TBL_DIR / f"{stem}.tex"

    df.to_csv(tsv_path, sep="\t", index=False, **tsv_kwargs)

    log.info("Wrote %s  (%d rows x %d cols) -> %s", stem, len(df), df.shape[1], tsv_path)


def _fmt_enrich(e, se, p):
    """'X.X [+/-Y.Y]***' - enrichment, SE, stars for p."""
    if pd.isna(e):
        return "-"
    stars = ""
    if p is not None and not pd.isna(p):
        if   p < 0.001: stars = "***"
        elif p < 0.01:  stars = "**"
        elif p < 0.05:  stars = "*"
    # Handle negative enrichments cleanly (unusual but present in 4B stromal)
    return f"{e:.2f} +/- {se:.2f}{stars}"


def _fmt_p(p):
    if pd.isna(p):
        return "-"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"


# ═══════════════════════════════════════════════════════════════════════════
# Table 1 - sLDSC benchmark (headline result)
# ═══════════════════════════════════════════════════════════════════════════

def make_table1() -> None:
    """Compact benchmark table - compressed row labels, stars-only format."""
    df = pd.read_csv(ENRICH_TSV, sep="\t")

    # Human-readable method labels - collapse 5-binary and 5-continuous
    # (they are numerically indistinguishable) into one combined row
    method_rename = {
        "method1_distance":         "1. Distance window",
        "method2_cicero":           "2a. Co-access (Pearson)",
        "method2b_glasso":          "2b. Co-access (GLASSO)",
        "method3_abc":              "3. ABC (power-law)",
        "method4a_paired":          "4A. Paired Multiome r",
        "method4b_crosscohort":     "4B. Cross-cohort r",
        "method5_re2g":             "5. rE2G (K562; binary ≈ continuous)",
        "method5_re2g_continuous":  None,   # skip - duplicate of binary
    }
    df["method_label"] = df["method"].map(method_rename)
    df = df.dropna(subset=["method_label"]).copy()

    # Ordered list of unique method labels (drop None)
    method_order = [m for m in method_rename.values() if m is not None]
    comp_order = ["epithelial", "immune", "stromal"]
    trait_order = ["IBD", "Height", "EA"]

    # Short-format enrichment: "X.X*" (stars indicate significance)
    def _short_enrich(e, p):
        if pd.isna(e):
            return "-"
        stars = ""
        if p is not None and not pd.isna(p):
            if   p < 0.001: stars = "***"
            elif p < 0.01:  stars = "**"
            elif p < 0.05:  stars = "*"
        return f"{e:.1f}{stars}"

    # Build compact wide table: one row per (method, compartment), short cells
    rows = []
    for method in method_order:
        for comp in comp_order:
            row = {"Method": method, "Compartment": comp}
            for trait in trait_order:
                sub = df[(df["method_label"] == method)
                         & (df["compartment"] == comp)
                         & (df["trait"] == trait)]
                if sub.empty:
                    row[trait] = "-"
                else:
                    r = sub.iloc[0]
                    row[trait] = _short_enrich(r["enrichment"], r["enrichment_p"])
            rows.append(row)

    clean = pd.DataFrame(rows)
    clean = clean.rename(columns={
        "IBD":    "IBD",
        "Height": "Height",
        "EA":     "EA",
    })

    _write_all(clean, "table1_ldsc_benchmark")

    # Also keep the extended version (with SE + p-values) as supplementary
    full_rows = []
    for method in method_order:
        for comp in comp_order:
            row = {"Method": method, "Compartment": comp}
            for trait in trait_order:
                sub = df[(df["method_label"] == method)
                         & (df["compartment"] == comp)
                         & (df["trait"] == trait)]
                if sub.empty:
                    row[trait] = row[f"{trait} p"] = "-"
                else:
                    r = sub.iloc[0]
                    row[trait] = _fmt_enrich(r["enrichment"], r["enrichment_se"], r["enrichment_p"])
                    row[f"{trait} p"] = _fmt_p(r["enrichment_p"])
            full_rows.append(row)
    extended = pd.DataFrame(full_rows)[
        ["Method", "Compartment", "IBD", "IBD p", "Height", "Height p", "EA", "EA p"]
    ]
    _write_all(extended, "table1_ldsc_benchmark_extended")


# ═══════════════════════════════════════════════════════════════════════════
# Table 2 - Method design reference
# ═══════════════════════════════════════════════════════════════════════════

def make_table2() -> None:
    """Methods reference - everything a reader needs to know about each method."""
    rows = [
        {
            "Method":            "1. Distance window",
            "Inputs":            "ATAC peak BED + GTF",
            "Algorithm":         "Peaks within +/-500 kb of any protein-coding TSS",
            "Key parameter":     "window = +/-500 kb",
            "Threshold":         "none (baseline)",
            "Canonical source":  "GTEx, ABC paper use +/-5 Mb",
            "Our setting":       "+/-500 kb (consistent with Methods 2/4)",
            "Defensibility":     "✓ standard",
        },
        {
            "Method":            "2a. Co-access (Pearson)",
            "Inputs":            "ATAC binary matrix",
            "Algorithm":         "LSI -> KNN meta-cells -> per-window Pearson r between peaks",
            "Key parameter":     "k=50, 2000 meta-cells, r >= 0.25, promoter +/-2 kb",
            "Threshold":         "r >= 0.25 (raw, no graphical lasso)",
            "Canonical source":  "ArchR getCoAccessibility (Granja 2021)",
            "Our setting":       "r >= 0.25 - pragmatic (no graphical lasso)",
            "Defensibility":     "⚠ includes indirect correlations",
        },
        {
            "Method":            "2b. Co-access (GLASSO)",
            "Inputs":            "ATAC binary matrix",
            "Algorithm":         "LSI -> KNN meta-cells -> per-window graphical-lasso partial corr",
            "Key parameter":     "k=50, 2000 meta-cells, alpha=0.5, pcor >= 0.05, promoter +/-2 kb",
            "Threshold":         "pcor >= 0.05 (Cicero shrunk-score scale)",
            "Canonical source":  "Pliner 2018 Mol Cell (Cicero with distance-penalised GLASSO)",
            "Our setting":       "pcor >= 0.05 with distance-dependent penalty",
            "Defensibility":     "✓ faithful to Cicero; removes indirect edges",
        },
        {
            "Method":            "3. ABC (power-law)",
            "Inputs":            "ATAC pseudobulk + distance",
            "Algorithm":         "Activity x Contact, per-gene normalised",
            "Key parameter":     "window = +/-5 Mb, γ = 0.87, ABC >= 0.02",
            "Threshold":         "ABC >= 0.02 (Fulco 2019)",
            "Canonical source":  "Fulco 2019 Nat Genet, Nasser 2021 Nature",
            "Our setting":       "power-law (no H3K27ac, no real Hi-C)",
            "Defensibility":     "✓ per Nasser 2021 fallback guidance (∼10% AUPRC loss)",
        },
        {
            "Method":            "4A. Paired Multiome r",
            "Inputs":            "paired RNA + ATAC (same cells, Multiome cohort)",
            "Algorithm":         "Meta-cell Pearson r between peak access and gene expr.",
            "Key parameter":     "+/-500 kb, 2000 meta-cells, |r| >= 0.45, FDR <= 0.05",
            "Threshold":         "|r| >= 0.45 + FDR (per-gene BH)",
            "Canonical source":  "ArchR addPeak2GeneLinks defaults",
            "Our setting":       "ArchR defaults, 100% paired cells",
            "Defensibility":     "✓ gold standard (true paired)",
        },
        {
            "Method":            "4B. Cross-cohort r",
            "Inputs":            "RNA cohort (ref) + Multiome ATAC (target)",
            "Algorithm":         "gene-activity KNN -> anchor transfer -> per-gene Pearson r",
            "Key parameter":     "k=30 anchors, |r| >= 0.45, FDR <= 0.05",
            "Threshold":         "|r| >= 0.45 + FDR",
            "Canonical source":  "Seurat FindTransferAnchors / ArchR addGeneIntegrationMatrix",
            "Our setting":       "simplified: gene-activity cosine KNN",
            "Defensibility":     "⚠ Finding E: cell-type anchors ≈ 0% for most types",
        },
        {
            "Method":            "5. rE2G (supervised)",
            "Inputs":            "8 ABC-style features + pretrained K562 LR model",
            "Algorithm":         "log-transformed features -> logistic regression",
            "Key parameter":     "threshold = 0.179 (K562 70% recall calibration)",
            "Threshold":         "0.179 binary; continuous variant tested",
            "Canonical source":  "Gschwind 2025 Nature (EngreitzLab/ENCODE_rE2G)",
            "Our setting":       "atac_megamap model, power-law Hi-C substitute",
            "Defensibility":     "⚠ K562 threshold not calibrated for gut cell types",
        },
    ]
    df = pd.DataFrame(rows)
    _write_all(df, "table2_method_design")


# ═══════════════════════════════════════════════════════════════════════════
# Table 3 - Anchor-agreement diagnostic (Finding E)
# ═══════════════════════════════════════════════════════════════════════════

def make_table3() -> None:
    """Per-compartment breakdown of anchor-agreement per cell type."""
    rows_all = []
    for comp in ("epithelial", "immune", "stromal"):
        f = ANCHOR_DIR / f"{comp}_per_cell_type.tsv"
        if not f.exists():
            log.warning("Missing %s", f)
            continue
        df = pd.read_csv(f, sep="\t")
        df["compartment"] = comp
        rows_all.append(df)
    if not rows_all:
        log.error("No anchor-agreement files found")
        return
    df = pd.concat(rows_all, ignore_index=True)

    # Keep the 3 most abundant cell types per compartment + any with agreement > 0.25
    keep_rows = []
    for comp, grp in df.groupby("compartment", sort=False):
        top3 = grp.nlargest(3, "n")
        winners = grp[grp["mean_agreement"] >= 0.25]  # "winners" with >= 25% anchor agreement
        combined = pd.concat([top3, winners]).drop_duplicates(subset=["cell_type_atac"])
        combined = combined.sort_values("mean_agreement", ascending=False)
        keep_rows.append(combined)
    out = pd.concat(keep_rows, ignore_index=True)
    out = out[["compartment", "cell_type_atac", "n", "mean_agreement", "median_agreement"]]
    out["mean_agreement"]   = out["mean_agreement"].apply(lambda x: f"{x:.3f}")
    out["median_agreement"] = out["median_agreement"].apply(lambda x: f"{x:.3f}")
    out = out.rename(columns={
        "compartment":       "Compartment",
        "cell_type_atac":    "Cell type",
        "n":                 "n cells",
        "mean_agreement":    "Mean agreement",
        "median_agreement":  "Median agreement",
    })

    # Also add a summary row per compartment (overall mean across ALL cells)
    all_rows = []
    for comp, grp in df.groupby("compartment", sort=False):
        n_all = int(grp["n"].sum())
        weighted_mean = float((grp["mean_agreement"] * grp["n"]).sum() / max(n_all, 1))
        all_rows.append({
            "Compartment":       comp,
            "Cell type":         f"[ ALL - {grp.shape[0]} cell types ]",
            "n cells":           n_all,
            "Mean agreement":    f"{weighted_mean:.3f}",
            "Median agreement":  "(see row-level)",
        })
    summary = pd.DataFrame(all_rows)

    # Interleave: one summary row followed by that compartment's top rows
    interleaved = []
    for comp in ("epithelial", "immune", "stromal"):
        summary_row = summary[summary["Compartment"] == comp]
        detail_rows = out[out["Compartment"] == comp]
        interleaved.append(pd.concat([summary_row, detail_rows], ignore_index=True))
    final = pd.concat(interleaved, ignore_index=True)

    _write_all(final, "table3_anchor_agreement")


# ═══════════════════════════════════════════════════════════════════════════
# Table 4 - Total heritability captured (prop_h2) - supplementary
# ═══════════════════════════════════════════════════════════════════════════

def make_table4() -> None:
    """Total h^2 captured per (method x compartment x trait), with annotation
    size - the 'portfolio coverage' view that complements per-SNP enrichment
    and tau\\*. Lets the reader see: a small annotation may be highly enriched
    per SNP yet capture little total heritability; a larger annotation may be
    less concentrated but cover more of the signal."""
    df = pd.read_csv(ENRICH_TSV, sep="\t")

    method_rename = {
        "method1_distance":         "1. Distance window",
        "method2_cicero":           "2a. Co-access (Pearson)",
        "method2b_glasso":          "2b. Co-access (GLASSO)",
        "method3_abc":              "3. ABC (power-law)",
        "method4a_paired":          "4A. Paired Multiome r",
        "method4b_crosscohort":     "4B. Cross-cohort r",
        "method5_re2g":             "5. rE2G (K562)",
        "method5_re2g_continuous":  None,
    }
    df["method_label"] = df["method"].map(method_rename)
    df = df.dropna(subset=["method_label"]).copy()

    method_order = [m for m in method_rename.values() if m is not None]
    comp_order = ["epithelial", "immune", "stromal"]
    trait_order = ["IBD", "Height", "EA"]

    def _fmt_pct(x):
        if pd.isna(x):
            return "-"
        return f"{100*x:.2f}%"

    def _fmt_h2(p, se):
        if pd.isna(p):
            return "-"
        return f"{100*p:.1f} +/- {100*se:.1f}%"

    rows = []
    for method in method_order:
        for comp in comp_order:
            row = {"Method": method, "Compartment": comp}
            # prop_snps is the same across traits within a (method, comp) cell
            sub_any = df[(df["method_label"] == method) & (df["compartment"] == comp)]
            if sub_any.empty:
                row["prop_snps"] = "-"
            else:
                row["prop_snps"] = _fmt_pct(sub_any.iloc[0]["prop_snps"])
            for trait in trait_order:
                sub = sub_any[sub_any["trait"] == trait]
                if sub.empty:
                    row[f"{trait} prop_h2"] = "-"
                else:
                    r = sub.iloc[0]
                    row[f"{trait} prop_h2"] = _fmt_h2(r["prop_h2"], r["prop_h2_se"])
            rows.append(row)

    out = pd.DataFrame(rows)[
        ["Method", "Compartment", "prop_snps",
         "IBD prop_h2", "Height prop_h2", "EA prop_h2"]
    ]
    _write_all(out, "table4_prop_h2_supplement")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    make_table1()
    make_table2()
    make_table3()
    make_table4()
    log.info("All tables written to %s", TBL_DIR)


if __name__ == "__main__":
    main()
