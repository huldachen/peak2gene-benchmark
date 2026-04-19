"""
Generate the writeup tables as TSV files.

Produces three tables under results/tables/: the sLDSC benchmark headline
table, the method-design reference, and the cross-cohort anchor-agreement
diagnostic.

Usage
-----
    python scripts/figures/make_tables.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

ROOT    = Path(__file__).resolve().parents[2]
TBL_DIR = ROOT / "results/tables"
TBL_DIR.mkdir(parents=True, exist_ok=True)

ENRICH_TSV = ROOT / "results/ldsc/enrichment_table.tsv"
ANCHOR_DIR = ROOT / "results/linkage/method4b_anchor_diagnostic"


def _write_all(df: pd.DataFrame, stem: str) -> None:
    """Write one DataFrame to a TSV file."""
    tsv_path = TBL_DIR / f"{stem}.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    log.info("Wrote %s  (%d rows x %d cols)  %s", stem, len(df), df.shape[1], tsv_path)


def _fmt_enrich(e, se, p):
    """Format enrichment +/- SE with stars for p."""
    if pd.isna(e):
        return "—"
    stars = ""
    if p is not None and not pd.isna(p):
        if   p < 0.001: stars = "***"
        elif p < 0.01:  stars = "**"
        elif p < 0.05:  stars = "*"
    return f"{e:.2f} ± {se:.2f}{stars}"


def _fmt_p(p):
    if pd.isna(p):
        return "—"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"


def make_table1() -> None:
    """Compact benchmark table."""
    df = pd.read_csv(ENRICH_TSV, sep="\t")

    # method5 binary and continuous are numerically indistinguishable; collapse to one row
    method_rename = {
        "method1_distance":         "1. Distance window",
        "method2_cicero":           "2. Cicero co-access",
        "method3_abc":              "3. ABC (power-law)",
        "method4a_paired":          "4A. Paired Multiome r",
        "method4b_crosscohort":     "4B. Cross-cohort r",
        "method5_re2g":             "5. rE2G (K562; binary ≈ continuous)",
        "method5_re2g_continuous":  None,
    }
    df["method_label"] = df["method"].map(method_rename)
    df = df.dropna(subset=["method_label"]).copy()

    method_order = [m for m in method_rename.values() if m is not None]
    comp_order = ["epithelial", "immune", "stromal"]
    trait_order = ["IBD", "Height", "EA"]

    def _short_enrich(e, p):
        if pd.isna(e):
            return "—"
        stars = ""
        if p is not None and not pd.isna(p):
            if   p < 0.001: stars = "***"
            elif p < 0.01:  stars = "**"
            elif p < 0.05:  stars = "*"
        return f"{e:.1f}{stars}"

    rows = []
    for method in method_order:
        for comp in comp_order:
            row = {"Method": method, "Compartment": comp}
            for trait in trait_order:
                sub = df[(df["method_label"] == method)
                         & (df["compartment"] == comp)
                         & (df["trait"] == trait)]
                if sub.empty:
                    row[trait] = "—"
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

    full_rows = []
    for method in method_order:
        for comp in comp_order:
            row = {"Method": method, "Compartment": comp}
            for trait in trait_order:
                sub = df[(df["method_label"] == method)
                         & (df["compartment"] == comp)
                         & (df["trait"] == trait)]
                if sub.empty:
                    row[trait] = row[f"{trait} p"] = "—"
                else:
                    r = sub.iloc[0]
                    row[trait] = _fmt_enrich(r["enrichment"], r["enrichment_se"], r["enrichment_p"])
                    row[f"{trait} p"] = _fmt_p(r["enrichment_p"])
            full_rows.append(row)
    extended = pd.DataFrame(full_rows)[
        ["Method", "Compartment", "IBD", "IBD p", "Height", "Height p", "EA", "EA p"]
    ]
    _write_all(extended, "table1_ldsc_benchmark_extended")


def make_table2() -> None:
    """Methods reference — one row per method."""
    rows = [
        {
            "Method":            "1. Distance window",
            "Inputs":            "ATAC peak BED + GTF",
            "Algorithm":         "Peaks within ±500 kb of any protein-coding TSS",
            "Key parameter":     "window = ±500 kb",
            "Threshold":         "none (baseline)",
            "Canonical source":  "GTEx, ABC paper use ±5 Mb",
            "Our setting":       "±500 kb (consistent with Methods 2/4)",
            "Defensibility":     "✓ standard",
        },
        {
            "Method":            "2. Cicero co-access",
            "Inputs":            "ATAC binary matrix",
            "Algorithm":         "LSI → KNN meta-cells → per-window Pearson r between peaks",
            "Key parameter":     "k=50, 2000 meta-cells, r ≥ 0.25, promoter ±2 kb",
            "Threshold":         "r ≥ 0.25 (raw, no graphical lasso)",
            "Canonical source":  "Pliner 2018 Mol Cell (r ≥ 0.05 on shrunk scores)",
            "Our setting":       "r ≥ 0.25 — pragmatic (no graphical lasso)",
            "Defensibility":     "⚠ threshold not calibrated against CRISPR truth",
        },
        {
            "Method":            "3. ABC (power-law)",
            "Inputs":            "ATAC pseudobulk + distance",
            "Algorithm":         "Activity × Contact, per-gene normalised",
            "Key parameter":     "window = ±5 Mb, γ = 0.87, ABC ≥ 0.02",
            "Threshold":         "ABC ≥ 0.02 (Fulco 2019)",
            "Canonical source":  "Fulco 2019 Nat Genet, Nasser 2021 Nature",
            "Our setting":       "power-law (no H3K27ac, no real Hi-C)",
            "Defensibility":     "✓ per Nasser 2021 fallback guidance (∼10% AUPRC loss)",
        },
        {
            "Method":            "4A. Paired Multiome r",
            "Inputs":            "paired RNA + ATAC (same cells, Multiome cohort)",
            "Algorithm":         "Meta-cell Pearson r between peak access and gene expr.",
            "Key parameter":     "±500 kb, 2000 meta-cells, |r| ≥ 0.45, FDR ≤ 0.05",
            "Threshold":         "|r| ≥ 0.45 + FDR (per-gene BH)",
            "Canonical source":  "ArchR addPeak2GeneLinks defaults",
            "Our setting":       "ArchR defaults, 100% paired cells",
            "Defensibility":     "✓ gold standard (true paired)",
        },
        {
            "Method":            "4B. Cross-cohort r",
            "Inputs":            "RNA cohort (ref) + Multiome ATAC (target)",
            "Algorithm":         "gene-activity KNN → anchor transfer → per-gene Pearson r",
            "Key parameter":     "k=30 anchors, |r| ≥ 0.45, FDR ≤ 0.05",
            "Threshold":         "|r| ≥ 0.45 + FDR",
            "Canonical source":  "Seurat FindTransferAnchors / ArchR addGeneIntegrationMatrix",
            "Our setting":       "simplified: gene-activity cosine KNN",
            "Defensibility":     "⚠ Finding E: cell-type anchors ≈ 0% for most types",
        },
        {
            "Method":            "5. rE2G (supervised)",
            "Inputs":            "8 ABC-style features + pretrained K562 LR model",
            "Algorithm":         "log-transformed features → logistic regression",
            "Key parameter":     "threshold = 0.179 (K562 70% recall calibration)",
            "Threshold":         "0.179 binary; continuous variant tested",
            "Canonical source":  "Gschwind 2025 Nature (EngreitzLab/ENCODE_rE2G)",
            "Our setting":       "atac_megamap model, power-law Hi-C substitute",
            "Defensibility":     "⚠ K562 threshold not calibrated for gut cell types",
        },
    ]
    df = pd.DataFrame(rows)
    _write_all(df, "table2_method_design")


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

    keep_rows = []
    for comp, grp in df.groupby("compartment", sort=False):
        top3 = grp.nlargest(3, "n")
        winners = grp[grp["mean_agreement"] >= 0.25]
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

    all_rows = []
    for comp, grp in df.groupby("compartment", sort=False):
        n_all = int(grp["n"].sum())
        weighted_mean = float((grp["mean_agreement"] * grp["n"]).sum() / max(n_all, 1))
        all_rows.append({
            "Compartment":       comp,
            "Cell type":         f"[ ALL — {grp.shape[0]} cell types ]",
            "n cells":           n_all,
            "Mean agreement":    f"{weighted_mean:.3f}",
            "Median agreement":  "(see row-level)",
        })
    summary = pd.DataFrame(all_rows)

    interleaved = []
    for comp in ("epithelial", "immune", "stromal"):
        summary_row = summary[summary["Compartment"] == comp]
        detail_rows = out[out["Compartment"] == comp]
        interleaved.append(pd.concat([summary_row, detail_rows], ignore_index=True))
    final = pd.concat(interleaved, ignore_index=True)

    _write_all(final, "table3_anchor_agreement")


def main():
    make_table1()
    make_table2()
    make_table3()
    log.info("All tables written to %s", TBL_DIR)


if __name__ == "__main__":
    main()
