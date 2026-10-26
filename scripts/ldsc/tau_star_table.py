#!/usr/bin/env python3
"""
tau* as a reader-facing table (the size-independent ranking).

Enrichment = prop_h2 / prop_snps mechanically favours compact annotations, so a
ranking by enrichment alone invites the objection that paired multiome wins on
compactness rather than biology. tau* -- the coefficient's contribution beyond
the baseline-LD model -- is invariant to annotation size and is the metric that
answers it.

The repo already computed tau* and discusses it, but it appeared only in
Figure 3 Row B and in prose; no committed table surfaced it, so a reader
checking the ranking saw only the size-sensitive metric. This emits it.

tau* is reported as a z-score (coefficient / SE), which already carries its own
uncertainty: |z| > 1.96 is the 95% test, so no separate interval is needed.

Input : results/ldsc/enrichment_table.tsv
Output: results/tables/table6_tau_star.tsv + stdout
"""
from __future__ import annotations

import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
SRC = os.path.join(ROOT, "results", "ldsc", "enrichment_table.tsv")
OUT = os.path.join(ROOT, "results", "tables", "table6_tau_star.tsv")
LABEL = {
    "method1_distance": "1. Distance window",
    "method2_cicero": "2a. Co-access (Pearson)",
    "method2b_glasso": "2b. Co-access (GLASSO)",
    "method3_abc": "3. ABC (power-law)",
    "method4a_paired": "4A. Paired Multiome r",
    "method4b_crosscohort": "4B. Cross-cohort r",
    "method5_re2g": "5. rE2G (K562)",
    "method5_re2g_continuous": "5c. rE2G (continuous)",
}
Z = 1.96


def main() -> int:
    d = pd.read_csv(SRC, sep="\t")
    d = d[d.method.isin(LABEL)].copy()
    d["Method"] = d.method.map(LABEL)
    w = d.pivot_table(index=["Method", "compartment"], columns="trait",
                      values="coefficient_z").reset_index()
    cov = (d[d.trait == "IBD"].set_index(["Method", "compartment"])["prop_snps"] * 100)
    w["prop_snps_%"] = [round(cov.get((m, c), float("nan")), 2)
                        for m, c in zip(w.Method, w.compartment)]
    w = w[["Method", "compartment", "prop_snps_%", "IBD", "Height", "EA"]]
    w = w.sort_values(["Method", "compartment"])
    w.to_csv(OUT, sep="\t", index=False, float_format="%.2f")

    print("tau* (coefficient z-score) — size-independent; |z| > 1.96 is 95%\n")
    print(f"{'Method':<24}{'compartment':<12}{'cover%':>7}{'IBD':>8}{'Height':>8}{'EA':>8}")
    print("-" * 67)
    for _, r in w.iterrows():
        star = lambda v: f"{v:+.2f}" + ("*" if abs(v) > Z else " ")
        print(f"{r.Method:<24}{r.compartment:<12}{r['prop_snps_%']:>7.2f}"
              f"{star(r.IBD):>9}{star(r.Height):>9}{star(r.EA):>9}")
    sig = ((w[["IBD", "Height", "EA"]].abs() > Z).sum().sum())
    print(f"\n  cells with |tau* z| > 1.96: {int(sig)} / {w.shape[0]*3}")
    print(f"  annotation coverage spans {w['prop_snps_%'].min():.2f}%–"
          f"{w['prop_snps_%'].max():.2f}% of SNPs ({w['prop_snps_%'].max()/w['prop_snps_%'].min():.0f}x)")
    print(f"\nSaved -> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
