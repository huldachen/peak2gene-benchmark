#!/usr/bin/env python3
"""
Full per-cell-type anchor-agreement table (Claim 2).

The README summarised this as "out of 44 analysed cell types, only 3 clear 25%;
the other 41 have under 10%" -- a count, which invites the objection that the
result is an artifact of where the threshold sits. Tabulating the whole
distribution answers that, and it also surfaces two things the summary hid:

  * the committed per-cell-type data contains 47 cell types, not 44
  * 10 of them have n < 100 cells, the smallest n = 2, where a mean agreement
    is not an estimate of anything

So the table reports n beside every value, and the summary statistics are given
both over all cell types and over those with n >= 100, so a reader can see the
conclusion does not depend on the tiny strata.

Input : results/linkage/method4b_anchor_diagnostic/*_per_cell_type.tsv
Output: results/tables/table3_anchor_agreement_full.tsv + stdout summary
"""
from __future__ import annotations

import glob
import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
SRC = os.path.join(ROOT, "results", "linkage", "method4b_anchor_diagnostic")
OUT = os.path.join(ROOT, "results", "tables", "table3_anchor_agreement_full.tsv")
MIN_N = 100          # below this a per-cell-type mean is not interpretable
CLEAR = 0.25
LOW = 0.10


def main() -> int:
    frames = []
    for f in sorted(glob.glob(os.path.join(SRC, "*_per_cell_type.tsv"))):
        comp = os.path.basename(f).split("_per_cell_type")[0]
        d = pd.read_csv(f, sep="\t")
        d.insert(0, "compartment", comp)
        frames.append(d)
    a = pd.concat(frames, ignore_index=True)
    a = a.rename(columns={"cell_type_atac": "cell_type", "n": "n_cells"})
    a["interpretable"] = a.n_cells >= MIN_N
    a = a.sort_values(["compartment", "mean_agreement"], ascending=[True, False])
    a.to_csv(OUT, sep="\t", index=False,
             float_format="%.4f", columns=["compartment", "cell_type", "n_cells",
                                           "mean_agreement", "median_agreement",
                                           "interpretable"])

    def summarise(df, label):
        n = len(df)
        hi = int((df.mean_agreement >= CLEAR).sum())
        lo = int((df.mean_agreement < LOW).sum())
        mid = n - hi - lo
        print(f"  {label:<28} n={n:>3}   >=25%: {hi}   10-25%: {mid}   <10%: {lo}")

    print(f"Anchor agreement, per cell type (n >= {MIN_N} = interpretable)\n")
    summarise(a, "all cell types")
    summarise(a[a.interpretable], f"n >= {MIN_N} only")
    print(f"\n  cell types with n < {MIN_N}: {int((~a.interpretable).sum())} "
          f"(smallest n = {int(a.n_cells.min())})")
    print(f"\n  clearing {CLEAR:.0%}:")
    for _, r in a[a.mean_agreement >= CLEAR].iterrows():
        print(f"    {r.compartment:<11} {r.cell_type:<28} n={int(r.n_cells):>5}  "
              f"mean={r.mean_agreement:.3f}")
    print(f"\nSaved -> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
