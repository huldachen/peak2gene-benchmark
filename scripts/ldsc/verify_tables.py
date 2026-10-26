#!/usr/bin/env python3
"""
Consistency guard: do the published tables still match the raw sLDSC output?

Written after a prose claim in results/DISCUSSION.md ("IBD enrichment moves
upward across all three compartments") was found to contradict the committed
numbers -- it moves DOWN in two of three. The tables and the raw output were
both correct; the narrative had drifted. Numbers quoted in prose are the layer
most likely to rot, because nothing recomputes them.

This script cannot re-run sLDSC (that is ~6 h of LD scores). What it can do is
verify every link in the chain that is checkable in seconds:

  1. enrichment == prop_h2 / prop_snps              (definition of enrichment)
  2. coefficient_z == coefficient / coefficient_se  (definition of tau* z)
  3. table1_ldsc_benchmark_extended.tsv enrichments == raw
  4. table4_prop_h2_supplement.tsv prop_snps/prop_h2 == raw

A failure means a table was hand-edited or a pipeline step changed without the
tables being regenerated. Exits nonzero so it can gate a commit.

Run:  python scripts/ldsc/verify_tables.py
"""

from __future__ import annotations

import csv
import os
import re
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
RAW = os.path.join(PROJECT_ROOT, "results", "ldsc", "enrichment_table.tsv")
T1 = os.path.join(PROJECT_ROOT, "results", "tables", "table1_ldsc_benchmark_extended.tsv")
T4 = os.path.join(PROJECT_ROOT, "results", "tables", "table4_prop_h2_supplement.tsv")

METHOD_LABELS = {
    "1. Distance window": "method1_distance",
    "2a. Co-access (Pearson)": "method2_cicero",
    "2b. Co-access (GLASSO)": "method2b_glasso",
    "3. ABC (power-law)": "method3_abc",
    "4A. Paired Multiome r": "method4a_paired",
    "4B. Cross-cohort r": "method4b_crosscohort",
    "5. rE2G (K562; binary ≈ continuous)": "method5_re2g",
    "5. rE2G (K562)": "method5_re2g",
}
TRAITS = ("IBD", "Height", "EA")
REL_TOL = 1e-6      # invariants should hold to machine precision
ABS_TOL_ENR = 0.02  # published tables are rounded to 2 dp
ABS_TOL_PCT = 0.1


def _lead_number(s: str) -> float | None:
    m = re.match(r"\s*(-?[\d.]+)", s or "")
    return float(m.group(1)) if m else None


def main() -> int:
    raw_rows = list(csv.DictReader(open(RAW), delimiter="\t"))
    raw = {(r["trait"], r["compartment"], r["method"]): r for r in raw_rows}
    failures: list[str] = []

    # --- 1 & 2: definitional invariants on the raw output ---
    for r in raw_rows:
        tag = f"{r['trait']}/{r['method']}/{r['compartment']}"
        ps, ph, e = float(r["prop_snps"]), float(r["prop_h2"]), float(r["enrichment"])
        if ps and abs(ph / ps - e) > REL_TOL * max(abs(e), 1.0):
            failures.append(f"[inv1] {tag}: enrichment {e} != prop_h2/prop_snps {ph / ps}")
        c, cse, cz = float(r["coefficient"]), float(r["coefficient_se"]), float(r["coefficient_z"])
        if cse and abs(c / cse - cz) > REL_TOL * max(abs(cz), 1.0):
            failures.append(f"[inv2] {tag}: coefficient_z {cz} != coef/se {c / cse}")

    # --- 3: table1 enrichments match raw ---
    n1 = 0
    for row in csv.DictReader(open(T1), delimiter="\t"):
        key = METHOD_LABELS.get(row["Method"].strip())
        if key is None:
            failures.append(f"[t1] unmapped method label: {row['Method']!r}")
            continue
        for trait in TRAITS:
            rr = raw.get((trait, row["Compartment"].strip(), key))
            if rr is None:
                continue
            n1 += 1
            got, exp = _lead_number(row[trait]), float(rr["enrichment"])
            if got is None or abs(got - exp) > ABS_TOL_ENR:
                failures.append(f"[t1] {trait}/{row['Method']}/{row['Compartment']}: "
                                f"table {got} != raw {exp:.2f}")

    # --- 4: table4 coverage/h2 match raw ---
    n4 = 0
    for row in csv.DictReader(open(T4), delimiter="\t"):
        key = METHOD_LABELS.get(row["Method"].strip())
        if key is None:
            failures.append(f"[t4] unmapped method label: {row['Method']!r}")
            continue
        rr = raw.get(("IBD", row["Compartment"].strip(), key))
        if rr is None:
            continue
        n4 += 1
        ps, exp_ps = _lead_number(row["prop_snps"]), float(rr["prop_snps"]) * 100
        ph, exp_ph = _lead_number(row["IBD prop_h2"]), float(rr["prop_h2"]) * 100
        if ps is None or abs(ps - exp_ps) > ABS_TOL_PCT:
            failures.append(f"[t4] {row['Method']}/{row['Compartment']}: prop_snps {ps} != {exp_ps:.2f}")
        if ph is None or abs(ph - exp_ph) > ABS_TOL_PCT:
            failures.append(f"[t4] {row['Method']}/{row['Compartment']}: prop_h2 {ph} != {exp_ph:.1f}")

    print(f"raw rows checked        : {len(raw_rows)} (invariants 1 & 2)")
    print(f"table1 cells checked    : {n1}")
    print(f"table4 cells checked    : {n4}")
    if failures:
        print(f"\nFAILED ({len(failures)}):")
        for f in failures:
            print("  " + f)
        return 1
    print("\nAll table values trace to the raw sLDSC output. OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
