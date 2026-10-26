#!/usr/bin/env python3
"""
Claim 1 with uncertainty: is ArchR's raw-Pearson shortcut (Method 2a) equivalent
to canonical Cicero GLASSO (Method 2b)?

The claim is a NON-equivalence claim, so it has to survive its own error bars.
This script attaches an interval to both halves of it:

  (a) the per-SNP enrichment shift (2a -> 2b), reported as a difference in
      enrichment units with a conservative SE, and
  (b) the tau* coefficient z-score, which is already uncertainty-scaled
      (z = coefficient / SE), so |z| > 1.96 is the significance test.

On the SE of the difference: the two annotations are built from the same peaks
and overlap heavily, so the two sLDSC estimates are positively correlated. The
variance of their difference is then sa^2 + sb^2 - 2*rho*sa*sb, which is SMALLER
than the independent sum. We do not know rho, so we report the independent
propagation sqrt(sa^2 + sb^2) and label it what it is: a CONSERVATIVE upper
bound. A shift that is not significant under this bound is not significant under
any plausible rho; a shift that IS significant here would need rho to be
confirmed before leaning on it.

Run:  python scripts/ldsc/claim1_equivalence.py
Out:  results/tables/table5_claim1_equivalence.tsv  + stdout summary
"""

from __future__ import annotations

import csv
import math
import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
SRC = os.path.join(PROJECT_ROOT, "results", "ldsc", "enrichment_table.tsv")
OUT = os.path.join(PROJECT_ROOT, "results", "tables", "table5_claim1_equivalence.tsv")

PEARSON = "method2_cicero"      # 2a - ArchR raw-Pearson shortcut
GLASSO = "method2b_glasso"      # 2b - canonical Cicero GLASSO
TRAITS = ("IBD", "Height", "EA")
COMPARTMENTS = ("epithelial", "immune", "stromal")
Z_CRIT = 1.96                   # two-sided 95%


def main() -> int:
    rows = list(csv.DictReader(open(SRC), delimiter="\t"))
    idx = {(r["trait"], r["compartment"], r["method"]): r for r in rows}

    out = []
    for trait in TRAITS:
        for comp in COMPARTMENTS:
            a = idx.get((trait, comp, PEARSON))
            b = idx.get((trait, comp, GLASSO))
            if not (a and b):
                continue
            ea, sa = float(a["enrichment"]), float(a["enrichment_se"])
            eb, sb = float(b["enrichment"]), float(b["enrichment_se"])
            delta = eb - ea
            delta_se = math.hypot(sa, sb)          # conservative: assumes rho = 0
            delta_z = delta / delta_se if delta_se else float("nan")
            za, zb = float(a["coefficient_z"]), float(b["coefficient_z"])

            enr_sig = abs(delta_z) > Z_CRIT
            # tau* non-equivalence: one annotation significant, the other not,
            # or the two significant with opposite sign.
            tau_flip = (abs(za) > Z_CRIT) != (abs(zb) > Z_CRIT) or (
                abs(za) > Z_CRIT and abs(zb) > Z_CRIT and za * zb < 0
            )
            out.append({
                "trait": trait,
                "compartment": comp,
                "enr_2a_pearson": f"{ea:.2f}",
                "enr_2a_se": f"{sa:.2f}",
                "enr_2b_glasso": f"{eb:.2f}",
                "enr_2b_se": f"{sb:.2f}",
                "enr_delta": f"{delta:+.2f}",
                "enr_delta_se_conservative": f"{delta_se:.2f}",
                "enr_delta_z": f"{delta_z:+.2f}",
                "enr_shift_significant": str(enr_sig),
                "tau_z_2a": f"{za:+.2f}",
                "tau_z_2b": f"{zb:+.2f}",
                "tau_nonequivalent": str(tau_flip),
            })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print("Claim 1 - Method 2a (raw Pearson) vs 2b (Cicero GLASSO), with intervals")
    print("enrichment shift = 2b - 2a, SE conservative (independent propagation)\n")
    hdr = f"{'trait':<7}{'compartment':<12}{'2a':>14}{'2b':>14}{'shift':>18}{'sig?':>6}   tau* 2a -> 2b"
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        print(f"{r['trait']:<7}{r['compartment']:<12}"
              f"{r['enr_2a_pearson']:>7}±{r['enr_2a_se']:<6}"
              f"{r['enr_2b_glasso']:>7}±{r['enr_2b_se']:<6}"
              f"{r['enr_delta']:>8} ±{r['enr_delta_se_conservative']:<7}"
              f"{'YES' if r['enr_shift_significant']=='True' else 'no':>5}   "
              f"{r['tau_z_2a']} -> {r['tau_z_2b']}"
              f"{'   [tau* NON-EQUIV]' if r['tau_nonequivalent']=='True' else ''}")

    n_enr = sum(r["enr_shift_significant"] == "True" for r in out)
    n_tau = sum(r["tau_nonequivalent"] == "True" for r in out)
    print(f"\n  enrichment shifts significant at 95%: {n_enr}/{len(out)}")
    print(f"  tau* non-equivalent cells:            {n_tau}/{len(out)}")
    print(f"\nSaved -> {os.path.relpath(OUT, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
