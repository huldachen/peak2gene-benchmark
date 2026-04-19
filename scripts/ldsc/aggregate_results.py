"""
Aggregate LDSC h² results across all (trait, annotation) runs into one TSV.

For each ``{trait}__{annot_name}.results`` / ``.log`` pair in the per-trait
subdirectories, extracts the custom annotation's row (LDSC category
``L2_0``) plus top-level heritability stats from the log. Output columns:
trait, method, compartment, enrichment, enrichment_se, enrichment_p,
prop_snps, prop_h2, coefficient, coefficient_z, total_h2, intercept.

Usage
-----
    python scripts/ldsc/aggregate_results.py --out results/ldsc/enrichment_table.tsv
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

log = logging.getLogger(__name__)

ROOT   = Path(__file__).resolve().parents[2]
H2_DIR = ROOT / "data/processed/ldsc/results"


def parse_results_row(results_path: Path) -> Optional[dict]:
    """Parse the LDSC ``.results`` file; return the L2_0 (custom annotation) row."""
    if not results_path.exists():
        return None
    df = pd.read_csv(results_path, sep="\t")
    if df.empty:
        return None
    row = df.iloc[0].to_dict()
    return {
        "category":       row.get("Category"),
        "prop_snps":      float(row.get("Prop._SNPs", "nan")),
        "prop_h2":        float(row.get("Prop._h2",   "nan")),
        "prop_h2_se":     float(row.get("Prop._h2_std_error", "nan")),
        "enrichment":     float(row.get("Enrichment", "nan")),
        "enrichment_se":  float(row.get("Enrichment_std_error", "nan")),
        "enrichment_p":   float(row.get("Enrichment_p", "nan")),
        "coefficient":    float(row.get("Coefficient", "nan")),
        "coefficient_se": float(row.get("Coefficient_std_error", "nan")),
        "coefficient_z":  float(row.get("Coefficient_z-score", "nan")),
    }


def parse_log_stats(log_path: Path) -> dict:
    """Extract total h², intercept, λ_GC, SNP count from an LDSC log file."""
    out = {"total_h2": float("nan"), "intercept": float("nan"),
           "lambda_gc": float("nan"), "n_snps": float("nan")}
    if not log_path.exists():
        return out
    text = log_path.read_text()
    patterns = {
        "total_h2":  r"Total Observed scale h2:\s*([\-\d.e]+)",
        "intercept": r"Intercept:\s*([\-\d.e]+)",
        "lambda_gc": r"Lambda GC:\s*([\-\d.e]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            out[key] = float(m.group(1))
    m = re.search(r"After merging with reference.*?(\d+)\s*SNPs remain", text, re.DOTALL)
    if m:
        out["n_snps"] = float(m.group(1))
    return out


def aggregate(h2_dir: Path) -> pd.DataFrame:
    """Walk h2_dir, parse every (trait, annotation) output, return a DataFrame."""
    rows = []
    for trait_dir in sorted(h2_dir.iterdir()):
        if not trait_dir.is_dir():
            continue
        trait = trait_dir.name
        for results_file in sorted(trait_dir.glob("*.results")):
            prefix = f"{trait}__"
            if not results_file.stem.startswith(prefix):
                continue
            annot_name = results_file.stem[len(prefix):]
            method, compartment = (annot_name.rsplit("__", 1)
                                   if "__" in annot_name else (annot_name, "all"))
            row = parse_results_row(results_file)
            if row is None:
                log.warning("Empty or missing: %s", results_file)
                continue
            rows.append({
                "trait": trait, "method": method, "compartment": compartment,
                **row, **parse_log_stats(results_file.with_suffix(".log")),
            })
    df = pd.DataFrame(rows)
    return df.sort_values(["trait", "method", "compartment"]).reset_index(drop=True) if not df.empty else df


def _format_for_print(df: pd.DataFrame) -> str:
    cols = ["trait", "method", "compartment", "enrichment", "enrichment_se",
            "enrichment_p", "prop_snps", "prop_h2", "coefficient_z"]
    sub = df[cols].copy()
    for c in ("enrichment", "enrichment_se", "prop_snps", "prop_h2", "coefficient_z"):
        sub[c] = sub[c].map(lambda x: f"{x:.3f}" if pd.notna(x) else "NA")
    sub["enrichment_p"] = sub["enrichment_p"].map(
        lambda x: f"{x:.2e}" if pd.notna(x) else "NA"
    )
    return sub.to_string(index=False)


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--h2-dir", default=str(H2_DIR))
    p.add_argument("--out", default=str(ROOT / "results/ldsc/enrichment_table.tsv"))
    p.add_argument("--print", action="store_true", dest="do_print")
    args = p.parse_args(argv)

    df = aggregate(Path(args.h2_dir))
    if df.empty:
        sys.exit(f"No .results files found in {args.h2_dir}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    log.info("Wrote %d rows → %s", len(df), args.out)

    if args.do_print:
        print(_format_for_print(df))


if __name__ == "__main__":
    main()
