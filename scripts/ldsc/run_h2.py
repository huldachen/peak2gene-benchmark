"""
Run partitioned heritability (sLDSC --h2) per (trait, annotation).

For each combination, fits sLDSC jointly with our custom annotation and the
baseline-LD v2.2 model (97 categories). Outputs a LDSC log + .results file.

Usage
-----
    python scripts/ldsc/run_h2.py --all
    python scripts/ldsc/run_h2.py --trait IBD --annot-name method1_distance__immune
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

ROOT            = Path(__file__).resolve().parents[2]
SUMSTATS_DIR    = ROOT / "data/raw/ldsc/gwas_sumstats/sumstats_107"
BASELINE_PREFIX = ROOT / "data/raw/ldsc/reference/baselineLD_v2.2/baselineLD."
WEIGHTS_PREFIX  = ROOT / "data/raw/ldsc/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FRQ_PREFIX      = ROOT / "data/raw/ldsc/reference/1000G_Phase3_frq/1000G.EUR.QC."
ANNOT_DIR       = ROOT / "data/processed/ldsc/annot"
H2_OUT_DIR      = ROOT / "data/processed/ldsc/results"
LDSC_ENV        = "ldsc"

TRAITS = {
    "IBD":    "PASS.Inflammatory_Bowel_Disease.deLange2017.sumstats.gz",
    "Height": "PASS.Height.Yengo2022.sumstats.gz",
    "EA":     "PASS.Education_Years.Rietveld2013.sumstats.gz",
}


def _list_annotations() -> list[str]:
    return sorted(p.name for p in ANNOT_DIR.iterdir()
                  if p.is_dir() and any(p.glob("*.1.l2.ldscore.gz")))


def run_h2_one(trait: str, annot_name: str) -> tuple[str, str, bool, str]:
    """Run ldsc.py --h2 for one (trait, annotation)."""
    sumstats = SUMSTATS_DIR / TRAITS[trait]
    annot_prefix = ANNOT_DIR / annot_name / f"{annot_name}."

    if not sumstats.exists():
        return trait, annot_name, False, f"missing sumstats: {sumstats}"
    if not (ANNOT_DIR / annot_name / f"{annot_name}.1.l2.ldscore.gz").exists():
        return trait, annot_name, False, "LD scores not yet computed"

    out_dir = H2_OUT_DIR / trait
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / f"{trait}__{annot_name}"

    # Skip if .results already exists (makes --all properly idempotent)
    if out_prefix.with_suffix(".results").exists():
        return trait, annot_name, True, "already exists"

    cmd = [
        "conda", "run", "-n", LDSC_ENV, "ldsc.py",
        "--h2",           str(sumstats),
        "--ref-ld-chr",   f"{annot_prefix},{BASELINE_PREFIX}",
        "--w-ld-chr",     str(WEIGHTS_PREFIX),
        "--frqfile-chr",  str(FRQ_PREFIX),
        "--overlap-annot",
        "--print-coefficients",
        "--print-delete-vals",
        "--out",          str(out_prefix),
    ]
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return trait, annot_name, False, f"exit {result.returncode}: {result.stderr[-400:]}"
    return trait, annot_name, True, f"{time.time() - t0:.0f}s"


def run_all(parallel: int = 4) -> None:
    traits, annots = list(TRAITS), _list_annotations()
    log.info("Running %d combinations (%d traits × %d annotations), parallel=%d",
             len(traits) * len(annots), len(traits), len(annots), parallel)
    jobs = [(t, a) for t in traits for a in annots]
    with ThreadPoolExecutor(max_workers=parallel) as ex:
        for fut in as_completed({ex.submit(run_h2_one, t, a): (t, a) for t, a in jobs}):
            t, a, ok, msg = fut.result()
            log.info("  %-6s %s  %-50s  %s", t, "OK" if ok else "FAIL", a, msg)


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--trait", choices=list(TRAITS))
    p.add_argument("--annot-name")
    p.add_argument("--all", action="store_true")
    p.add_argument("--parallel", type=int, default=4)
    args = p.parse_args(argv)

    if args.all:
        run_all(parallel=args.parallel)
    else:
        if not args.trait or not args.annot_name:
            sys.exit("--trait and --annot-name required unless --all")
        t, a, ok, msg = run_h2_one(args.trait, args.annot_name)
        log.info("%s / %s → %s (%s)", t, a, "OK" if ok else "FAIL", msg)


if __name__ == "__main__":
    main()
