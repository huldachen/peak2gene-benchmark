"""
Compute LD scores for custom sLDSC annotations.

For each annotation directory under ``data/processed/ldsc/annot/``, runs
``ldsc.py --l2`` over chr1–22 using the 1000G EUR plink files and HapMap3
SNP list, producing ``.l2.ldscore.gz``, ``.l2.M``, ``.l2.M_5_50`` files
alongside the input annotations.

Parallelism
-----------
By default, chromosomes within one annotation run in parallel (4 threads).
Pass ``--serial`` if memory-constrained.

Usage
-----
    python scripts/ldsc/compute_ldscores.py --all
    python scripts/ldsc/compute_ldscores.py --annot-name method1_distance__immune
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

ROOT        = Path(__file__).resolve().parents[2]
PLINK_DIR   = ROOT / "data/raw/ldsc/reference/1000G_EUR_Phase3_plink"
HM3_SNPLIST = ROOT / "data/raw/ldsc/reference/hm3_no_MHC.list.txt"
ANNOT_DIR   = ROOT / "data/processed/ldsc/annot"
LDSC_ENV    = "ldsc"


def _list_annotations() -> list[str]:
    return sorted(p.name for p in ANNOT_DIR.iterdir() if p.is_dir())


def compute_ldscore_one_chrom(annot_name: str, chrom: int) -> tuple[int, bool, str]:
    """Run ldsc.py --l2 for one (annot, chrom)."""
    annot_dir  = ANNOT_DIR / annot_name
    annot_file = annot_dir / f"{annot_name}.{chrom}.annot.gz"
    if not annot_file.exists():
        return chrom, False, f"missing: {annot_file.name}"
    if (annot_dir / f"{annot_name}.{chrom}.l2.ldscore.gz").exists():
        return chrom, True, "already exists"

    cmd = [
        "conda", "run", "-n", LDSC_ENV, "ldsc.py",
        "--l2",
        "--bfile",       str(PLINK_DIR / f"1000G.EUR.QC.{chrom}"),
        "--ld-wind-cm",  "1",
        "--annot",       str(annot_file),
        "--thin-annot",
        "--print-snps",  str(HM3_SNPLIST),
        "--out",         str(annot_dir / f"{annot_name}.{chrom}"),
    ]
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return chrom, False, f"exit {result.returncode}: {result.stderr[-400:]}"
    return chrom, True, f"{time.time() - t0:.0f}s"


def compute_ldscore_all_chroms(annot_name: str, parallel: int = 4) -> None:
    """Run ldsc.py --l2 for chr1–22 of one annotation."""
    t0 = time.time()
    log.info("══ %s (parallel=%d) ══", annot_name, parallel)

    if parallel > 1:
        with ThreadPoolExecutor(max_workers=parallel) as ex:
            for fut in as_completed({ex.submit(compute_ldscore_one_chrom, annot_name, c): c
                                     for c in range(1, 23)}):
                c, ok, msg = fut.result()
                log.info("  chr%-2d %s  %s", c, "OK" if ok else "FAIL", msg)
    else:
        for c in range(1, 23):
            c, ok, msg = compute_ldscore_one_chrom(annot_name, c)
            log.info("  chr%-2d %s  %s", c, "OK" if ok else "FAIL", msg)

    dt = time.time() - t0
    log.info("══ done %s: %.0f s (%.1f min) ══", annot_name, dt, dt / 60)


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--annot-name")
    p.add_argument("--all", action="store_true",
                   help="run for every annotation directory")
    p.add_argument("--parallel", type=int, default=4,
                   help="threads across chromosomes (default 4)")
    p.add_argument("--serial", action="store_true", help="overrides --parallel")
    args = p.parse_args(argv)

    parallel = 1 if args.serial else args.parallel
    if args.all:
        for n in _list_annotations():
            compute_ldscore_all_chroms(n, parallel=parallel)
    else:
        if not args.annot_name:
            sys.exit("--annot-name required unless --all")
        compute_ldscore_all_chroms(args.annot_name, parallel=parallel)


if __name__ == "__main__":
    main()
