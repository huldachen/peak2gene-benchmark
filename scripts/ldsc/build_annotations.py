"""
Build sLDSC annotation files from the peak-gene linkage methods.

For each (method, compartment), produces a BED file of linked peaks in hg19
coordinates and then a per-chromosome ``.annot.gz`` file compatible with
LDSC's ``ldsc.py --l2`` workflow. A SNP receives value 1 if it falls within
any peak that has at least one peak-gene link for that compartment, else 0.

Usage
-----
    python scripts/ldsc/build_annotations.py --all
    python scripts/ldsc/build_annotations.py --method method1 --compartment immune
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
from pyliftover import LiftOver

log = logging.getLogger(__name__)

ROOT          = Path(__file__).resolve().parents[2]
CHAIN         = ROOT / "data/raw/reference/hg38ToHg19.over.chain.gz"
PLINK_DIR     = ROOT / "data/raw/ldsc/reference/1000G_EUR_Phase3_plink"
ANNOT_BED_DIR = ROOT / "data/raw/ldsc/annotations"
ANNOT_OUT_DIR = ROOT / "data/processed/ldsc/annot"
MAKE_ANNOT_PY = ROOT / "external/ldsc/make_annot.py"
LDSC_ENV      = "ldsc"

METHOD_SOURCES = {
    "method1_distance":     "results/linkage/method1",
    "method2_cicero":       "results/linkage/method2",
    "method2b_glasso":      "results/linkage/method2b",
    "method3_abc":          "results/linkage/method3",
    "method4a_paired":      "results/linkage/method4_paired",
    "method4b_crosscohort": "results/linkage/method4_crosscohort",
    "method5_re2g":         "results/linkage/method5",
}
COMPARTMENTS = ["epithelial", "immune", "stromal"]

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)}


_LO_SINGLETON: Optional[LiftOver] = None

def _get_liftover() -> LiftOver:
    global _LO_SINGLETON
    if _LO_SINGLETON is None:
        log.info("Loading liftover chain: %s", CHAIN)
        _LO_SINGLETON = LiftOver(str(CHAIN))
    return _LO_SINGLETON


def liftover_peak(lo: LiftOver, chrom: str, start: int, end: int):
    """Liftover a peak [start, end) on ``chrom``. Returns (chrom19, start19, end19) or None."""
    s = lo.convert_coordinate(chrom, start)
    e = lo.convert_coordinate(chrom, end - 1)   # end-1 because BED is half-open
    if not s or not e:
        return None
    s_chrom, s_pos = s[0][0], s[0][1]
    e_chrom, e_pos = e[0][0], e[0][1]
    if s_chrom != e_chrom:
        return None
    lo_start = min(s_pos, e_pos)
    lo_end   = max(s_pos, e_pos) + 1
    if lo_end <= lo_start:
        return None
    return (s_chrom, lo_start, lo_end)


def method_compartment_to_bed(method: str, compartment: str) -> pd.DataFrame:
    """Return hg19 BED DataFrame of unique linked peaks for (method, compartment)."""
    tsv = ROOT / METHOD_SOURCES[method] / f"{compartment}.tsv"
    if not tsv.exists():
        raise FileNotFoundError(f"Missing linkage TSV: {tsv}")
    df = pd.read_csv(tsv, sep="\t", usecols=["peak_id"])
    peak_ids = df["peak_id"].dropna().unique()
    log.info("[%s/%s] %d unique linked peaks", method, compartment, len(peak_ids))

    records = []
    for pid in peak_ids:
        try:
            c, coords = pid.split(":", 1)
            s, e = coords.split("-")
            records.append((c, int(s), int(e)))
        except Exception:
            continue
    hg38 = pd.DataFrame(records, columns=["chrom", "start", "end"])
    hg38 = hg38[hg38["chrom"].isin(CANONICAL_CHROMS)].copy()

    lo = _get_liftover()
    out = []
    drops = 0
    for _, row in hg38.iterrows():
        lifted = liftover_peak(lo, row["chrom"], int(row["start"]), int(row["end"]))
        if lifted is None:
            drops += 1; continue
        out.append(lifted)
    log.info("[%s/%s] liftover hg38->hg19: kept %d, dropped %d (%.2f%%)",
             method, compartment, len(out), drops,
             100 * drops / max(len(hg38), 1))
    bed = pd.DataFrame(out, columns=["chrom", "start", "end"])
    bed = bed.drop_duplicates().sort_values(["chrom", "start"]).reset_index(drop=True)
    return bed


def run_make_annot(bed_file: Path, bimfile: Path, annot_out: Path) -> None:
    """Call ldsc's make_annot.py on a single-chromosome bim + genome-wide BED."""
    annot_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", LDSC_ENV,
        "python", str(MAKE_ANNOT_PY),
        "--bed-file",  str(bed_file),
        "--bimfile",   str(bimfile),
        "--annot-file", str(annot_out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("make_annot.py failed for %s / %s: %s",
                  bed_file.name, bimfile.name, result.stderr)
        raise RuntimeError(f"make_annot.py failed: {result.stderr}")


def build_annotation(method: str, compartment: str) -> None:
    """Full pipeline for one (method, compartment)."""
    t0 = time.time()
    annot_name = f"{method}__{compartment}"
    log.info("== Building annotation: %s ==", annot_name)

    bed_dir = ANNOT_BED_DIR / annot_name
    bed_dir.mkdir(parents=True, exist_ok=True)
    bed_path = bed_dir / "peaks_hg19.bed"
    if bed_path.exists():
        log.info("[%s] BED already exists, skipping liftover", annot_name)
    else:
        bed = method_compartment_to_bed(method, compartment)
        bed.to_csv(bed_path, sep="\t", header=False, index=False)
        log.info("[%s] wrote %d peaks -> %s", annot_name, len(bed), bed_path)

    out_dir = ANNOT_OUT_DIR / annot_name
    out_dir.mkdir(parents=True, exist_ok=True)
    for chrom in range(1, 23):
        bim = PLINK_DIR / f"1000G.EUR.QC.{chrom}.bim"
        annot = out_dir / f"{annot_name}.{chrom}.annot.gz"
        if annot.exists():
            continue
        run_make_annot(bed_path, bim, annot)

    log.info("== Done: %s  (%.1f s) ==", annot_name, time.time() - t0)


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--method", choices=list(METHOD_SOURCES.keys()))
    p.add_argument("--compartment", choices=COMPARTMENTS)
    p.add_argument("--all", action="store_true",
                   help="Run all methods x compartments.")
    args = p.parse_args(argv)

    if args.all:
        for m in METHOD_SOURCES:
            for c in COMPARTMENTS:
                try:
                    build_annotation(m, c)
                except FileNotFoundError as e:
                    log.warning("Skipping %s / %s: %s", m, c, e)
    else:
        if not args.method or not args.compartment:
            sys.exit("--method and --compartment required unless --all")
        build_annotation(args.method, args.compartment)


if __name__ == "__main__":
    main()
