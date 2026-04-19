"""
Build a continuous-valued sLDSC annotation for Method 5 (ENCODE-rE2G).

Each SNP's annotation value is the maximum rE2G probability among peaks
containing it, where each peak's score is the max rE2G probability across
its gene candidates. Output format is a thin-annot ``.annot.gz`` with a
single ``ANNOT`` column of floats in [0, 1] per chromosome.

Usage
-----
    python scripts/ldsc/build_method5_continuous.py --compartment immune
    python scripts/ldsc/build_method5_continuous.py --all
"""
from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from pyliftover import LiftOver

log = logging.getLogger(__name__)

ROOT          = Path(__file__).resolve().parents[2]
CHAIN         = ROOT / "data/raw/reference/hg38ToHg19.over.chain.gz"
PLINK_DIR     = ROOT / "data/raw/ldsc/reference/1000G_EUR_Phase3_plink"
ANNOT_OUT_DIR = ROOT / "data/processed/ldsc/annot"
METHOD5_DIR   = ROOT / "results/linkage/method5"
COMPARTMENTS  = ["epithelial", "immune", "stromal"]

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)}


def _get_liftover(path: Path = CHAIN) -> LiftOver:
    return LiftOver(str(path))


def load_method5_best_per_peak(compartment: str) -> pd.DataFrame:
    """Load Method 5's TSV and return (peak_id, best_probability) per peak."""
    tsv = METHOD5_DIR / f"{compartment}.tsv"
    if not tsv.exists():
        raise FileNotFoundError(tsv)
    df = pd.read_csv(tsv, sep="\t", usecols=["peak_id", "rE2G_probability"])
    best = (df.groupby("peak_id", sort=False)["rE2G_probability"]
              .max()
              .reset_index()
              .rename(columns={"rE2G_probability": "score"}))
    log.info("[%s] %d unique peaks with rE2G probability", compartment, len(best))
    log.info("[%s]   score distribution: min=%.4f p25=%.4f p50=%.4f p75=%.4f max=%.4f",
             compartment, best["score"].min(), best["score"].quantile(0.25),
             best["score"].quantile(0.5), best["score"].quantile(0.75),
             best["score"].max())
    return best


def liftover_peaks_with_scores(peaks: pd.DataFrame, lo: LiftOver) -> pd.DataFrame:
    """Liftover peak_id (hg38) to hg19, keeping the rE2G score."""
    records = []
    drops = 0
    for _, row in peaks.iterrows():
        pid = row["peak_id"]
        try:
            c, coords = pid.split(":", 1)
            s, e = [int(x) for x in coords.split("-")]
        except Exception:
            drops += 1; continue
        if c not in CANONICAL_CHROMS:
            drops += 1; continue
        s_map = lo.convert_coordinate(c, s)
        e_map = lo.convert_coordinate(c, e - 1)
        if not s_map or not e_map:
            drops += 1; continue
        sc, sp = s_map[0][0], s_map[0][1]
        ec, ep = e_map[0][0], e_map[0][1]
        if sc != ec:
            drops += 1; continue
        start19 = min(sp, ep)
        end19   = max(sp, ep) + 1
        if end19 <= start19:
            drops += 1; continue
        records.append((sc, start19, end19, float(row["score"])))
    log.info("liftover: %d kept, %d dropped (%.2f%%)",
             len(records), drops, 100 * drops / max(len(peaks), 1))
    return pd.DataFrame(records, columns=["chrom", "start", "end", "score"])


def _build_snp_annotations_for_chrom(
    peaks_hg19: pd.DataFrame,
    bim_path: Path,
    chrom: str,
) -> pd.DataFrame:
    """Return a DataFrame matching the bim file's SNPs with ANNOT column."""
    bim = pd.read_csv(bim_path, sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    chrom_peaks = peaks_hg19[peaks_hg19["chrom"] == chrom].copy()
    if chrom_peaks.empty:
        bim["ANNOT"] = 0.0
        return bim[["CHR", "BP", "SNP", "CM", "ANNOT"]]

    # Two-pointer sweep: for each peak, SNPs in [start, end) take max(annot, score)
    chrom_peaks = chrom_peaks.sort_values("start").reset_index(drop=True)
    snps = bim.sort_values("BP").reset_index(drop=True)
    annot = np.zeros(len(snps), dtype=np.float64)
    bps = snps["BP"].values

    for _, peak in chrom_peaks.iterrows():
        lo = np.searchsorted(bps, peak["start"], side="left")
        hi = np.searchsorted(bps, peak["end"], side="left")
        if hi <= lo:
            continue
        annot[lo:hi] = np.maximum(annot[lo:hi], peak["score"])

    snps["ANNOT"] = annot
    # Restore bim original order (ldsc expects bim order)
    result = (snps.set_index("SNP").reindex(bim["SNP"])
                  .reset_index())
    return result[["CHR", "BP", "SNP", "CM", "ANNOT"]]


def build_compartment(compartment: str) -> None:
    t0 = time.time()
    annot_name = f"method5_re2g_continuous__{compartment}"
    out_dir = ANNOT_OUT_DIR / annot_name
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("== Building %s ==", annot_name)

    peaks = load_method5_best_per_peak(compartment)
    lo = _get_liftover()
    peaks19 = liftover_peaks_with_scores(peaks, lo)

    bed_with_score = out_dir / "peaks_hg19_with_score.bed"
    peaks19.to_csv(bed_with_score, sep="\t", header=False, index=False)
    log.info("[%s] wrote %d lifted peaks -> %s", annot_name, len(peaks19), bed_with_score)

    for chrom_i in range(1, 23):
        chrom = f"chr{chrom_i}"
        bim = PLINK_DIR / f"1000G.EUR.QC.{chrom_i}.bim"
        annot_df = _build_snp_annotations_for_chrom(peaks19, bim, chrom)

        # Thin-annot format matches the rest of the pipeline (compute_ldscores uses --thin-annot)
        out = out_dir / f"{annot_name}.{chrom_i}.annot.gz"
        with gzip.open(out, "wt") as fh:
            fh.write("ANNOT\n")
            for v in annot_df["ANNOT"]:
                fh.write(f"{v:.6f}\n")

        nonzero = int((annot_df["ANNOT"] > 0).sum())
        log.info("[%s] chr%d: %d SNPs, %d non-zero (%.2f%%), mean=%.4f",
                 annot_name, chrom_i, len(annot_df), nonzero,
                 100 * nonzero / max(len(annot_df), 1),
                 float(annot_df["ANNOT"].mean()))

    log.info("== Done %s (%.1f s) ==", annot_name, time.time() - t0)


def main(argv: Optional[list[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        datefmt="%H:%M:%S")
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--compartment", choices=COMPARTMENTS)
    p.add_argument("--all", action="store_true")
    args = p.parse_args(argv)

    if args.all:
        for c in COMPARTMENTS:
            build_compartment(c)
    else:
        if not args.compartment:
            sys.exit("--compartment required unless --all")
        build_compartment(args.compartment)


if __name__ == "__main__":
    main()
