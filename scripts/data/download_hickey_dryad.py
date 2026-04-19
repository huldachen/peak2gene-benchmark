"""
Download Hickey et al. 2023 processed data from Dryad (CC0 license).

Dryad uses Cloudflare bot protection that blocks plain Python requests;
this script uses a browser-exported cookies.txt to authenticate.

Setup (one time):
    1. Install "Get cookies.txt LOCALLY" (Chrome) or "Cookie Quick Manager" (Firefox).
    2. Open https://datadryad.org/stash/dataset/doi:10.5061/dryad.8pk0p2ns8
    3. Export cookies -> save as data/raw/hickey2023/dryad_cookies.txt

Usage
-----
    python scripts/data/download_hickey_dryad.py metadata
    python scripts/data/download_hickey_dryad.py rna --samples B001-A-301 B001-A-401
    python scripts/data/download_hickey_dryad.py atac --samples B001-A-301
    python scripts/data/download_hickey_dryad.py peaks
    python scripts/data/download_hickey_dryad.py list
    python scripts/data/download_hickey_dryad.py check-cookies
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

BASE      = "https://datadryad.org"
PROJ_ROOT = Path(__file__).resolve().parents[2]
COOKIE_F  = PROJ_ROOT / "data/raw/hickey2023/dryad_cookies.txt"

# {filename: (file_id, dest_subdir)}
METADATA_FILES = {
    "sample_location_metadata.csv":          (2121314, "metadata"),
    "rna_README.md":                         (2124085, "metadata"),
    "Colon_Epithelial_UMAP_CellType.tsv":    (2121305, "cell_types"),
    "Immune_UMAP_CellType.tsv":              (2124088, "cell_types"),
    "Stromal_UMAP_CellType.tsv":             (2124089, "cell_types"),
    "atac_sample_location_metadata.csv":     (2186702, "metadata"),
    "atac_README.md":                        (2186706, "metadata"),
    "scATAC_multiome_cell_types_epithelial_colon.tsv": (2186701, "cell_types"),
    "scATAC_multiome_cell_types_immune.tsv": (2186699, "cell_types"),
    "scATAC_multiome_cell_types_stromal.tsv":(2186700, "cell_types"),
    "colon_epithelial_peaks.bed":            (2186666, "atac"),
    "immune_peaks.bed":                      (2186656, "atac"),
    "stromal_peaks.bed":                     (2186664, "atac"),
}

# Colon RNA samples only (RNA cohort, donors B001/B004/B005). Small-intestine
# samples (Duodenum/Ileum/Jejunum) are excluded.
RNA_SAMPLE_FILES: dict[str, dict[str, int]] = {
    "B001-A-301": {"barcodes.tsv.gz": 2121080, "features.tsv.gz": 2121079, "matrix.mtx.gz": 2121081},
    "B001-A-401": {"barcodes.tsv.gz": 2121078, "features.tsv.gz": 2121076, "matrix.mtx.gz": 2121077},
    "B004-A-204": {"barcodes.tsv.gz": 2121115, "features.tsv.gz": 2121107, "matrix.mtx.gz": 2121114},
    "B005-A-201": {"barcodes.tsv.gz": 2121174, "features.tsv.gz": 2121169, "matrix.mtx.gz": 2121172},
}

# Colon ATAC samples (RNA cohort fragment files). Multiome-cohort fragment
# IDs are not enumerated here; Multiome analyses use the aggregated peak
# matrices below.
ATAC_SAMPLE_FILES: dict[str, dict[str, int]] = {
    "B001-A-301": {"atac_fragments.tsv.gz": 2186636, "atac_fragments.tsv.gz.tbi": 2186633},
    "B001-A-401": {"atac_fragments.tsv.gz": 2186639, "atac_fragments.tsv.gz.tbi": 2186638},
    "B004-A-204": {"atac_fragments.tsv.gz": 2186645, "atac_fragments.tsv.gz.tbi": 2186644},
    "B005-A-201": {"atac_fragments.tsv.gz": 2186609, "atac_fragments.tsv.gz.tbi": 2186607},
}

PEAK_MATRIX_FILES = {
    "colon_epithelial_peak_matrix_cells.tsv": (2186667, "atac"),
    "colon_epithelial_peak_matrix.mtx":       (2186668, "atac"),
    "immune_peak_matrix_cells.tsv":           (2186658, "atac"),
    "immune_peak_matrix.mtx":                 (2186659, "atac"),
    "stromal_peak_matrix_cells.tsv":          (2186663, "atac"),
    "stromal_peak_matrix.mtx":               (2186665, "atac"),
}


def load_cookies(cookie_file: Path) -> dict[str, str]:
    """Parse Netscape-format cookies.txt into a dict."""
    cookies: dict[str, str] = {}
    if not cookie_file.exists():
        return cookies
    with open(cookie_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 7:
                cookies[parts[5]] = parts[6]
    return cookies


def download_file(file_id: int, dest: Path, cookies: dict, label: str = "") -> bool:
    """Download a single file from Dryad by file ID."""
    url = f"{BASE}/api/v2/files/{file_id}/download"
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists():
        print(f"  [skip] {label or dest.name} (already exists, {dest.stat().st_size/1e6:.1f} MB)")
        return True

    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
        "Accept": "*/*",
        "Referer": "https://datadryad.org/",
    }

    try:
        r = requests.get(url, headers=headers, cookies=cookies,
                         stream=True, timeout=60, allow_redirects=True)
        if r.status_code != 200:
            print(f"  [ERROR] {label}: HTTP {r.status_code}")
            if r.status_code == 401:
                print("         -> Need Dryad cookies. See script header for setup.")
            return False

        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded / total
                    mb_done = downloaded / 1e6
                    mb_total = total / 1e6
                    print(f"\r  {label or dest.name}: {pct:.1f}%  ({mb_done:.0f}/{mb_total:.0f} MB)", end="", flush=True)
        print(f"\r  OK {label or dest.name}  ({downloaded/1e6:.1f} MB)                    ")
        return True
    except Exception as e:
        print(f"  [ERROR] {label}: {e}")
        if dest.exists():
            dest.unlink()
        return False


def cmd_metadata(args):
    cookies = load_cookies(COOKIE_F)
    dest_root = PROJ_ROOT / "data/raw/hickey2023"
    print(f"Downloading metadata files ({len(METADATA_FILES)} files) -> {dest_root}")
    ok = 0
    for fname, (fid, subdir) in METADATA_FILES.items():
        dest = dest_root / subdir / fname
        if download_file(fid, dest, cookies, fname):
            ok += 1
        time.sleep(0.3)
    print(f"\n{ok}/{len(METADATA_FILES)} metadata files downloaded.")


def cmd_rna(args):
    cookies = load_cookies(COOKIE_F)
    samples = args.samples or list(RNA_SAMPLE_FILES.keys())
    dest_root = PROJ_ROOT / "data/raw/hickey2023/rna"
    print(f"Downloading scRNA MTX trios for {len(samples)} samples -> {dest_root}")
    for sample in samples:
        if sample not in RNA_SAMPLE_FILES:
            print(f"  [unknown] {sample}")
            continue
        print(f"\n  Sample: {sample}")
        for fname, fid in RNA_SAMPLE_FILES[sample].items():
            dest = dest_root / sample / fname
            download_file(fid, dest, cookies, f"{sample}/{fname}")
            time.sleep(0.3)


def cmd_atac(args):
    cookies = load_cookies(COOKIE_F)
    samples = args.samples or list(ATAC_SAMPLE_FILES.keys())
    dest_root = PROJ_ROOT / "data/raw/hickey2023/atac"
    print(f"Downloading scATAC fragment files for {len(samples)} samples -> {dest_root}")
    for sample in samples:
        if sample not in ATAC_SAMPLE_FILES:
            print(f"  [unknown] {sample}")
            continue
        print(f"\n  Sample: {sample}")
        for fname, fid in ATAC_SAMPLE_FILES[sample].items():
            dest = dest_root / sample / fname
            download_file(fid, dest, cookies, f"{sample}/{fname}")
            time.sleep(0.3)


def cmd_peaks(args):
    cookies = load_cookies(COOKIE_F)
    dest_root = PROJ_ROOT / "data/raw/hickey2023"
    print(f"Downloading aggregated peak matrices ({len(PEAK_MATRIX_FILES)} files) -> {dest_root}")
    for fname, (fid, subdir) in PEAK_MATRIX_FILES.items():
        dest = dest_root / subdir / fname
        download_file(fid, dest, cookies, fname)
        time.sleep(0.3)


def cmd_check_cookies(args):
    """Verify cookie file exists, is parseable, and test one download."""
    if not COOKIE_F.exists():
        print(f"Cookie file NOT found: {COOKIE_F}")
        print("\nSetup instructions:")
        print("  1. Install 'Get cookies.txt LOCALLY' Chrome extension")
        print("  2. Open: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8pk0p2ns8")
        print(f"  3. Export cookies -> save as: {COOKIE_F}")
        return

    cookies = load_cookies(COOKIE_F)
    print(f"Cookie file found: {COOKIE_F}")
    print(f"Parsed {len(cookies)} cookies")
    dryad_cookies = {k: v for k, v in cookies.items()
                     if any(x in k.lower() for x in ["session", "dryad", "dash", "aws"])}
    print(f"Dryad-relevant cookies: {list(dryad_cookies.keys())}")

    print("\nTesting download of README.md (tiny file)...")
    test_dest = PROJ_ROOT / "data/raw/hickey2023/metadata/rna_README.md"
    if test_dest.exists():
        test_dest.unlink()
    ok = download_file(2124085, test_dest, cookies, "rna_README.md")
    if ok and test_dest.exists() and test_dest.stat().st_size > 100:
        print(f"\nOK Cookie authentication working. File: {test_dest.stat().st_size} bytes")
        print(f"  Content preview: {test_dest.read_text()[:200]}")
    else:
        print("\nDownload failed - cookies may be expired or invalid.")
        print("  Try re-exporting cookies after refreshing the Dryad page.")


def cmd_list(args):
    print("=== Metadata files ===")
    for fn, (fid, _) in METADATA_FILES.items():
        print(f"  {fn:<55} https://datadryad.org/downloads/file_stream/{fid}")
    print("\n=== Available colon RNA samples ===")
    for s in RNA_SAMPLE_FILES:
        print(f"  {s}")
    print("\n=== Available colon ATAC samples ===")
    for s in ATAC_SAMPLE_FILES:
        print(f"  {s}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="step")

    sub.add_parser("metadata", help="Download metadata and cell-type annotation files")
    sub.add_parser("list",     help="Print all available samples and download URLs")
    sub.add_parser("peaks",    help="Download aggregated peak matrices by compartment")

    p_rna = sub.add_parser("rna",  help="Download scRNA MTX trios for specified samples")
    p_rna.add_argument("--samples", nargs="+", default=None,
                       help="Sample IDs (e.g. B001-A-201). Default: all colon samples.")

    p_atac = sub.add_parser("atac", help="Download scATAC fragment files for specified samples")
    p_atac.add_argument("--samples", nargs="+", default=None,
                        help="Sample IDs (e.g. B001-A-301). Default: all colon samples.")

    sub.add_parser("check-cookies", help="Verify cookie file is valid and test a download")

    args = parser.parse_args()

    dispatch = {"metadata": cmd_metadata, "rna": cmd_rna,
                "atac": cmd_atac, "peaks": cmd_peaks, "list": cmd_list,
                "check-cookies": cmd_check_cookies}

    if args.step not in dispatch:
        parser.print_help()
        sys.exit(1)

    dispatch[args.step](args)


if __name__ == "__main__":
    main()
