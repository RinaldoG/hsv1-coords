#!/usr/bin/env python3
"""Fetch and checksum the HSV1-coords inputs.

Every input is recorded in data/raw/MANIFEST.tsv with its sha256, byte count and the
date fetched, so a release can state exactly what it was built from. Re-running is
idempotent: an input whose checksum already matches the manifest is not re-downloaded.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import pathlib
import sys
import time
import urllib.error
import urllib.request

RAW = pathlib.Path(__file__).resolve().parents[1] / "data" / "raw"
MANIFEST = RAW / "MANIFEST.tsv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# A generic client name on purpose: no personal or project identifying details.
UA = "hsv1-coords/0.1 (data fetch script)"

# accession -> (filename, why we want it)
GENBANK = {
    "NC_001806.2": ("ref_NC_001806.2.gb", "RefSeq reference, 77 CDS, provisional, 2018"),
    "BK012101.1": ("alt_BK012101.1.gb", "300 CDS over the identical 152,222 bases"),
    "JN555585.1": ("base_JN555585.1.gb", "coordinate base of the Whisnant 2020 annotation"),
    "ON783214.1": ("isomer_ON783214.1.gb", "isomer-resolved genome 1 of 4 (PRJNA846339)"),
    "ON783215.1": ("isomer_ON783215.1.gb", "isomer-resolved genome 2 of 4"),
    "ON783216.1": ("isomer_ON783216.1.gb", "isomer-resolved genome 3 of 4"),
    "ON783217.1": ("isomer_ON783217.1.gb", "isomer-resolved genome 4 of 4"),
}

ZENODO = {
    "HSV1_Annotation.gtf": (
        "https://zenodo.org/records/3724023/files/HSV1_Annotation.gtf?download=1",
        "Whisnant et al. 2020 re-annotation, seqname JN555585",
    ),
}


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def get(url: str, dest: pathlib.Path, tries: int = 3) -> None:
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as resp:
                dest.write_bytes(resp.read())
            return
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == tries:
                raise
            print(f"    attempt {attempt} failed ({exc}), retrying", file=sys.stderr)
            time.sleep(2 * attempt)


def fetch_all(force: bool = False) -> list[tuple[str, str, int, str, str]]:
    RAW.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, int, str, str]] = []
    today = dt.date.today().isoformat()

    for acc, (fname, why) in GENBANK.items():
        dest = RAW / fname
        if dest.exists() and not force:
            print(f"  have {fname}")
        else:
            print(f"  fetching {acc}")
            url = f"{EUTILS}?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
            get(url, dest)
            time.sleep(0.4)  # courtesy rate limit for E-utilities
        rows.append((fname, acc, dest.stat().st_size, sha256(dest), f"{today}\t{why}"))

    for fname, (url, why) in ZENODO.items():
        dest = RAW / fname
        if dest.exists() and not force:
            print(f"  have {fname}")
        else:
            print(f"  fetching {fname}")
            get(url, dest)
        rows.append((fname, url, dest.stat().st_size, sha256(dest), f"{today}\t{why}"))

    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download everything")
    args = ap.parse_args()

    rows = fetch_all(force=args.force)
    with MANIFEST.open("w") as fh:
        fh.write("file\tsource\tbytes\tsha256\tfetched\twhy\n")
        for name, src, size, digest, tail in rows:
            fh.write(f"{name}\t{src}\t{size}\t{digest}\t{tail}\n")
    print(f"\n{len(rows)} inputs, manifest written to {MANIFEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
