#!/usr/bin/env python3
"""Fetch and checksum the HSV1-coords inputs.

Every input is recorded in data/raw/MANIFEST.tsv with its sha256, byte count and the
date this content was first recorded, so a release can state exactly what it was built
from. Re-running is idempotent: an input already present is not re-downloaded, and a row
changes only when the content behind it changes, never because the script ran again.
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
    # The previous version of the reference. Published work is defined against both,
    # and the offset between them is not constant, so a coordinate can be wrong purely
    # because two papers cited two versions of one accession.
    "NC_001806.1": ("ref_NC_001806.1.gb", "the previous reference version, 152,261 bp"),
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


def recorded() -> dict[str, tuple[str, str]]:
    """What the committed manifest already says about each file: (sha256, date).

    Round 11: this script wrote today's date into every row whether or not it downloaded
    anything, so the manifest changed on the clock alone. CI's step "Fail if any input
    changed upstream" is `git diff --exit-code` on this file immediately after running
    this script, so from the day after the last commit it would have failed on every
    push, pull request and weekly cron, and the five steps after it, including the whole
    test suite, would never have run. An actual upstream change would also have been
    indistinguishable from the date moving.

    Round 11's fix keyed the date on whether the file was downloaded, and round 12 found
    that this does not help where it was needed. `data/raw/` is gitignored, so a CI
    checkout has no inputs at all, so every input IS genuinely fetched, so every row got
    today's date and the diff gate still failed on the clock. Measured on a fresh clone
    with the network stubbed: all nine rows rewritten, `git diff --exit-code` exit 1.

    The date is therefore keyed on the CONTENT, not on whether a download happened: a
    file whose checksum matches the committed row keeps that row's date, whoever fetched
    it and whenever. A file whose checksum differs is new content, so it gets today and
    the row changes, which is exactly the upstream change the gate exists to catch. The
    column's meaning is now the date this content was first recorded.
    """
    manifest = RAW / "MANIFEST.tsv"
    if not manifest.exists():
        return {}
    out: dict[str, tuple[str, str]] = {}
    for line in manifest.read_text().splitlines()[1:]:
        cells = line.split("\t")
        if len(cells) >= 5:
            out[cells[0]] = (cells[3], cells[4])
    return out


def fetch_all(force: bool = False) -> list[tuple[str, str, int, str, str]]:
    RAW.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, int, str, str]] = []
    today = dt.date.today().isoformat()
    known = recorded()

    def first_recorded(fname: str, digest: str) -> str:
        """Today only if this content is not what the committed manifest already holds."""
        was = known.get(fname)
        return was[1] if was and was[0] == digest else today

    for acc, (fname, why) in GENBANK.items():
        dest = RAW / fname
        if dest.exists() and not force:
            print(f"  have {fname}")
        else:
            print(f"  fetching {acc}")
            url = f"{EUTILS}?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
            get(url, dest)
            time.sleep(0.4)  # courtesy rate limit for E-utilities
        digest = sha256(dest)
        rows.append((fname, acc, dest.stat().st_size, digest,
                     f"{first_recorded(fname, digest)}\t{why}"))

    for fname, (url, why) in ZENODO.items():
        dest = RAW / fname
        if dest.exists() and not force:
            print(f"  have {fname}")
        else:
            print(f"  fetching {fname}")
            get(url, dest)
        digest = sha256(dest)
        rows.append((fname, url, dest.stat().st_size, digest,
                     f"{first_recorded(fname, digest)}\t{why}"))

    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download everything")
    args = ap.parse_args()

    rows = fetch_all(force=args.force)
    with MANIFEST.open("w") as fh:
        fh.write("file\tsource\tbytes\tsha256\tfirst_recorded\twhy\n")
        for name, src, size, digest, tail in rows:
            fh.write(f"{name}\t{src}\t{size}\t{digest}\t{tail}\n")
    print(f"\n{len(rows)} inputs, manifest written to {MANIFEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
