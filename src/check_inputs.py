#!/usr/bin/env python3
"""Verify data/raw against MANIFEST.tsv, and assert the facts the layer is built on.

This is the M0 gate. It fails loudly rather than warning, so continuous integration
catches an input that moved or an assumption that stopped holding.
"""
from __future__ import annotations

import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fetch_inputs import MANIFEST, RAW, sha256  # noqa: E402

GENOME_BP = 152222


def sequence(path: pathlib.Path) -> str:
    text = path.read_text()
    body = text.split("ORIGIN")[1].split("//")[0]
    return "".join(re.findall(r"[acgtnACGTN]", body.replace("\n", ""))).upper()


def count_feature(path: pathlib.Path, key: str) -> int:
    return sum(1 for line in path.read_text().splitlines()
               if line.startswith(f"     {key} ") or line.startswith(f"     {key}\t"))


def main() -> int:
    failures: list[str] = []

    if not MANIFEST.exists():
        print("MANIFEST.tsv missing, run `make inputs` first", file=sys.stderr)
        return 1

    rows = [l.split("\t") for l in MANIFEST.read_text().splitlines()[1:] if l.strip()]
    for name, _src, size, digest, *_ in rows:
        path = RAW / name
        if not path.exists():
            failures.append(f"{name}: missing")
            continue
        if path.stat().st_size != int(size):
            failures.append(f"{name}: size {path.stat().st_size}, manifest says {size}")
        if sha256(path) != digest:
            failures.append(f"{name}: sha256 does not match the manifest")
    print(f"checked {len(rows)} inputs against the manifest")

    # The assumptions the layer rests on. Each is asserted, not assumed.
    ref = RAW / "ref_NC_001806.2.gb"
    alt = RAW / "alt_BK012101.1.gb"
    base = RAW / "base_JN555585.1.gb"

    checks: list[tuple[str, bool, str]] = []

    ref_seq, base_seq = sequence(ref), sequence(base)
    checks.append(("reference is 152,222 bp", len(ref_seq) == GENOME_BP, str(len(ref_seq))))
    checks.append(("BK012101.1 covers the same length",
                   len(sequence(alt)) == GENOME_BP, str(len(sequence(alt)))))
    # Confirmed on 13 November 2026: these two accessions hold the identical sequence, so
    # the Whisnant 2020 coordinate base differs from the reference only in accession, not
    # in coordinates. If this ever fails, M2 and M3 both need rethinking.
    checks.append(("NC_001806.2 and JN555585.1 are the same sequence",
                   ref_seq == base_seq, "differ"))
    ref_cds, alt_cds = count_feature(ref, "CDS"), count_feature(alt, "CDS")
    checks.append(("reference carries 77 CDS", ref_cds == 77, str(ref_cds)))
    checks.append(("BK012101.1 carries 300 CDS", alt_cds == 300, str(alt_cds)))
    text = ref.read_text()
    checks.append(("TRL is annotated at 1..9213", "repeat_region   1..9213" in text, "absent"))
    checks.append(("LAT's first copy is truncated at the sequence start",
                   "complement(<1..7569)" in text, "absent"))
    checks.append(("LAT is annotated twice", text.count('/gene="LAT"') >= 2,
                   str(text.count('/gene="LAT"'))))

    for label, ok, got in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {label}" + ("" if ok else f" (got {got})"))
        if not ok:
            failures.append(f"assumption broken: {label} (got {got})")

    if failures:
        print("\n" + "\n".join(failures), file=sys.stderr)
        return 1
    print("\nall inputs and assumptions verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
