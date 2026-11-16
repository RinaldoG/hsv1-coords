#!/usr/bin/env python3
"""Verify data/raw against MANIFEST.tsv, and assert the facts the layer is built on.

This is the M0 gate. It fails loudly rather than warning, so continuous integration
catches an input that moved or an assumption that stopped holding.

Why it is shaped as two tables. Round 12 found that nothing in the test suite executed
this file, so every line in it could be made inert without a single test noticing:
changing `ref_cds == 77` to `ref_cds == ref_cds` left both the suite green AND `make
check` printing "ok   reference carries 77 CDS" and "all inputs and assumptions
verified". The gate that exists to catch an assumption breaking reported success while
asserting nothing, and the same edit worked on any of the eight assumptions.

So the measurement and the number it is compared against are now separated.
`MEASUREMENTS` holds a function per assumption and nothing else; `EXPECTED` holds a
literal per assumption and nothing else; `assumptions()` is the only thing that compares
them, once, by key. A tautology now requires editing the literal table, which is a
visible diff of a named number, and `tests/test_m5_conventions.py` runs every
measurement against the real records, asserts each literal is a plain int or string,
asserts no measurement reads the expectation table, and asserts that each measurement
disagrees with a deliberately wrong expectation.
"""
from __future__ import annotations

import pathlib
import re
import sys
from typing import Callable, NamedTuple

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fetch_inputs import MANIFEST, RAW, sha256  # noqa: E402

GENOME_BP = 152222

# data/raw files that are deliberately not in MANIFEST.tsv, each with the artefact that
# pins it instead. Round 12: this gate's loop iterated the manifest's own rows, so a
# file sitting in data/raw with no row was never checksummed and never missed, which is
# how genomes.fasta came to be outside this gate entirely.
PINNED_ELSEWHERE = {
    "genomes.fasta": pathlib.Path("data") / "genomes.sha256",
}
NOT_INPUTS = {".gitkeep", "MANIFEST.tsv"}


def sequence(path: pathlib.Path) -> str:
    text = path.read_text()
    body = text.split("ORIGIN")[1].split("//")[0]
    return "".join(re.findall(r"[acgtnACGTN]", body.replace("\n", ""))).upper()


def count_feature(path: pathlib.Path, key: str) -> int:
    return sum(1 for line in path.read_text().splitlines()
               if line.startswith(f"     {key} ") or line.startswith(f"     {key}\t"))


def _ref() -> pathlib.Path:
    return RAW / "ref_NC_001806.2.gb"


def _alt() -> pathlib.Path:
    return RAW / "alt_BK012101.1.gb"


def _base() -> pathlib.Path:
    return RAW / "base_JN555585.1.gb"


# ---- the two tables ----------------------------------------------------------------
# One measurement per assumption. Each returns a measured value and compares nothing.

MEASUREMENTS: dict[str, Callable[[], object]] = {
    "ref_length": lambda: len(sequence(_ref())),
    "alt_length": lambda: len(sequence(_alt())),
    "ref_base_relationship": lambda: ("identical" if sequence(_ref()) == sequence(_base())
                                      else "different"),
    "ref_cds": lambda: count_feature(_ref(), "CDS"),
    "alt_cds": lambda: count_feature(_alt(), "CDS"),
    "trl_annotation": lambda: ("repeat_region   1..9213" in _ref().read_text()),
    "lat_first_copy_truncated": lambda: ("complement(<1..7569)" in _ref().read_text()),
    "lat_gene_qualifiers": lambda: _ref().read_text().count('/gene="LAT"'),
}

# Which record each measurement reads. Round 14: the sensitivity tests moved the input
# and asked only whether the value CHANGED, never which file produced it.
# `base_JN555585.1.gb` holds the identical sequence to the reference, 152,222 bp, 77
# CDS, `repeat_region   1..9213` and 6 `/gene="LAT"`, so swapping `_ref()` for
# `_base()` in four measurements left `make check` certifying four assumptions about
# NC_001806.2 against a different accession, with all 301 tests green and an upstream
# change to NC_001806.2 no longer able to move them.
#
# This is a declaration, so the suite checks it rather than trusting it: stubbing one
# record must move exactly the measurements that name it and no others.
READS: dict[str, tuple[str, ...]] = {
    "ref_length": ("ref_NC_001806.2.gb",),
    "alt_length": ("alt_BK012101.1.gb",),
    "ref_base_relationship": ("ref_NC_001806.2.gb", "base_JN555585.1.gb"),
    "ref_cds": ("ref_NC_001806.2.gb",),
    "alt_cds": ("alt_BK012101.1.gb",),
    "trl_annotation": ("ref_NC_001806.2.gb",),
    "lat_first_copy_truncated": ("ref_NC_001806.2.gb",),
    "lat_gene_qualifiers": ("ref_NC_001806.2.gb",),
}

# One literal per assumption, and the sentence it is stated as. No expression here may
# mention a measurement: that is the whole point of the split.
EXPECTED: dict[str, tuple[object, str]] = {
    "ref_length": (GENOME_BP, "the reference is 152,222 bp"),
    "alt_length": (GENOME_BP, "BK012101.1 covers the same length"),
    # Confirmed on 13 November 2026: these two accessions hold the identical sequence, so
    # the Whisnant 2020 coordinate base differs from the reference only in accession, not
    # in coordinates. If this ever fails, M2 and M3 both need rethinking.
    "ref_base_relationship": ("identical",
                              "NC_001806.2 and JN555585.1 are the same sequence"),
    "ref_cds": (77, "the reference carries 77 CDS"),
    "alt_cds": (300, "BK012101.1 carries 300 CDS"),
    "trl_annotation": (True, "TRL is annotated at 1..9213"),
    "lat_first_copy_truncated": (True,
                                 "LAT's first copy is truncated at the sequence start"),
    # Round 12: this was `>= 2` against a true count of 6, which is a budget of four
    # qualifiers. The count is pinned exactly, like every other number here.
    "lat_gene_qualifiers": (6, "LAT carries 6 /gene qualifiers across its copies"),
}


class Assumption(NamedTuple):
    key: str
    label: str
    got: object
    want: object

    @property
    def ok(self) -> bool:
        return self.got == self.want


def assumptions() -> list[Assumption]:
    """Every assumption, measured and compared. The only comparison in this file."""
    if set(MEASUREMENTS) != set(EXPECTED) or set(MEASUREMENTS) != set(READS):
        missing = (set(MEASUREMENTS) ^ set(EXPECTED)) | (set(MEASUREMENTS) ^ set(READS))
        raise RuntimeError(f"the measurement, expectation and reads tables disagree: "
                           f"{missing}")
    out = []
    for key, measure in MEASUREMENTS.items():
        want, label = EXPECTED[key]
        out.append(Assumption(key, label, measure(), want))
    return out


def manifest_rows() -> list[list[str]]:
    return [l.split("\t") for l in MANIFEST.read_text().splitlines()[1:] if l.strip()]


def check_manifest() -> list[str]:
    """Every input against its recorded checksum, and every file against the manifest."""
    failures: list[str] = []
    rows = manifest_rows()
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

    # The other direction. A file in data/raw that no row mentions is an input nothing
    # checksums, so it is named here unless something else is stated to pin it.
    listed = {r[0] for r in rows}
    for path in sorted(RAW.iterdir()):
        if path.name in NOT_INPUTS or path.name in listed:
            continue
        pin = PINNED_ELSEWHERE.get(path.name)
        if pin is None:
            failures.append(f"{path.name}: in data/raw, pinned by nothing. Add it to "
                            f"MANIFEST.tsv or to PINNED_ELSEWHERE with what pins it.")
        elif not (RAW.parents[1] / pin).exists():
            failures.append(f"{path.name}: said to be pinned by {pin}, which is missing")
        else:
            print(f"  ok   {path.name} is pinned by {pin}")
    return failures


def main() -> int:
    if not MANIFEST.exists():
        print("MANIFEST.tsv missing, run `make inputs` first", file=sys.stderr)
        return 1

    failures = check_manifest()

    for a in assumptions():
        print(f"  {'ok  ' if a.ok else 'FAIL'} {a.label}"
              + ("" if a.ok else f" (measured {a.got!r}, expected {a.want!r})"))
        if not a.ok:
            failures.append(f"assumption broken: {a.label} "
                            f"(measured {a.got!r}, expected {a.want!r})")

    if failures:
        print("\n" + "\n".join(failures), file=sys.stderr)
        return 1
    print("\nall inputs and assumptions verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
