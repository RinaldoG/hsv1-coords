"""Copy the built artefacts into the package, so a wheel carries its own data.

Round 11: `pip install .` produced an importable library in which every loader raised
FileNotFoundError, because each module resolved its data as `parents[2] / "data"`, which
is the repository layout and not the installed one. `pip install -e .` worked, which is
what README.md recommends, so the breakage hit only the user who installed this as a
dependency, which is the user least able to diagnose it.

`hsv1coords._paths` prefers `hsv1coords/data` and falls back to the checkout, so this
script is what makes the first of those exist. It is listed in `make all` and the staged
copy is gitignored: it is a build product, not a second copy of the data to keep in step.
"""
from __future__ import annotations

import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data"
TARGET = ROOT / "src" / "hsv1coords" / "data"

WANTED = ("conventions.tsv", "genomes.sha256", "accessions.txt")
DERIVED = ("aliases_literature.tsv", "annotation_v0.1.0.tsv", "repeat_map.tsv",
           "diploid_pairs.tsv", "version_crosswalk.tsv", "crosswalk_coverage.tsv",
           "base_agreement.tsv", "strain_inventory.tsv")


def main() -> int:
    if not SOURCE.is_dir():
        print(f"no {SOURCE}; run the build first", file=sys.stderr)
        return 1
    if TARGET.exists():
        shutil.rmtree(TARGET)
    (TARGET / "derived" / "crosswalk").mkdir(parents=True, exist_ok=True)
    copied = 0
    for name in WANTED:
        src = SOURCE / name
        if src.exists():
            shutil.copy2(src, TARGET / name)
            copied += 1
    for name in DERIVED:
        src = SOURCE / "derived" / name
        if src.exists():
            shutil.copy2(src, TARGET / "derived" / name)
            copied += 1
    blocks = sorted((SOURCE / "derived" / "crosswalk").glob("*.blocks.tsv"))
    for src in blocks:
        shutil.copy2(src, TARGET / "derived" / "crosswalk" / src.name)
    copied += len(blocks)
    missing = [n for n in DERIVED if not (TARGET / "derived" / n).exists()]
    print(f"staged {copied} files into {TARGET.relative_to(ROOT)} "
          f"({len(blocks)} block files)")
    if missing or not blocks:
        print(f"missing from the staged copy: {missing or 'the block files'}. "
              f"Run `make all` first, or a wheel built from this tree will import and "
              f"then fail on every call.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
