#!/usr/bin/env python3
"""Turn the prose test-case document into a machine-readable table.

`data/derived/alias_testcases.md` is written for a human reader. The M1 gate has to run
on data, so this extracts one assertion per (paper, gene, name string) into
`tests/testcases_synonymy.tsv`. Regenerating it is part of the build, so the prose and
the tests cannot drift apart.

Run: python3 -I src/build_testcases.py
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "derived" / "alias_testcases.md"
OUT = ROOT / "tests" / "testcases_synonymy.tsv"

ROW = re.compile(r"^\|\s*([A-Za-z0-9_.\-/() ]+?)\s*\|\s*(.+?)\s*\|")

# The prose tables use a first column that is sometimes a gene and sometimes a grouping
# label for a list of genes ("others", "miRNAs"). A grouping label is not a canonical
# name, so those rows assert only that each name resolves to something, with no expected
# target. Recorded with an empty `expected_canonical` rather than dropped, because the
# names in them are still real published usage.
GROUPING = {"others", "mirnas", "other", "mirna", "various",
            "ie transcripts", "sequence"}


def main() -> int:
    if not SRC.exists():
        print(f"missing {SRC}", file=sys.stderr)
        return 1
    blocks = re.split(r"^## (Test case .*)$", SRC.read_text(), flags=re.M)
    rows: list[tuple[str, str, str, str, str]] = []
    for i in range(1, len(blocks), 2):
        case = blocks[i].split(":")[0].strip()
        body = blocks[i + 1]
        m = re.search(r"PMID (\d+)", body)
        pmid = m.group(1) if m else ""
        if "Access: full text read" in body:
            access = "full text read"
        elif "Access: abstract only" in body:
            access = "abstract only"
        else:
            access = ""
        for line in body.splitlines():
            r = ROW.match(line)
            if not r:
                continue
            gene, names = r.group(1), r.group(2)
            if gene.lower() in ("gene", "---") or set(gene) <= set("-"):
                continue
            expected = "" if gene.lower() in GROUPING else gene
            for name in re.findall(r"`([^`]+)`", names):
                rows.append((case, pmid, access, expected, name))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("case\tpmid\taccess\texpected_canonical\tname_string\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    papers = len({r[0] for r in rows})
    print(f"{len(rows)} assertions across {papers} papers -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
