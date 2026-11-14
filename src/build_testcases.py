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

ROW = re.compile(r"^\|\s*([A-Za-z0-9_.\-/(), ]+?)\s*\|\s*(.+?)\s*\|")
"""A row of a prose table: the gene column, then the names column.

Round 11: the gene column's character class excluded the comma, so every row naming more
than one gene was skipped without a word. Three rows and nine name strings went missing
from the generated test cases, seven of which (ICP22, gE, gI, gM, UL25, UL17, UL36)
appeared nowhere else in them, while the builder's docstring said regenerating the file is
part of the build "so the prose and the tests cannot drift apart". The comma is allowed
now, a gene column holding several genes is split on it, and any row this parser cannot
read is counted and reported rather than passed over.
"""

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
    unread: list[tuple[str, str]] = []
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
            if not line.lstrip().startswith("|"):
                continue
            r = ROW.match(line)
            if not r:
                # A table row this parser cannot read. Counted, because silently
                # skipping these is what lost nine name strings.
                first = line.strip().strip("|").split("|")[0].strip()
                # Only a row that looks like a gene row. The prose also holds citation
                # tables, whose first cell is an author list with a year, and those are
                # not test cases. Deliberately narrow: a gene column is short, carries
                # no year and no PMID.
                looks_like_a_gene = (
                    "`" in line and len(first) <= 40
                    and not re.search(r"\b(19|20)\d\d\b|PMID|et al", first))
                if looks_like_a_gene:
                    unread.append((case, line.strip()[:70]))
                continue
            gene, names = r.group(1), r.group(2)
            if gene.lower() in ("gene", "---") or set(gene) <= set("-"):
                continue
            genes = [g.strip() for g in gene.split(",") if g.strip()]
            found = re.findall(r"`([^`]+)`", names)
            if not found:
                # A row naming genes but quoting no alias, such as "| UL20, UL45 |
                # systematic only |". There is nothing to assert, which is different
                # from a row this parser could not read.
                continue
            if len(genes) > 1:
                # "US1, US12" against "`ICP22`, `ICP47`": paired in order, which is how
                # the prose reads. Any other shape is reported rather than guessed at.
                if len(genes) != len(found):
                    unread.append((case, f"{len(genes)} genes against {len(found)} "
                                         f"names: {line.strip()[:60]}"))
                    continue
                for expected_gene, name in zip(genes, found):
                    rows.append((case, pmid, access,
                                 "" if expected_gene.lower() in GROUPING
                                 else expected_gene, name))
                continue
            expected = "" if gene.lower() in GROUPING else gene
            for name in found:
                rows.append((case, pmid, access, expected, name))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("case\tpmid\taccess\texpected_canonical\tname_string\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    papers = len({r[0] for r in rows})
    print(f"{len(rows)} assertions across {papers} papers -> {OUT}")
    if unread:
        print(f"\n{len(unread)} table rows could not be read and are NOT in the test "
              f"cases:", file=sys.stderr)
        for case, line in unread:
            print(f"  {case}: {line}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
