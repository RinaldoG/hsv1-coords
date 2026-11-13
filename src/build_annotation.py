#!/usr/bin/env python3
"""Build the canonical versioned annotation (roadmap M2).

The point of this step is not to arbitrate between annotations. `FINDINGS.md` established
that the two records do not conflict: `NC_001806.2` describes whole genes, and
`BK012101.1` is the 2020 re-annotation deposited as third-party annotation on the same
sequence, describing nested reading frames and multi-gene transcript units as well. The
counts differ because the scopes differ.

So this build does the thing that was actually missing: it labels every feature with the
CLASS of thing it is, and with where it came from, so that two annotations are never again
compared as though they measured the same quantity.

Output:
  data/derived/annotation_v0.1.0.tsv  one row per feature, classed and sourced
  data/derived/annotation_summary.md  counts per class and per source

Run: python3 -I src/build_annotation.py
"""
from __future__ import annotations

import collections
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
VERSION = "0.1.0"

sys.path.insert(0, str(ROOT / "src"))
from extract_names import parse_genbank, parse_gtf, squeeze  # noqa: E402

SOURCES = [
    ("NC_001806.2", "ref_NC_001806.2.gb", "genbank",
     "RefSeq reference annotation, provisional, last modified 13 August 2018"),
    ("BK012101.1", "alt_BK012101.1.gb", "genbank",
     "Whisnant et al. 2020 re-annotation, third-party annotation on JN555585.1, "
     "Nat Commun 11:2038, PMID 32341360"),
    ("whisnant2020-gtf", "HSV1_Annotation.gtf", "gtf",
     "Whisnant et al. 2020 annotation as distributed on Zenodo, 10.5281/zenodo.3724023"),
]

# How a feature's name and qualifiers decide what CLASS of thing it is. Order matters:
# the first match wins, so the more specific patterns come first.
CLASS_RULES: list[tuple[str, str, str]] = [
    ("sORF", r"\bsORF\b",
     "A small open reading frame named by the 2020 re-annotation"),
    ("uORF", r"\buORF\b",
     "An upstream open reading frame, 5-prime of a main coding sequence"),
    ("uoORF", r"\buoORF\b",
     "An upstream open reading frame overlapping the main coding sequence"),
    ("iORF", r"\biORF\b",
     "An internal open reading frame, inside another coding sequence"),
    ("dORF", r"\bdORF\b",
     "A downstream open reading frame, 3-prime of a main coding sequence"),
    ("transcript_span", r"^[A-Za-z0-9.]+-[A-Za-z0-9.]+$",
     "A transcript unit spanning more than one named gene, so not a gene itself"),
    ("novel_orf_2020", r"^(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+\.\d+$",
     "An open reading frame resolved by the 2020 re-annotation and numbered by "
     "position, with no independent name in the literature"),
    ("primary_mirna", r"pri-?miRNA",
     "A primary microRNA transcript"),
    ("named_gene", r"^(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+$",
     "A classically named gene in the systematic positional scheme"),
]
DEFAULT_CLASS = ("other", "Not matched by any class rule; see the notes column")


def classify(name: str) -> str:
    for label, pattern, _ in CLASS_RULES:
        if re.search(pattern, name):
            return label
    return DEFAULT_CLASS[0]


def main() -> int:
    DERIVED.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []

    for source, fname, kind, provenance in SOURCES:
        path = RAW / fname
        if not path.exists():
            print(f"missing input {path}, run `make inputs` first", file=sys.stderr)
            return 1
        feats = parse_genbank(path) if kind == "genbank" else parse_gtf(path)
        for f in feats:
            name = (f.get("gene") or "").strip()
            if not name:
                continue
            rows.append({
                "name": name,
                "feature_type": f["_type"],
                "location": f["_location"],
                "feature_class": classify(name),
                "source": source,
                "provenance": provenance,
                "product": squeeze(f.get("product") or ""),
                "notes": squeeze(f.get("note") or ""),
            })
        print(f"  {source}: {sum(1 for r in rows if r['source'] == source)} named features")

    out = DERIVED / f"annotation_v{VERSION}.tsv"
    cols = ["name", "feature_class", "feature_type", "location", "source", "product",
            "notes", "provenance"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["source"], r["name"], r["feature_type"])):
            fh.write("\t".join(r[c] for c in cols) + "\n")

    # The summary is the artefact that makes the point: counts are only comparable
    # within a class.
    per_source_class = collections.Counter((r["source"], r["feature_class"]) for r in rows)
    classes = sorted({r["feature_class"] for r in rows})
    summary = DERIVED / "annotation_summary.md"
    with summary.open("w") as fh:
        fh.write(f"# Annotation v{VERSION}: features by class and source\n\n")
        fh.write("Produced by `src/build_annotation.py`. Every feature carries the class "
                 "of thing it is and the record it came from.\n\n")
        fh.write("**Why this table exists.** The two deposited annotations of this genome "
                 "are often described as disagreeing, with 77 coding sequences against "
                 "300 on the identical 152,222 bases. They do not disagree. They describe "
                 "different classes of feature, and a count is only comparable within a "
                 "class. See `FINDINGS.md`.\n\n")
        fh.write("| Class | " + " | ".join(s for s, _, _, _ in SOURCES) + " | What it is |\n")
        fh.write("|---" * (len(SOURCES) + 2) + "|\n")
        descriptions = {label: desc for label, _, desc in CLASS_RULES}
        descriptions[DEFAULT_CLASS[0]] = DEFAULT_CLASS[1]
        for cls in classes:
            counts = [str(per_source_class[(s, cls)]) for s, _, _, _ in SOURCES]
            fh.write(f"| `{cls}` | " + " | ".join(counts) + f" | {descriptions[cls]} |\n")
        fh.write("\n## Provenance\n\n")
        for source, _, _, provenance in SOURCES:
            fh.write(f"- **`{source}`**: {provenance}\n")

    print(f"\n{len(rows)} feature rows, {len(classes)} classes")
    print(f"wrote {out.name}, {summary.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
