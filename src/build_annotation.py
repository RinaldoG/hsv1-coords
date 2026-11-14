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
# Features that carry no /gene qualifier but must not be dropped. The inverted repeats
# are the reason this repository exists, and the first version of this build silently
# omitted all four of them because it only kept features with a gene name.
UNNAMED_CLASSES = {
    "repeat_region": ("repeat_region",
                      "A repeat, including the inverted repeats and the reiteration "
                      "sets nested inside them"),
    "stem_loop": ("stem_loop", "A stem-loop structure, such as an origin palindrome"),
    "rep_origin": ("rep_origin", "An origin of replication"),
    "misc_feature": ("misc_feature", "A feature the source record did not type further"),
    "regulatory": ("regulatory",
                   "A regulatory element. In this reference record 53 are "
                   "polyadenylation signals and 19 are TATA boxes, which the `subtype` "
                   "column now distinguishes; the earlier description said \"such as a "
                   "TATA box\" of a set that is three quarters polyA signals"),
    "misc_RNA": ("misc_RNA", "An RNA the source record did not type further"),
    "source": (None, None),          # the whole-genome record, deliberately not a feature
}

# A feature type that is not a gene or a transcript is classed by what it IS, whatever
# its name. Round 10: `classify` took only the name, so 47 six-base regulatory elements
# carrying a `/gene` qualifier were classed `named_gene`, under a description reading "A
# classically named gene in the systematic positional scheme", while the 23 regulatory
# features that name their gene in `/note` instead went to class `regulatory`. The same
# polyadenylation signal landed in a different class depending on which qualifier the
# depositor used, and `feature_class == "named_gene"` is the obvious filter for "give me
# the genes".
GENE_LIKE_TYPES = {"gene", "CDS", "mRNA", "transcript", "exon", "ncRNA", "misc_RNA"}

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
    ("lat_region", r"^(LAT|AL-RNA|ORF-O|ORF-P|LAT-ORF|ORF-O/P)",
     "A named transcript or open reading frame in the LAT region, which is named by "
     "function rather than by position"),
    ("origin_locus", r"^Ori-[SL]",
     "A locus named after an origin of replication. Note that `Ori-S` here is the 2020 "
     "re-annotation's locus spanning the TRS copy and TRS1, not the oriS palindrome "
     "itself, which is a `stem_loop` feature"),
    ("transcript_span", r"^(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+(\.\d+)?"
                        r"-(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+(\.\d+)?$",
     "A transcript unit spanning more than one named gene, so not a gene itself. Both "
     "halves must be systematic names: the earlier rule matched any hyphenated name, "
     "so AL-RNA, ORF-O, ORF-P, Ori-S and gG-1 were all declared multi-gene spans"),
    ("novel_orf_2020", r"^(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+\.\d+$",
     "An open reading frame resolved by the 2020 re-annotation and numbered by "
     "position, with no independent name in the literature. A decimal name that the "
     "2018 reference record already carries is a named gene instead, which is what "
     "UL26.5 (ICP35, VP22a, described in 1974) is"),
    ("primary_mirna", r"pri-?miRNA",
     "A primary microRNA transcript"),
    ("named_gene", r"^(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+(\.\d+)?[A-Z]?$",
     "A classically named gene in the systematic positional scheme. The trailing letter "
     "is part of the scheme: UL49A (gN) and US8A are genes, and were classed `other`"),
]
DEFAULT_CLASS = ("other", "Not matched by any class rule; see the notes column")


def subtype_of(feature: dict) -> str:
    """The source's own sub-typing qualifier for this feature, if it has one."""
    for qualifier in SUBTYPE_QUALIFIERS:
        value = (feature.get(qualifier) or "").strip()
        if value:
            return value
    return ""


def classify(name: str, feature_type: str = "gene",
             pre_2020_names: frozenset[str] = frozenset()) -> str:
    """What class of thing a feature is, from its name, its type and its provenance.

    `feature_type` decides first for anything that is not gene-like: a regulatory
    element is a regulatory element whether or not it carries a gene name.

    `pre_2020_names` are the names the 2018 reference record already uses. A positional
    decimal name in that set was not resolved by the 2020 re-annotation, whatever its
    shape, which is the UL26.5 case.
    """
    if feature_type and feature_type not in GENE_LIKE_TYPES:
        if feature_type in UNNAMED_CLASSES and UNNAMED_CLASSES[feature_type][0]:
            return UNNAMED_CLASSES[feature_type][0]
        # A type that is neither gene-like nor named here. Round 11, as a latent case:
        # this fell through to the name rules, so a record carrying `mat_peptide`,
        # `3'UTR` or `polyA_site` with a `/gene` qualifier would be classed
        # `named_gene`, which is the defect the type-first branch was added to fix. The
        # three records in use carry no such type today, so this is reached by nothing;
        # a new source would reach it immediately.
        return DEFAULT_CLASS[0]
    for label, pattern, _ in CLASS_RULES:
        if re.search(pattern, name):
            if label == "novel_orf_2020" and name in pre_2020_names:
                continue
            return label
    return DEFAULT_CLASS[0]


SUBTYPE_QUALIFIERS = ("regulatory_class", "rpt_type", "ncRNA_class", "biotype")
"""The source's own sub-typing qualifier, kept in a column of its own.

Round 10: `regulatory_class` was dropped entirely, so TATA_box and
polyA_signal_sequence were unrecoverable from the table, and the GTF's `gene_type` was
being written into `product`, where the constant "protein_coding" read as a product for
LAT and the pri-miRNAs.
"""


def pre_2020_names(sources, raw) -> frozenset[str]:
    """The gene names the 2018 reference record already carries."""
    names = set()
    for source, filename, kind, _ in sources:
        if source != "NC_001806.2":
            continue
        path = raw / filename
        if not path.exists():
            return frozenset()
        for f in (parse_genbank(path) if kind == "genbank" else parse_gtf(path)):
            name = (f.get("gene") or "").strip()
            if name:
                names.add(name)
    return frozenset(names)


def main() -> int:
    DERIVED.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    known = pre_2020_names(SOURCES, RAW)

    for source, fname, kind, provenance in SOURCES:
        path = RAW / fname
        if not path.exists():
            print(f"missing input {path}, run `make inputs` first", file=sys.stderr)
            return 1
        feats = parse_genbank(path) if kind == "genbank" else parse_gtf(path)
        for f in feats:
            name = (f.get("gene") or "").strip()
            if not name:
                ftype = f["_type"]
                if ftype not in UNNAMED_CLASSES or UNNAMED_CLASSES[ftype][0] is None:
                    continue
                note = squeeze(f.get("note") or "")
                # Name it from its own note where it has one, so TRL and IRL are
                # findable rather than appearing as identical "repeat_region" rows.
                label = note.split(";")[0].strip() if note else ""
                rows.append({
                    "name": label or f"{ftype}@{f['_location']}",
                    "feature_type": ftype,
                    "location": f["_location"],
                    "feature_class": UNNAMED_CLASSES[ftype][0],
                    "source": source,
                    "provenance": provenance,
                    "product": squeeze(f.get("product") or ""),
                    "subtype": subtype_of(f),
                    "transcript_id": squeeze(f.get("transcript") or ""),
                    "notes": note,
                })
                continue
            rows.append({
                "name": name,
                "feature_type": f["_type"],
                "location": f["_location"],
                "feature_class": classify(name, f["_type"], known),
                "source": source,
                "provenance": provenance,
                "product": squeeze(f.get("product") or ""),
                "subtype": subtype_of(f),
                "transcript_id": squeeze(f.get("transcript") or ""),
                "notes": squeeze(f.get("note") or ""),
            })
        print(f"  {source}: {sum(1 for r in rows if r['source'] == source)} named features")

    out = DERIVED / f"annotation_v{VERSION}.tsv"
    cols = ["name", "feature_class", "feature_type", "location", "source", "product",
            "subtype", "transcript_id", "notes", "provenance"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["source"], r["name"], r["feature_type"])):
            fh.write("\t".join(r[c] for c in cols) + "\n")

    # The summary is the artefact that makes the point: counts are only comparable
    # within a class AND within a feature type. Round 10: it counted rows, and the
    # feature-type mix differs per source, so reading across a row compared different
    # quantities. It reported 30 against 479 in one class for two formats of the same
    # annotation by the same authors, because one format states a gene and an mRNA per
    # locus and the other states a transcript and its exons. Distinct names are
    # comparable, so that is what the table counts, with the row counts and the
    # feature-type mix stated separately rather than blended.
    names_by = collections.defaultdict(set)
    rows_by = collections.Counter()
    types_by = collections.defaultdict(collections.Counter)
    for r in rows:
        key = (r["source"], r["feature_class"])
        names_by[key].add(r["name"])
        rows_by[key] += 1
        types_by[key][r["feature_type"]] += 1
    classes = sorted({r["feature_class"] for r in rows})
    type_rows = collections.Counter((r["source"], r["feature_type"]) for r in rows)
    feature_types = sorted({r["feature_type"] for r in rows})
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
        fh.write("**How to read it.** The first table counts **distinct names**, which is "
                 "the only figure comparable across these sources. A row count is not: "
                 "`NC_001806.2` states a gene and a CDS and sometimes an mRNA for one "
                 "locus, `BK012101.1` states a gene and an mRNA, and the GTF states a "
                 "transcript and each of its exons, so the same locus contributes a "
                 "different number of rows to each. The second table gives the rows by "
                 "feature type, where a count is comparable because the type is fixed. "
                 "Until round 10 this file blended the two, and reported 30 against 479 "
                 "in one class for two formats of the same annotation by the same "
                 "authors.\n\n")
        fh.write("## Distinct names by class\n\n")
        fh.write("| Class | " + " | ".join(s for s, _, _, _ in SOURCES)
                 + " | What it is |\n")
        fh.write("|---" * (len(SOURCES) + 2) + "|\n")
        descriptions = {label: desc for label, _, desc in CLASS_RULES}
        descriptions[DEFAULT_CLASS[0]] = DEFAULT_CLASS[1]
        for cls, desc in UNNAMED_CLASSES.values():
            if cls:
                descriptions[cls] = desc
        for cls in classes:
            counts = [str(len(names_by[(s, cls)])) for s, _, _, _ in SOURCES]
            fh.write(f"| `{cls}` | " + " | ".join(counts)
                     + f" | {descriptions[cls]} |\n")
        fh.write("\n## Rows by feature type\n\n")
        fh.write("| Feature type | " + " | ".join(s for s, _, _, _ in SOURCES)
                 + " |\n")
        fh.write("|---" * (len(SOURCES) + 1) + "|\n")
        for ftype in feature_types:
            counts = [str(type_rows[(s, ftype)]) for s, _, _, _ in SOURCES]
            fh.write(f"| `{ftype}` | " + " | ".join(counts) + " |\n")
        fh.write("\n## Rows by class, and the feature types in them\n\n")
        fh.write("Stated so that nobody has to infer the mix from the counts above.\n\n")
        fh.write("| Class | Source | Rows | Feature types |\n|---|---|---|---|\n")
        for cls in classes:
            for source, _, _, _ in SOURCES:
                key = (source, cls)
                if not rows_by[key]:
                    continue
                mix = ", ".join(f"{t} {n}" for t, n in
                                sorted(types_by[key].items()))
                fh.write(f"| `{cls}` | `{source}` | {rows_by[key]} | {mix} |\n")
        fh.write("\n## Provenance\n\n")
        for source, _, _, provenance in SOURCES:
            fh.write(f"- **`{source}`**: {provenance}\n")

    print(f"\n{len(rows)} feature rows, {len(classes)} classes")
    print(f"wrote {out.name}, {summary.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
