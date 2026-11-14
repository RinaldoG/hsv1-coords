#!/usr/bin/env python3
"""Extract gene names from the deposited records and reconcile them across sources.

This is the mechanical half of milestone 1. It reads what the three records and the
Whisnant 2020 annotation actually say, and produces:

  data/derived/names_by_source.tsv   one row per (source, gene), with products and notes
  data/derived/name_reconciliation.tsv  one row per gene name, present in which sources
  data/derived/mechanical_report.md  what the reconciliation shows

It does not invent synonymy. A name appearing in two sources is recorded as appearing in
two sources, which is a fact. Whether two different names mean the same gene is a
literature question, handled separately and with a citation per alias.

Run: python3 -I src/extract_names.py
"""
from __future__ import annotations

import collections
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"

# source label -> (filename, kind)
SOURCES = {
    "NC_001806.2": ("ref_NC_001806.2.gb", "genbank"),
    "BK012101.1": ("alt_BK012101.1.gb", "genbank"),
    "JN555585.1": ("base_JN555585.1.gb", "genbank"),
    "whisnant2020": ("HSV1_Annotation.gtf", "gtf"),
}

FEATURE_START = re.compile(r"^     (\S+)\s+(\S.*)$")
QUALIFIER = re.compile(r"^                     /(\w+)=(.*)$")


def parse_genbank(path: pathlib.Path) -> list[dict]:
    """Return one dict per feature, with its qualifiers folded to single strings."""
    features: list[dict] = []
    current: dict | None = None
    key: str | None = None
    in_features = False

    for line in path.read_text().splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN"):
            break

        m = FEATURE_START.match(line)
        if m:
            if current:
                features.append(current)
            current = {"_type": m.group(1), "_location": m.group(2)}
            # A location can wrap onto the next line, and GenBank wraps it with no
            # space. Setting the continuation key here is what stops those lines being
            # dropped; leaving it None silently truncated every wrapped location,
            # including RL2 (ICP0) and the IRL2-IRL1 span.
            key = "_location"
            continue
        if current is None:
            continue

        q = QUALIFIER.match(line)
        if q:
            key = q.group(1)
            val = q.group(2).strip().strip('"')
            current[key] = val
        elif key is not None and line.startswith(" " * 21):
            # A wrapped continuation line. Locations wrap with no separator; qualifier
            # text wraps at a space, so the two are joined differently.
            joiner = "" if key == "_location" else " "
            current[key] = (current[key] + joiner + line.strip().strip('"')).strip()

    if current:
        features.append(current)
    return features


def parse_gtf(path: pathlib.Path) -> list[dict]:
    out: list[dict] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', parts[8]))
        # Column 7 is the strand, and it was not read. Round 10: 423 of this file's
        # lines are on the minus strand and every one of them reached the annotation
        # table as a plain `start..end` range, while the same features from the
        # GenBank-format version of the same annotation carry `complement(...)`. The
        # absence of `complement` therefore read as a positive assertion of the plus
        # strand for 83 genes including UL23, UL27, UL29, UL48 and IRS1, and anyone
        # pulling codons for them would have got the reverse complement of a CDS.
        span = f"{parts[3]}..{parts[4]}"
        location = f"complement({span})" if parts[6] == "-" else span
        out.append({"_type": parts[2], "_location": location,
                    "strand": parts[6],
                    "gene": attrs.get("gene_id", ""),
                    # `gene_type` is a biotype, not a product: it is the constant
                    # "protein_coding" on all 910 lines, LAT and the pri-miRNAs
                    # included. It was being written into the `product` column, where a
                    # reader would take it for a product.
                    "biotype": attrs.get("gene_type", ""),
                    # The functional annotation this file does carry, which was being
                    # discarded: `comment "(ICP4; IE175;)"`, `comment "(=pri-miRNA-H1)"`.
                    "note": attrs.get("comment", "").strip("()"),
                    "transcript": attrs.get("transcript_id", "")})
    return out


def squeeze(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def main() -> int:
    DERIVED.mkdir(parents=True, exist_ok=True)

    # gene name -> source -> {products, notes, types, locus_tags}
    index: dict[str, dict[str, dict[str, set[str]]]] = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(set)))

    for label, (fname, kind) in SOURCES.items():
        path = RAW / fname
        if not path.exists():
            print(f"missing input {path}, run `make inputs` first", file=sys.stderr)
            return 1
        feats = parse_genbank(path) if kind == "genbank" else parse_gtf(path)
        for f in feats:
            gene = (f.get("gene") or "").strip()
            if not gene:
                continue
            slot = index[gene][label]
            slot["types"].add(f["_type"])
            for field in ("product", "note", "locus_tag", "acronym", "function"):
                v = f.get(field)
                if v:
                    slot[field].add(squeeze(v))
        print(f"  {label}: {len(feats)} features, "
              f"{sum(1 for g in index if label in index[g])} named genes so far")

    # names_by_source.tsv
    by_source = DERIVED / "names_by_source.tsv"
    with by_source.open("w") as fh:
        fh.write("gene\tsource\tfeature_types\tproducts\tlocus_tags\tnotes\n")
        for gene in sorted(index):
            for label in SOURCES:
                if label not in index[gene]:
                    continue
                s = index[gene][label]
                fh.write("\t".join([
                    gene, label,
                    ";".join(sorted(s["types"])),
                    " | ".join(sorted(s["product"])) or "",
                    ";".join(sorted(s["locus_tag"])) or "",
                    " | ".join(sorted(s["note"])) or "",
                ]) + "\n")

    # name_reconciliation.tsv
    recon = DERIVED / "name_reconciliation.tsv"
    counts = collections.Counter()
    with recon.open("w") as fh:
        fh.write("gene\t" + "\t".join(SOURCES) + "\tn_sources\tproduct_agreement\n")
        for gene in sorted(index):
            present = [("yes" if label in index[gene] else "no") for label in SOURCES]
            n = present.count("yes")
            prods = {label: " | ".join(sorted(index[gene][label]["product"]))
                     for label in SOURCES if label in index[gene]
                     and index[gene][label]["product"]}
            distinct = set(prods.values())
            if len(distinct) <= 1:
                agree = "single_or_absent"
            else:
                agree = "DIFFERS: " + " // ".join(f"{k}={v}" for k, v in sorted(prods.items()))
            counts[n] += 1
            fh.write(f"{gene}\t" + "\t".join(present) + f"\t{n}\t{agree}\n")

    # The report
    only = {label: sorted(g for g in index if list(index[g]) == [label])
            for label in SOURCES}
    all_four = sorted(g for g in index if len(index[g]) == len(SOURCES))
    differs = sum(1 for line in recon.read_text().splitlines()[1:] if "DIFFERS" in line)

    report = DERIVED / "mechanical_report.md"
    with report.open("w") as fh:
        fh.write("# Gene names in the deposited records\n\n")
        fh.write("Produced by `src/extract_names.py` from the files in `data/raw`. Facts "
                 "only: a name appearing in two sources is recorded as appearing in two "
                 "sources. Whether two different names mean the same gene is a literature "
                 "question and is not decided here.\n\n")
        fh.write(f"- Distinct gene names across all sources: **{len(index)}**\n")
        fh.write(f"- Present in all {len(SOURCES)} sources: **{len(all_four)}**\n")
        fh.write(f"- Genes whose `product` text differs between sources: **{differs}**\n\n")
        fh.write("## Names per source\n\n| Source | Named genes | Unique to it |\n|---|---|---|\n")
        for label in SOURCES:
            n = sum(1 for g in index if label in index[g])
            fh.write(f"| `{label}` | {n} | {len(only[label])} |\n")
        fh.write("\n## How many sources each name appears in\n\n")
        fh.write("| Sources | Gene names |\n|---|---|\n")
        for k in sorted(counts, reverse=True):
            fh.write(f"| {k} | {counts[k]} |\n")
        for label in SOURCES:
            if not only[label]:
                continue
            fh.write(f"\n## Only in `{label}` ({len(only[label])})\n\n")
            fh.write(", ".join(f"`{g}`" for g in only[label]) + "\n")
    print(f"\n{len(index)} distinct gene names")
    print(f"wrote {by_source.name}, {recon.name}, {report.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
