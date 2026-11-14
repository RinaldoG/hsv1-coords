#!/usr/bin/env python3
"""Fetch the complete HSV-1 genomes and build a strain inventory (roadmap M3 input).

The crosswalk needs the population, not one reference. This fetches every complete HSV-1
genome in the expected length range and records, per record, its accession, strain label,
length, and the metadata a reader needs to judge whether two records claiming the same
strain actually agree.

The output inventory is the evidence for the project's motivating example: several
deposited genomes labelled with the same laboratory strain differ in length, from
different groups using different assemblers, with no mapping between them.

Outputs:
  data/raw/genomes.fasta            the sequences
  data/derived/strain_inventory.tsv one row per record

Run: python3 -I src/fetch_genomes.py
"""
from __future__ import annotations

import argparse
import hashlib
import collections
import datetime as dt
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "hsv1-coords/0.1 (data fetch script)"

QUERY = ('"Human alphaherpesvirus 1"[Organism] AND complete genome[Title] '
         'AND 140000:160000[SLEN]')
BATCH = 25
# The accession list is pinned and committed. A live NCBI query returns a different set
# as new genomes are deposited, which would make every downstream artefact
# irreproducible and would silently change the crosswalk. `--refresh` re-queries and
# rewrites the list, which is a deliberate act that shows up as a diff.
ACCESSION_LIST = pathlib.Path(__file__).resolve().parents[1] / "data" / "accessions.txt"


def get(url: str, tries: int = 3) -> bytes:
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001 - retried, then raised
            if attempt == tries:
                raise
            print(f"    attempt {attempt} failed ({exc}), retrying", file=sys.stderr)
            time.sleep(3 * attempt)
    raise RuntimeError("unreachable")


def accessions() -> list[str]:
    url = (f"{EUTILS}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(QUERY)}"
           f"&retmax=1000&retmode=json")
    import json
    data = json.loads(get(url))
    ids = data["esearchresult"]["idlist"]
    print(f"  {data['esearchresult']['count']} records match; {len(ids)} ids retrieved")
    return ids


def fetch_fasta(ids: list[str], dest: pathlib.Path) -> None:
    chunks: list[str] = []
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(batch)}"
               f"&rettype=fasta&retmode=text")
        chunks.append(get(url).decode("utf8", "replace"))
        print(f"    fetched {min(i + BATCH, len(ids))}/{len(ids)}")
        time.sleep(0.4)
    dest.write_text("".join(chunks))


def parse_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    header, seq = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                out.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        out.append((header, "".join(seq)))
    return out


# Strain labels contain characters a plain word class drops: `B^3x1.1A` became `B`,
# which then looked like four records of one strain disagreeing on length.
STRAIN = re.compile(
    r"\bstrain\s+(.+?)(?:,|\s+(?:nearly\s+)?(?:complete|partial)\b|"
    r"\s+genome\b|\s+DNA\b|\)|$)", re.I)
r"""The strain label in a FASTA defline.

Round 11: the old anchor was `(?:,|\s+complete\b|$)`, which fails on "nearly complete
genome", so `AB618031.1`'s label became its whole 69-character description and
`X14112.1`'s became "complete genome". 41 of the 183 headers were mis-parsed.
"""

QUALIFIERS = re.compile(r"^(?:isolate|clone|mutant|variant|sub-?strain|str\.?)\s+",
                        re.I)
TRAILING = re.compile(r"\s+(?:genome\s+assembly|genome|assembly|DNA|sequence)$", re.I)


def strain_group(label: str) -> str:
    """The laboratory strain a label belongs to, for counting.

    Round 11: the four KOS-family records are labelled `KOS`, `KOS`, `isolate KOS` and
    `isolate KOS 1.1`, so grouping on the label split them three ways and the inventory
    reported one label whose records disagree on length. The project's motivating claim
    is about one laboratory strain deposited more than once at different lengths, and
    this is the unit that claim is about: leading qualifiers are dropped and the first
    token is taken, so all four group as KOS while the full label is kept beside it.

    A derivative such as KOS 1.1 is a real distinction, which is why the label stays in
    its own column and only the grouping is coarsened.
    """
    bare = TRAILING.sub("", label.strip()).strip()
    # Repeatedly, not once: "Mutant  clone McKrae-dsncRNA2" has two qualifiers, and
    # stripping one left "clone" as the group, so three unrelated mutant clones were
    # grouped together under CLONE and reported as one strain of disagreeing length.
    while True:
        stripped = QUALIFIERS.sub("", bare).strip()
        if stripped == bare:
            break
        bare = stripped
    token = bare.split()[0] if bare.split() else bare
    return token.strip(",;:()").upper() or "(UNLABELLED)"


def strain_of(header: str) -> str:
    m = STRAIN.search(header)
    if m:
        label = TRAILING.sub("", m.group(1).strip().rstrip(",")).strip()
        if label:
            return label
    # Many headers read "Human alphaherpesvirus 1 <strain>, complete genome"
    body = header.split(" ", 1)[1] if " " in header else header
    for organism in ("Human alphaherpesvirus 1", "Human herpesvirus 1",
                     "Herpes simplex virus", "(type 1", "type 1"):
        body = body.replace(organism, "")
    body = re.split(r",\s*(?:nearly\s+)?(?:complete|partial)\b", body)[0]
    body = TRAILING.sub("", body.strip(" ,/()")).strip(" ,/()")
    # "complete genome" alone is not a strain label; it is what is left when the defline
    # names no strain at all. Round 11 found it recorded as X14112.1's strain.
    if not body or re.fullmatch(r"(?:nearly\s+)?(?:complete|partial)?\s*(?:genome|"
                                r"sequence|DNA)?", body, re.I):
        return "(unlabelled)"
    return body


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the FASTA is present")
    ap.add_argument("--refresh", action="store_true",
                    help="re-query NCBI and rewrite the pinned accession list")
    args = ap.parse_args()

    RAW.mkdir(parents=True, exist_ok=True)
    DERIVED.mkdir(parents=True, exist_ok=True)
    fasta = RAW / "genomes.fasta"

    if args.refresh or not ACCESSION_LIST.exists():
        ids = accessions()
        ACCESSION_LIST.parent.mkdir(parents=True, exist_ok=True)
        ACCESSION_LIST.write_text("\n".join(ids) + "\n")
        print(f"  pinned {len(ids)} accessions to {ACCESSION_LIST.name}")
    else:
        ids = [x.strip() for x in ACCESSION_LIST.read_text().split() if x.strip()]
        print(f"  using {len(ids)} pinned accessions from {ACCESSION_LIST.name}")

    if fasta.exists() and not args.force and not args.refresh:
        print(f"  have {fasta.name}")
    else:
        fetch_fasta(ids, fasta)

    records = parse_fasta(fasta)
    print(f"\n  {len(records)} sequences parsed")

    # Round 11: 133 of the 183 records were deleted from the FASTA and this script
    # printed "183 pinned accessions" and "50 sequences parsed" two lines apart, exited
    # 0, and rewrote strain_inventory.md so the project's own motivating figure went
    # from 1 to 0. `check_inputs.py` passed too, because the FASTA has no manifest row.
    # The count the script already knows is compared now, and the checksum is recorded.
    if len(records) != len(ids):
        print(f"the pinned list holds {len(ids)} accessions and the FASTA holds "
              f"{len(records)} sequences. Re-fetch with `make genomes FORCE=1` or "
              f"`python3 -I src/fetch_genomes.py --force`; a partial FASTA silently "
              f"shrinks every artefact built from it.", file=sys.stderr)
        return 1

    # Written to data/, not data/raw/, because data/raw is gitignored: a digest nobody
    # can see after a clone pins nothing.
    digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
    (ACCESSION_LIST.parent / "genomes.sha256").write_text(
        f"{digest}\t{len(records)}\tgenomes.fasta\n")
    print(f"  genomes.fasta: {len(records)} sequences, sha256 {digest[:16]}...")

    rows = []
    for header, seq in records:
        acc = header.split()[0]
        label = strain_of(header)
        rows.append({
            "accession": acc,
            "strain": label,
            "strain_group": strain_group(label),
            "length_bp": str(len(seq)),
            "ambiguous_bases": str(sum(1 for c in seq.upper() if c not in "ACGT")),
            "gc_percent": f"{100 * sum(1 for c in seq.upper() if c in 'GC') / len(seq):.2f}",
            "description": header,
        })

    out = DERIVED / "strain_inventory.tsv"
    cols = ["accession", "strain", "strain_group", "length_bp", "ambiguous_bases",
            "gc_percent", "description"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["strain"].lower(), r["accession"])):
            fh.write("\t".join(r[c] for c in cols) + "\n")

    # The motivating evidence: strains deposited more than once, and whether they agree.
    # Grouped by laboratory strain rather than by the exact label. Round 11: grouping on
    # the label split the four KOS-family records across `KOS`, `isolate KOS` and
    # `isolate KOS 1.1`, so this reported one disputed label and missed the 987 bp spread
    # that is the example the project leads with.
    by_strain: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for r in rows:
        by_strain[r["strain_group"]].append(r)
    disputed = {s: rs for s, rs in by_strain.items()
                if len(rs) > 1 and len({r["length_bp"] for r in rs}) > 1
                and s != "(UNLABELLED)"}
    by_label: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for r in rows:
        by_label[r["strain"]].append(r)

    report = DERIVED / "strain_inventory.md"
    with report.open("w") as fh:
        fh.write("# Complete HSV-1 genomes: strain inventory\n\n")
        fh.write(f"Fetched {dt.date.today().isoformat()} by `src/fetch_genomes.py`, "
                 f"NCBI nuccore query:\n\n```\n{QUERY}\n```\n\n")
        fh.write(f"- Records: **{len(rows)}**\n")
        fh.write(f"- Distinct strain labels, as deposited: **{len(by_label)}**\n")
        fh.write(f"- Distinct laboratory strains, after dropping qualifiers such as "
                 f"`isolate` and `clone`: **{len(by_strain)}**\n")
        fh.write(f"- Laboratory strains deposited more than once: "
                 f"**{sum(1 for s, rs in by_strain.items() if len(rs) > 1 and s != '(UNLABELLED)')}**\n")
        fh.write(f"- Of those, strains whose records **disagree on length**: "
                 f"**{len(disputed)}**\n\n")
        fh.write("The two counts differ because a depositor may write `KOS`, `isolate "
                 "KOS` and `isolate KOS 1.1` for records of the same laboratory strain. "
                 "Grouping on the label as deposited splits those three ways, which is "
                 "what this file did until round 11, and it reported one disputed strain "
                 "while missing the KOS family entirely. The `strain` column in the TSV "
                 "keeps the label as deposited; `strain_group` is the unit counted "
                 "here.\n\n"
                 "The grouping drops leading qualifiers only. A hyphenated construct "
                 "such as `McKrae-dsncRNA2` keeps its own identity rather than joining "
                 "`McKrae`, because merging constructs with their parents would also "
                 "merge isolate designations such as `HSV1/USA/WA-UW-2L9/2020`, which "
                 "are not laboratory strains. So these counts are a floor on how much "
                 "the same material has been deposited at different lengths, not a "
                 "complete reckoning.\n\n")
        lo = min(int(r["length_bp"]) for r in rows)
        hi = max(int(r["length_bp"]) for r in rows)
        fh.write(f"- Length range: {lo:,} to {hi:,} bp, a spread of {hi - lo:,} bp\n\n")
        if disputed:
            fh.write("## Same strain label, different length\n\n")
            fh.write("These are the records a coordinate crosswalk has to reconcile. A "
                     "position cited against one of them does not mean the same position "
                     "in another.\n\n")
            for strain, rs in sorted(disputed.items(),
                                     key=lambda kv: -(max(int(r['length_bp']) for r in kv[1])
                                                      - min(int(r['length_bp']) for r in kv[1]))):
                lens = sorted(int(r["length_bp"]) for r in rs)
                fh.write(f"### `{strain}` - {len(rs)} records, spread {lens[-1] - lens[0]:,} bp\n\n")
                fh.write("| Accession | Length | Ambiguous bases |\n|---|---|---|\n")
                for r in sorted(rs, key=lambda r: -int(r["length_bp"])):
                    fh.write(f"| `{r['accession']}` | {int(r['length_bp']):,} | "
                             f"{r['ambiguous_bases']} |\n")
                fh.write("\n")

    print(f"  {len(by_strain)} strain labels, {len(disputed)} disagree on length")
    print(f"  wrote {out.name}, {report.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
