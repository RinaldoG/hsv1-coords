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
STRAIN = re.compile(r"\bstrain\s+(\S+?)(?:,|\s+complete\b|$)", re.I)


def strain_of(header: str) -> str:
    m = STRAIN.search(header)
    if m:
        return m.group(1).rstrip(",")
    # Many headers read "Human alphaherpesvirus 1 <strain>, complete genome"
    body = header.split(" ", 1)[1] if " " in header else header
    body = body.replace("Human alphaherpesvirus 1", "").replace("Human herpesvirus 1", "")
    body = body.split(", complete genome")[0].strip(" ,")
    return body or "(unlabelled)"


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

    rows = []
    for header, seq in records:
        acc = header.split()[0]
        rows.append({
            "accession": acc,
            "strain": strain_of(header),
            "length_bp": str(len(seq)),
            "ambiguous_bases": str(sum(1 for c in seq.upper() if c not in "ACGT")),
            "gc_percent": f"{100 * sum(1 for c in seq.upper() if c in 'GC') / len(seq):.2f}",
            "description": header,
        })

    out = DERIVED / "strain_inventory.tsv"
    cols = ["accession", "strain", "length_bp", "ambiguous_bases", "gc_percent",
            "description"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["strain"].lower(), r["accession"])):
            fh.write("\t".join(r[c] for c in cols) + "\n")

    # The motivating evidence: strains deposited more than once, and whether they agree.
    by_strain: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for r in rows:
        by_strain[r["strain"]].append(r)
    disputed = {s: rs for s, rs in by_strain.items()
                if len(rs) > 1 and len({r["length_bp"] for r in rs}) > 1}

    report = DERIVED / "strain_inventory.md"
    with report.open("w") as fh:
        fh.write("# Complete HSV-1 genomes: strain inventory\n\n")
        fh.write(f"Fetched {dt.date.today().isoformat()} by `src/fetch_genomes.py`, "
                 f"NCBI nuccore query:\n\n```\n{QUERY}\n```\n\n")
        fh.write(f"- Records: **{len(rows)}**\n")
        fh.write(f"- Distinct strain labels: **{len(by_strain)}**\n")
        fh.write(f"- Strain labels deposited more than once: "
                 f"**{sum(1 for rs in by_strain.values() if len(rs) > 1)}**\n")
        fh.write(f"- Of those, labels whose records **disagree on length**: "
                 f"**{len(disputed)}**\n\n")
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
