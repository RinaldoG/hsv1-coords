#!/usr/bin/env python3
"""Build the repeat and diploid-locus map (roadmap M4).

The HSV-1 genome carries two large inverted repeats, so a substantial part of it exists
in two copies. A coordinate inside one of those copies is not fully determined: it names
a position in the sequence, but not which copy of the locus is meant. The field's usual
response is to mask the repeats with Ns, which discards the information rather than
describing it.

This builds the description instead:

  data/derived/repeat_map.tsv      every repeated or structurally ambiguous region
  data/derived/diploid_pairs.tsv   loci present in both copies, paired to each other
  data/derived/repeat_map.md       the summary

Pairs are derived from the annotation by matching equal-length features across the two
copies of a repeat, not hand-listed, so the pairing can be checked and cannot quietly
drift from the annotation.

Run: python3 -I src/build_repeat_map.py
"""
from __future__ import annotations

import collections
import csv
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from build_crosswalk import genbank_sequence  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
ANNOTATION = DERIVED / "annotation_v0.1.0.tsv"
REFERENCE = "NC_001806.2"
GENOME_BP = 152222

# The two inverted repeat pairs, as the reference record's own feature table gives them.
# TRL/IRL flank the long segment; IRS/TRS flank the short one. Each pair is one locus
# set present twice, and that is what makes a coordinate inside them under-determined.
REPEAT_PAIRS = [("TRL", "IRL"), ("IRS", "TRS")]

SPAN = re.compile(r"(\d+)\.\.>?(\d+)")


def spans_of(location: str) -> list[tuple[int, int]]:
    return [(int(a), int(b)) for a, b in SPAN.findall(location)]


def extent(location: str) -> tuple[int, int] | None:
    s = spans_of(location)
    if not s:
        return None
    return min(a for a, _ in s), max(b for _, b in s)


def load_annotation() -> list[dict[str, str]]:
    with ANNOTATION.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["source"] == REFERENCE]


def main() -> int:
    if not ANNOTATION.exists():
        print(f"missing {ANNOTATION}, run `make annotation` first", file=sys.stderr)
        return 1
    rows = load_annotation()

    # 1. Every repeated or structurally ambiguous region.
    regions: list[dict[str, object]] = []
    for r in rows:
        # A repeat that carries a /gene qualifier is classed `named_gene` upstream, so
        # filtering on class alone silently dropped three of the reference's 31
        # repeat_region features. Select on the source feature type as well.
        is_repeat = (r["feature_class"] in ("repeat_region", "stem_loop", "rep_origin")
                     or r["feature_type"] in ("repeat_region", "stem_loop", "rep_origin"))
        if not is_repeat:
            continue
        ext = extent(r["location"])
        if not ext:
            continue
        # The `/gene` qualifier wins in M2's name column, so three reiteration sets
        # arrived here called `RL1`, `RL2` and `US1`: those are ICP34.5, ICP0 and
        # ICP22, whole genes of up to 3.6 kb, and this row is a 53 to 159 bp repeat
        # inside one. Round 8: joining this column against the gene dictionary
        # therefore returned a 53 bp interval for ICP34.5. The repeat's own note names
        # it, so where a note exists for a repeat feature it is preferred, and the gene
        # qualifier is kept beside it in `within_gene`.
        name = r["name"].split(";")[0].strip()
        note_name = r["notes"].split(";")[0].strip()
        within_gene = ""
        if note_name and note_name != name and (
                "reiteration" in note_name or "sequence" in note_name
                or "palindrome" in note_name):
            within_gene, name = name, note_name
        regions.append({
            "name": name,
            "within_gene": within_gene,
            "kind": (r["feature_type"]
                     if r["feature_type"] in ("repeat_region", "stem_loop", "rep_origin")
                     else r["feature_class"]),
            "start": ext[0],
            "end": ext[1],
            "length": ext[1] - ext[0] + 1,
            "notes": r["notes"][:160],
        })
    regions.sort(key=lambda x: (x["start"], x["end"]))

    named = {str(x["name"]): x for x in regions}
    big = {}
    for a, b in REPEAT_PAIRS:
        for label in (a, b):
            hit = next((x for x in regions if str(x["name"]).startswith(label)), None)
            if hit is None:
                print(f"the reference does not annotate {label}", file=sys.stderr)
                return 1
            big[label] = hit

    # 2. Diploid pairs: features of equal length sitting inside the two copies of a
    #    repeat. Derived by position and length, never hand-listed.
    def inside(feature_extent: tuple[int, int], rep: dict[str, object]) -> bool:
        return feature_extent[0] >= int(rep["start"]) and feature_extent[1] <= int(rep["end"])

    def overlaps(feature_extent: tuple[int, int], rep: dict[str, object]) -> bool:
        """Touches the repeat without being contained by it.

        LAT is why this exists. Its copy in the long repeat runs 778 bp past the end of
        IRL into IRS, so containment never matches it and the locus looked unpaired
        while the roadmap claimed otherwise.
        """
        return (feature_extent[0] <= int(rep["end"])
                and feature_extent[1] >= int(rep["start"]))

    pairs: list[dict[str, object]] = []
    for a, b in REPEAT_PAIRS:
        ra, rb = big[a], big[b]
        in_a: list[tuple[str, tuple[int, int], str]] = []
        in_b: list[tuple[str, tuple[int, int], str]] = []
        for r in rows:
            ext = extent(r["location"])
            if not ext:
                continue
            # The four big repeats are the frame, not loci inside it.
            if any(r["name"].startswith(x) and r["feature_type"] == "repeat_region"
                   for pair in REPEAT_PAIRS for x in pair):
                continue
            if inside(ext, ra) or overlaps(ext, ra):
                in_a.append((r["name"], ext, r["feature_type"]))
            if inside(ext, rb) or overlaps(ext, rb):
                in_b.append((r["name"], ext, r["feature_type"]))
        # Match on equal length and equal feature type. Where several candidates have
        # the same length, match by offset within the repeat, which is what a true
        # second copy preserves.
        used: set[int] = set()
        for name_a, ext_a, type_a in in_a:
            len_a = ext_a[1] - ext_a[0]
            off_a = ext_a[0] - int(ra["start"])
            best, best_delta = None, None
            for i, (name_b, ext_b, type_b) in enumerate(in_b):
                if i in used or type_b != type_a:
                    continue
                if ext_b[1] - ext_b[0] != len_a:
                    continue
                # IRL is the reverse copy of TRL, so the offset runs from its far end.
                off_b = int(rb["end"]) - ext_b[1]
                delta = abs(off_b - off_a)
                if best_delta is None or delta < best_delta:
                    best, best_delta = i, delta
            if best is None:
                # Same gene name in both copies but different lengths. LAT is the case
                # that matters: one copy is truncated by the end of the linear record,
                # so equal-length matching never finds it, and the roadmap claimed LAT
                # was paired when it was not. Recorded as a partial pair, with the
                # length difference stated rather than hidden.
                for i, (name_b, ext_b, type_b) in enumerate(in_b):
                    if i in used or type_b != type_a or name_b != name_a:
                        continue
                    if ext_b == ext_a:
                        continue      # the same feature seen from both repeats
                    used.add(i)
                    pairs.append({
                        "locus": name_a,
                        "repeat_pair": f"{a}/{b}",
                        "copy_1_name": name_a,
                        "copy_1_start": ext_a[0],
                        "copy_1_end": ext_a[1],
                        "copy_2_name": name_b,
                        "copy_2_start": ext_b[0],
                        "copy_2_end": ext_b[1],
                        "feature_type": type_a,
                        "offset_delta": "partial: lengths differ by "
                                        f"{abs((ext_b[1] - ext_b[0]) - len_a)} bp",
                    })
                    break
            if best is not None:
                name_b, ext_b, _ = in_b[best]
                used.add(best)
                pairs.append({
                    "locus": name_a if name_a == name_b else f"{name_a} / {name_b}",
                    "repeat_pair": f"{a}/{b}",
                    "copy_1_name": name_a,
                    "copy_1_start": ext_a[0],
                    "copy_1_end": ext_a[1],
                    "copy_2_name": name_b,
                    "copy_2_start": ext_b[0],
                    "copy_2_end": ext_b[1],
                    "feature_type": type_a,
                    "offset_delta": best_delta,
                })

    # Orientation of each copy of a multi-copy region, measured against the sequence
    # rather than stated. Round 7: the 'a' sequence copies and their orientations were
    # a literal list inside `hsv1coords/repeats.py`, in a module whose docstring said
    # nothing was hand-listed, and two mutations to that list changed a hundred live
    # answers with the whole suite green. The orientations are derivable, so they are
    # derived.
    seq = genbank_sequence(RAW / f"ref_{REFERENCE}.gb").upper()
    comp = str.maketrans("ACGTN", "TGCAN")

    def piece(a: int, b: int) -> str:
        return seq[a - 1:b]

    by_name: dict[str, list[dict]] = collections.defaultdict(list)
    for x in regions:
        by_name[str(x["name"])].append(x)
    for name, copies in by_name.items():
        if len(copies) < 2:
            for x in copies:
                x["orientation"] = "single"
            continue
        first = piece(int(copies[0]["start"]), int(copies[0]["end"]))
        for x in copies:
            here = piece(int(x["start"]), int(x["end"]))
            if len(here) != len(first):
                x["orientation"] = "unequal_length"
                continue
            rc = here.translate(comp)[::-1]
            fwd = sum(1 for p, q in zip(here, first) if p == q)
            rev = sum(1 for p, q in zip(rc, first) if p == q)
            if fwd == rev:
                x["orientation"] = "undetermined"
            else:
                x["orientation"] = "forward" if fwd > rev else "reverse"

    rmap = DERIVED / "repeat_map.tsv"
    cols = ["name", "kind", "start", "end", "length", "orientation", "within_gene",
            "notes"]
    with rmap.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in regions:
            fh.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")

    dpairs = DERIVED / "diploid_pairs.tsv"
    pcols = ["locus", "repeat_pair", "feature_type", "copy_1_name", "copy_1_start",
             "copy_1_end", "copy_2_name", "copy_2_start", "copy_2_end", "offset_delta"]
    with dpairs.open("w") as fh:
        fh.write("\t".join(pcols) + "\n")
        for p in sorted(pairs, key=lambda p: (str(p["repeat_pair"]), int(p["copy_1_start"]))):
            fh.write("\t".join(str(p[c]) for c in pcols) + "\n")

    covered = set()
    for x in regions:
        covered.update(range(int(x["start"]), int(x["end"]) + 1))

    report = DERIVED / "repeat_map.md"
    with report.open("w") as fh:
        fh.write("# The repeat and diploid-locus map\n\n")
        fh.write("Built by `src/build_repeat_map.py` from the reference record's own "
                 "feature table. Nothing here is hand-listed.\n\n")
        fh.write(f"- Repeated or structurally ambiguous regions: **{len(regions)}**\n")
        fh.write(f"- Bases inside at least one of them: **{len(covered):,}** of "
                 f"{GENOME_BP:,}, which is {100 * len(covered) / GENOME_BP:.1f} per cent "
                 f"of the genome\n")
        fh.write(f"- Diploid locus pairs derived: **{len(pairs)}**\n\n")
        fh.write("## The four inverted repeats\n\n")
        fh.write("| Repeat | Start | End | Length |\n|---|---|---|---|\n")
        for a, b in REPEAT_PAIRS:
            for label in (a, b):
                x = big[label]
                fh.write(f"| `{label}` | {int(x['start']):,} | {int(x['end']):,} | "
                         f"{int(x['length']):,} |\n")
        fh.write("\n`TRL` and `IRL` are the same 9,213 bp of sequence in two copies, and "
                 "`IRS` and `TRS` likewise at 6,633 bp. A position inside any of them "
                 "names a place in the linear record but not which copy of the locus is "
                 "meant.\n\n")
        if pairs:
            fh.write("## Diploid loci\n\n")
            fh.write("Features present in both copies of a repeat, matched by feature "
                     "type, equal length and offset within the repeat.\n\n")
            fh.write("The `Pairing` column is the one to read before transferring a "
                     "coordinate. `exact` means the two copies have the same length and "
                     "the same offset within their repeat, so a position in one maps to "
                     "the other by the repeat's own arithmetic. `partial` means the same "
                     "locus is named in both copies and the two differ in length, which "
                     "happens where a copy is truncated by the end of the linear record "
                     "or runs past the repeat: the stated length difference is how far "
                     "the two disagree, and a constant offset between them is wrong by "
                     "that much. Round 9 found this table printing the seven partial "
                     "pairs in the same form as the exact ones, on LAT, which is the "
                     "locus this repository uses as its worked example.\n\n")
            fh.write("| Locus | Repeat pair | Type | Copy 1 | Copy 2 | Pairing |\n")
            fh.write("|---|---|---|---|---|---|\n")
            for p in sorted(pairs, key=lambda p: int(p["copy_1_start"]))[:40]:
                len_1 = int(p["copy_1_end"]) - int(p["copy_1_start"]) + 1
                len_2 = int(p["copy_2_end"]) - int(p["copy_2_start"]) + 1
                if len_1 == len_2:
                    pairing = f"exact, offset delta {int(p['offset_delta'])}"
                else:
                    pairing = (f"**partial**, lengths differ by "
                               f"{abs(len_1 - len_2):,} bp")
                fh.write(f"| `{p['locus']}` | {p['repeat_pair']} | {p['feature_type']} | "
                         f"{int(p['copy_1_start']):,}..{int(p['copy_1_end']):,} | "
                         f"{int(p['copy_2_start']):,}..{int(p['copy_2_end']):,} | "
                         f"{pairing} |\n")
        kinds = collections.Counter(str(x["kind"]) for x in regions)
        fh.write("\n## Regions by kind\n\n| Kind | Count |\n|---|---|\n")
        for k, n in kinds.most_common():
            fh.write(f"| `{k}` | {n} |\n")

    print(f"  {len(regions)} regions, {len(covered):,} bp covered "
          f"({100 * len(covered) / GENOME_BP:.1f}% of the genome)")
    print(f"  {len(pairs)} diploid pairs")
    print(f"  wrote {rmap.name}, {dpairs.name}, {report.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
