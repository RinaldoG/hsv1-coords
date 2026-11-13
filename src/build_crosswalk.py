#!/usr/bin/env python3
"""Align every complete genome to the reference and build the coordinate crosswalk (M3).

Each genome is aligned to `NC_001806.2` with minimap2 (the preset in `PRESET`,
with `--cs`), and the
alignment is reduced to a list of ungapped blocks: (query_start, target_start, length).
A lift-over is then a search over those blocks.

Design decisions, stated because they constrain what the output means:

- **One target, the reference.** Every genome maps to `NC_001806.2`, so any pair of
  genomes is related through it. Direct pairwise alignment of all 183 against each other
  would be 16,653 alignments and would not obviously be better: errors would differ by
  path, and two routes between the same pair could disagree.
- **Ungapped blocks only.** A position inside an insertion relative to the reference has
  no reference coordinate, and the lift-over returns nothing rather than the nearest
  match. Returning a nearby coordinate silently is the failure this repository exists to
  prevent.
- **Secondary alignments are discarded.** The inverted repeats mean parts of the genome
  align to two places. Keeping both would produce a lift-over that silently picks one.
  Those regions are what milestone 4 marks as un-callable; until then the crosswalk simply
  declines there, which is visible in the coverage figures.

Outputs:
  data/derived/crosswalk/<accession>.blocks.tsv   per genome, the ungapped blocks
  data/derived/crosswalk_coverage.tsv             per genome, how much is mappable
  data/derived/crosswalk_report.md                the summary

Run: python3 -I src/build_crosswalk.py [--limit N]
"""
from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
CROSSWALK = DERIVED / "crosswalk"
REFERENCE = "NC_001806.2"
# `asm5` is documented for about 0.1 per cent divergence. HSV-1 strains differ by 0.7 to
# 1.3 per cent in substitutions and more with indels, so asm5 was too strict and lost
# real alignments: measured on this corpus, moving to asm10 raises the median mappable
# fraction from 98.08 to 98.84 per cent and takes the genomes below 95 per cent from 13
# to 8.
PRESET = "asm10"


def read_fasta(path: pathlib.Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif name and line.strip():
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def genbank_sequence(path: pathlib.Path) -> str:
    body = path.read_text().split("ORIGIN")[1].split("//")[0]
    return "".join(re.findall(r"[acgtnACGTN]", body.replace("\n", ""))).upper()


CS_TOKEN = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def blocks_from_cs(q_start: int, t_start: int, cs: str) -> list[tuple[int, int, int]]:
    """Reduce a cs tag to ungapped (query_start, target_start, length) blocks, 0-based.

    Only the short form is understood. `--cs=long` uses `=` and `~`, which this parser
    would silently skip, producing a crosswalk that looks fine and is wrong everywhere,
    so anything not fully consumed raises instead.
    """
    consumed = sum(len(op) + len(val) for op, val in CS_TOKEN.findall(cs))
    if consumed != len(cs):
        raise ValueError(
            f"cs tag not fully understood ({consumed} of {len(cs)} characters); "
            f"this parser handles the short form only, not --cs=long")
    q, t = q_start, t_start
    blocks: list[tuple[int, int, int]] = []
    for op, val in CS_TOKEN.findall(cs):
        if op == ":":                       # identical run of val bases
            n = int(val)
            blocks.append((q, t, n))
            q += n
            t += n
        elif op == "*":                     # substitution, one base, still aligned
            blocks.append((q, t, 1))
            q += 1
            t += 1
        elif op == "+":                     # insertion in query, no target coordinate
            q += len(val)
        elif op == "-":                     # deletion from query, no query coordinate
            t += len(val)
    # Merge blocks that are contiguous in both sequences.
    merged: list[tuple[int, int, int]] = []
    for b in blocks:
        if merged and merged[-1][0] + merged[-1][2] == b[0] and \
           merged[-1][1] + merged[-1][2] == b[1]:
            qs, ts, ln = merged[-1]
            merged[-1] = (qs, ts, ln + b[2])
        else:
            merged.append(b)
    return merged


def ambiguous_intervals(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Half-open intervals covered more than once by `spans`."""
    events: list[tuple[int, int]] = []
    for start, end in spans:
        events.append((start, 1))
        events.append((end, -1))
    events.sort()
    out: list[tuple[int, int]] = []
    depth = 0
    run_start = None
    for pos, delta in events:
        was_deep = depth > 1
        depth += delta
        now_deep = depth > 1
        if not was_deep and now_deep:
            run_start = pos
        elif was_deep and not now_deep and run_start is not None:
            out.append((run_start, pos))
            run_start = None
    return out


def subtract(block: tuple[int, int, int], bad_q: list[tuple[int, int]],
             bad_t: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    """Remove ambiguous query and target intervals from one block."""
    qs, ts, ln = block
    cuts: list[tuple[int, int]] = []
    for lo, hi in bad_q:
        cuts.append((max(qs, lo), min(qs + ln, hi)))
    for lo, hi in bad_t:                 # translate target interval into query offsets
        cuts.append((max(qs, qs + (lo - ts)), min(qs + ln, qs + (hi - ts))))
    cuts = sorted((a, b) for a, b in cuts if b > a)
    pieces: list[tuple[int, int, int]] = []
    cursor = qs
    for lo, hi in cuts:
        if lo > cursor:
            pieces.append((cursor, ts + (cursor - qs), lo - cursor))
        cursor = max(cursor, hi)
    if cursor < qs + ln:
        pieces.append((cursor, ts + (cursor - qs), qs + ln - cursor))
    return pieces


def union_length(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by these half-open intervals, counting overlap once.

    Round 7: the dropped figure was `sum(block lengths) - sum(kept lengths)`, which
    counts a doubly-placed region twice, because both placements are in the input. On
    the real data 19 of 24 nonzero rows were exactly 2x the genome bases removed. The
    figure a reader wants is bases of this genome, so overlap is counted once.
    """
    total = 0
    end = None
    for a, b in sorted(intervals):
        if end is None or a > end:
            total += b - a
            end = b
        elif b > end:
            total += b - end
            end = b
    return total


def drop_ambiguous(blocks: list[tuple[int, int, int]]) -> tuple[list[tuple[int, int, int]], int]:
    """Remove every region this aligner placed in more than one way.

    The inverted repeats mean parts of the genome legitimately align in two places. A
    crosswalk that kept both would silently pick one, which is precisely the error this
    repository exists to prevent, so those regions are dropped and the lift-over declines
    there instead. The dropped total is reported, not hidden.

    **This filter sees only what the aligner emits.** `align` runs minimap2 with
    `--secondary=no`, so a second placement is discarded upstream and never reaches
    here: round 7 measured roughly 92 per cent of positions inside the inverted repeats
    still receiving one confident answer. Repeat ambiguity is therefore handled by
    `liftover.lift` consulting the repeat map, not by this function, and the count below
    is only about regions a single alignment pass placed twice.
    """
    bad_q = ambiguous_intervals([(q, q + ln) for q, _, ln in blocks])
    bad_t = ambiguous_intervals([(t, t + ln) for _, t, ln in blocks])
    if not bad_q and not bad_t:
        return blocks, 0
    kept: list[tuple[int, int, int]] = []
    for b in blocks:
        kept.extend(subtract(b, bad_q, bad_t))
    dropped = (union_length([(q, q + ln) for q, _, ln in blocks])
               - union_length([(q, q + ln) for q, _, ln in kept]))
    return sorted(kept), dropped


def drop_non_monotonic(blocks: list[tuple[int, int, int]]) -> tuple[list[tuple[int, int, int]], int]:
    """Keep the largest chain of blocks that is monotonic in both sequences.

    Within one strand the mapping must not jump backwards. A block whose target order
    contradicts its query order has been placed on the other copy of a repeat: review
    found a position 900 bp into one genome resolving to 151,879, the far copy of the
    terminal repeat, which is the error this layer exists to prevent.

    A greedy sweep is not enough, because replacing one block can create a fresh
    backwards jump against the block before it. This is the standard chaining problem,
    solved exactly here by a longest-increasing-subsequence weighted by block length.
    The genomes have a few hundred blocks each, so the quadratic solution is fine and
    is worth the certainty.
    """
    if not blocks:
        return blocks, 0
    ordered = sorted(blocks)
    n = len(ordered)
    best = [b[2] for b in ordered]          # best total length of a chain ending here
    prev = [-1] * n
    for i in range(n):
        qs_i, ts_i, len_i = ordered[i]
        for j in range(i):
            qs_j, ts_j, len_j = ordered[j]
            if qs_j + len_j <= qs_i and ts_j + len_j <= ts_i and best[j] + len_i > best[i]:
                best[i] = best[j] + len_i
                prev[i] = j
    end = max(range(n), key=lambda i: best[i])
    chain: list[tuple[int, int, int]] = []
    while end != -1:
        chain.append(ordered[end])
        end = prev[end]
    chain.reverse()
    dropped = sum(b[2] for b in ordered) - sum(b[2] for b in chain)
    return chain, dropped


def align(query_fa: pathlib.Path, ref_fa: pathlib.Path) -> list[list[str]]:
    proc = subprocess.run(
        ["minimap2", "-cx", PRESET, "--cs", "-t", "4", "--secondary=no",
         str(ref_fa), str(query_fa)],
        capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {proc.stderr[-300:]}")
    return [line.split("\t") for line in proc.stdout.splitlines() if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=0, help="only process N genomes")
    args = ap.parse_args()

    if not shutil.which("minimap2"):
        print("minimap2 not found; it is a declared dependency of this milestone",
              file=sys.stderr)
        return 1

    ref_seq = genbank_sequence(RAW / "ref_NC_001806.2.gb")
    genomes = read_fasta(RAW / "genomes.fasta")
    print(f"  reference {REFERENCE}: {len(ref_seq):,} bp")
    print(f"  {len(genomes)} genomes to map")

    CROSSWALK.mkdir(parents=True, exist_ok=True)
    tmp = pathlib.Path(tempfile.mkdtemp())
    ref_fa = tmp / "ref.fa"
    ref_fa.write_text(f">{REFERENCE}\n{ref_seq}\n")

    coverage: list[dict[str, str]] = []
    items = sorted(genomes.items())
    if args.limit:
        items = items[:args.limit]

    for i, (acc, seq) in enumerate(items, 1):
        q_fa = tmp / "q.fa"
        q_fa.write_text(f">{acc}\n{seq}\n")
        deposited_length = len(seq)
        rows = align(q_fa, ref_fa)
        orientation = "forward"
        if rows and all(f[4] == "-" for f in rows):
            # Deposited on the opposite strand. Reverse complement and realign, rather
            # than reporting zero coverage as though the assembly were unusable.
            seq = seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
            q_fa.write_text(f">{acc}\n{seq}\n")
            rows = align(q_fa, ref_fa)
            orientation = "reverse_complemented"
        blocks: list[tuple[int, int, int]] = []
        for f in rows:
            if f[4] != "+":          # reverse-strand hits are repeat artefacts here
                continue
            cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
            if not cs:
                continue
            blocks.extend(blocks_from_cs(int(f[2]), int(f[7]), cs))
        blocks.sort()
        blocks, dropped = drop_ambiguous(blocks)
        blocks, non_mono = drop_non_monotonic(blocks)
        mapped = sum(b[2] for b in blocks)
        out = CROSSWALK / f"{acc}.blocks.tsv"
        with out.open("w") as fh:
            # The frame these query offsets are in, written beside them. Round 7: the
            # three reverse-complemented genomes stored offsets into the reverse
            # complement while the only record of that sat in a different file, which
            # nothing read, so a lift was wrong by tens of kilobases. A coordinate's
            # frame now travels with the coordinate.
            fh.write(f"# orientation={orientation} "
                     f"deposited_length={deposited_length}\n")
            fh.write("query_start\ttarget_start\tlength\n")
            for qs, ts, ln in blocks:
                fh.write(f"{qs}\t{ts}\t{ln}\n")
        coverage.append({
            "accession": acc,
            "length_bp": str(deposited_length),
            "blocks": str(len(blocks)),
            "mapped_bp": str(mapped),
            "mapped_percent": f"{100 * mapped / deposited_length:.2f}",
            "alignments": str(len(rows)),
            "ambiguous_bp_dropped": str(dropped),
            "non_monotonic_bp_dropped": str(non_mono),
            "orientation": orientation,
        })
        if i % 25 == 0 or i == len(items):
            print(f"    {i}/{len(items)}")

    cov = DERIVED / "crosswalk_coverage.tsv"
    cols = ["accession", "length_bp", "blocks", "mapped_bp", "mapped_percent",
            "alignments", "ambiguous_bp_dropped", "non_monotonic_bp_dropped",
            "orientation"]
    with cov.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in coverage:
            fh.write("\t".join(r[c] for c in cols) + "\n")

    pcts = [float(r["mapped_percent"]) for r in coverage]
    pcts_mappable = [float(r["mapped_percent"]) for r in coverage
                     if not (float(r["mapped_percent"]) == 0
                             and int(r["alignments"]) > 0)]
    mappable = [r for r in coverage
                if not (float(r["mapped_percent"]) == 0 and int(r["alignments"]) > 0)]
    low = sorted(mappable, key=lambda r: float(r["mapped_percent"]))[:10]
    report = DERIVED / "crosswalk_report.md"
    with report.open("w") as fh:
        fh.write("# Coordinate crosswalk: coverage\n\n")
        fh.write(f"Every genome aligned to `{REFERENCE}` with minimap2 (`{PRESET}`, "
                 f"`--cs`, `--secondary=no`), reduced to ungapped blocks by "
                 f"`src/build_crosswalk.py`.\n\n")
        fh.write(f"- Genomes mapped: **{len(coverage)}**\n")
        fh.write(f"- Median mappable fraction: "
                 f"**{statistics.median(pcts_mappable):.2f}%**\n")
        fh.write(f"- Range: {min(pcts_mappable):.2f}% to {max(pcts_mappable):.2f}%\n")
        fh.write(f"- Genomes below 95% mappable: "
                 f"**{sum(1 for p in pcts if p < 95)}**\n\n")
        amb = [int(r["ambiguous_bp_dropped"]) for r in coverage]
        mono = [int(r["non_monotonic_bp_dropped"]) for r in coverage]
        fh.write(f"- Genomes with regions this alignment placed in more than one way: "
                 f"**{sum(1 for a in amb if a)}**, median "
                 f"{statistics.median([a for a in amb if a]) if any(amb) else 0:,.0f} bp "
                 f"removed where present\n")
        fh.write(f"- Genomes with blocks dropped by chaining: "
                 f"**{sum(1 for a in mono if a)}**, "
                 f"{sum(mono):,} bp in total\n\n")
        fh.write("**What this does not do.** minimap2 runs with `--secondary=no`, so a "
                 "second placement of the same query region is discarded by the aligner "
                 "before the ambiguity filter sees it. The counts above are therefore "
                 "about regions a single alignment pass placed twice, and they are not a "
                 "measure of repeat ambiguity: round 7 measured roughly 92 per cent of "
                 "positions inside the inverted repeats still receiving one confident "
                 "block. Repeat ambiguity is handled by `liftover.lift` consulting the "
                 "repeat map, which declines where a position has more than one valid "
                 "copy, and not by this build.\n\n")
        fh.write("A position that falls outside every block has no reference coordinate, "
                 "and the lift-over returns nothing for it rather than the nearest "
                 "match.\n\n")
        inverted = [r for r in coverage
                    if float(r["mapped_percent"]) == 0 and int(r["alignments"]) > 0]
        if inverted:
            fh.write("## Genomes in a different isomer orientation\n\n")
            fh.write(f"**{len(inverted)} genomes map at 0 per cent while having "
                     f"alignments.** These are not poor assemblies. The HSV-1 genome "
                     "exists as four isomers, because the long and short segments invert "
                     "independently, and a genome deposited in a different arrangement "
                     "from the reference has its alignments on the reverse strand, which "
                     "this build discards.\n\n")
            fh.write("Handling them needs per-segment orientation, which belongs with the "
                     "isomer-aware mask in milestone 4. Until then the crosswalk declines "
                     "for them entirely, which is correct but incomplete.\n\n")
            fh.write("| Accession | Length | Alignments found | Blocks kept |\n")
            fh.write("|---|---|---|---|\n")
            for r in inverted:
                fh.write(f"| `{r['accession']}` | {int(r['length_bp']):,} | "
                         f"{r['alignments']} | {r['blocks']} |\n")
            fh.write("\n")

        fh.write("## Ten least mappable genomes\n\n")
        fh.write("Excluding any genome that has alignments but maps at 0 per cent, "
                 "listed above where such genomes exist.\n\n")
        fh.write("| Accession | Length | Mapped | Blocks |\n|---|---|---|---|\n")
        for r in low:
            fh.write(f"| `{r['accession']}` | {int(r['length_bp']):,} | "
                     f"{r['mapped_percent']}% | {r['blocks']} |\n")

    print(f"\n  median mappable {statistics.median(pcts):.2f}%, "
          f"{sum(1 for p in pcts if p < 95)} below 95%")
    print(f"  wrote {cov.name}, {report.name}, and {len(coverage)} block files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
