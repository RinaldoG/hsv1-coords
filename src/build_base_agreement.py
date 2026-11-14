"""Measure every genome's base agreement with the reference and pin it to a file.

Why this exists. The agreement gate used to assert a floor: no genome below 95 per cent.
The corpus sits between 95.8 and 99.5, so every genome carried a three to four point
corruption budget, and the 1-in-37 stride on long blocks turns one undetected sample into
37 bp of damage. Review round 9 spent that budget deliberately: 1,764,625 bp shifted by
one base across 136 of the 183 genomes, with all 180 tests green, and 1,543,561 bp shifted
by up to 500 bp likewise. A floor set below what the data achieves is not a gate on the
data, it is a gate on how badly the data would have to be broken before anyone noticed.

So the measurement is recorded per genome, exactly, and the gate asserts equality rather
than a threshold. A single shifted block changes the `ok` count for its genome and the
gate names that genome. Regenerating this file to make a failure go away is a visible diff
against a committed artefact, which is the point: the laundering has to be done on the
record rather than inside a tolerance.

`deposited_length` is written here as well, and checked against the sequence. Nothing
validated it before: setting one genome's header to a thousandth of its true length left
the suite green while that genome lost 148,809 of its 149,490 liftable positions.

Run: `make agreement` (needs `data/raw/genomes.fasta` and the reference record).
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from build_crosswalk import genbank_sequence, read_fasta  # noqa: E402
from hsv1coords.liftover import Crosswalk  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
CROSSWALK = DERIVED / "crosswalk"
OUT = DERIVED / "base_agreement.tsv"

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
SHORT_BLOCK = 200
STRIDE = 37

COLUMNS = ["accession", "orientation", "deposited_length", "sequence_length",
           "blocks", "mapped_bp", "sampled_blocks", "unsampled_bp", "thin_blocks",
           "clustered_blocks", "worst_sample_gap", "sampled_positions", "scored",
           "matching", "short_scored", "short_matching"]

MIN_SAMPLES = 5
"""The fewest samples a block may contribute unless it is scored at every base.

Round 10 found this inert, which makes it the third version of this check that could not
fail. With `SHORT_BLOCK = 200` and `STRIDE = 37`, a block under 200 bp is scored at every
base, so the "unless" clause is always satisfied, and a block of 200 bp or more yields at
least ceil(200/37) = 6 samples, so the floor of 5 is always met. No block of any length
could be counted thin, and `MIN_SAMPLES` could be set anywhere from 2 to 6 with the suite
green.

It is kept, because it is the right condition if those parameters change, and the gate
now also measures the thing the count could not see: the largest gap between consecutive
sampled positions inside a block. That is what distinguishes real coverage from a sample
that has been quietly clustered, which was the reviewer's successful attack: replacing
the stride with six samples packed into the first six bases of every block left every
printed coverage figure identical.
"""

MAX_SAMPLE_GAP = STRIDE
"""The largest allowed distance between consecutive sampled positions inside a block.

A stride of `STRIDE` gives gaps of exactly `STRIDE`; exhaustive scoring gives gaps of 1.
Anything larger means part of the block is further from a sampled base than the sampling
claims, which is how a clustered sample looks.
"""


def sampled_offsets(length: int) -> range:
    """The positions the gate looks at inside a block of this length.

    One definition, used by the builder and by the gate, so the gate's reported coverage
    cannot drift from what it actually samples. Blocks under `SHORT_BLOCK` are scored at
    every base: a step-sampled 3 bp block gets one sample, and a single real substitution
    in it reads as 67 per cent agreement.
    """
    return range(0, length, 1 if length < SHORT_BLOCK else STRIDE)


def measure(accession: str, sequence: str, reference: str,
            directory: pathlib.Path | None = None) -> dict:
    """Agreement for one genome, as counts rather than a percentage.

    Counts, not a rate, because a rate hides its denominator and the denominator is what
    round 8 and round 9 both found wrong.
    """
    cw = Crosswalk.load(accession, directory or CROSSWALK)
    rc = cw.orientation == "reverse_complemented"
    seq = sequence.upper()
    row = {"accession": accession, "orientation": cw.orientation,
           "deposited_length": cw.deposited_length or len(seq),
           "sequence_length": len(seq), "blocks": len(cw.blocks), "mapped_bp": 0,
           "sampled_blocks": 0, "unsampled_bp": 0, "thin_blocks": 0,
           "clustered_blocks": 0, "worst_sample_gap": 0,
           "sampled_positions": 0, "scored": 0, "matching": 0,
           "short_scored": 0, "short_matching": 0}
    for b in cw.blocks:
        row["mapped_bp"] += b.length
        offsets = sampled_offsets(b.length)
        if len(offsets) == 0:
            row["unsampled_bp"] += b.length
            continue
        row["sampled_blocks"] += 1
        row["sampled_positions"] += len(offsets)
        if len(offsets) < min(MIN_SAMPLES, b.length) and len(offsets) != b.length:
            row["thin_blocks"] += 1
        # The spread, not just the count. A block sampled six times in its first six
        # bases has the same count as one sampled every 37th base and tells you nothing
        # about the rest of itself.
        positions = list(offsets)
        gaps = [second - first
                for first, second in zip(positions, positions[1:])]
        worst = max(gaps) if gaps else b.length
        tail = b.length - positions[-1]
        row["worst_sample_gap"] = max(row["worst_sample_gap"], worst, tail)
        if worst > MAX_SAMPLE_GAP or tail > MAX_SAMPLE_GAP:
            row["clustered_blocks"] += 1
        short = b.length < SHORT_BLOCK
        for off in offsets:
            dep = cw._to_deposited(b.query_start + off + 1)
            ref_pos = cw.to_reference(dep)
            if ref_pos is None or dep < 1 or dep > len(seq):
                continue
            a = seq[dep - 1]
            if rc:
                a = a.translate(COMPLEMENT)
            c = reference[ref_pos - 1].upper()
            if a not in "ACGT" or c not in "ACGT":
                continue    # an ambiguous base in a deposited assembly, not a lift error
            row["scored"] += 1
            row["matching"] += a == c
            if short:
                row["short_scored"] += 1
                row["short_matching"] += a == c
    return row


def main() -> None:
    reference = genbank_sequence(RAW / "ref_NC_001806.2.gb")
    genomes = read_fasta(RAW / "genomes.fasta")
    rows = []
    mismatched_lengths = []
    for acc in Crosswalk.available(CROSSWALK):
        seq = genomes.get(acc)
        if not seq:
            continue
        row = measure(acc, seq, reference)
        if row["deposited_length"] != row["sequence_length"]:
            mismatched_lengths.append(
                (acc, row["deposited_length"], row["sequence_length"]))
        rows.append(row)
    if mismatched_lengths:
        raise SystemExit(
            f"{len(mismatched_lengths)} block files state a deposited_length that is "
            f"not the length of the deposited sequence: {mismatched_lengths}. A "
            f"reverse-complemented genome's coordinates are computed from this number, "
            f"so a wrong one moves every lift in that genome. Rebuild with "
            f"`make crosswalk`.")
    OUT.write_text(
        "# every genome's measured agreement with NC_001806.2, pinned exactly.\n"
        f"# sampling: blocks under {SHORT_BLOCK} bp at every base, longer blocks every "
        f"{STRIDE}th. Ambiguous bases in either sequence are not scored.\n"
        "# the gate in tests/test_m3_crosswalk.py asserts equality against this file,\n"
        "# not a threshold, because a threshold below what the corpus achieves is a\n"
        "# corruption budget. Round 9 spent it: 1,764,625 bp shifted, suite green.\n"
        + "\t".join(COLUMNS) + "\n"
        + "".join("\t".join(str(r[c]) for c in COLUMNS) + "\n" for r in rows))
    scored = sum(r["scored"] for r in rows)
    matching = sum(r["matching"] for r in rows)
    worst = min(rows, key=lambda r: r["matching"] / max(1, r["scored"]))
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} genomes, {scored:,} positions "
          f"scored, {100 * matching / scored:.2f} per cent pooled")
    print(f"  worst genome: {worst['accession']} at "
          f"{100 * worst['matching'] / worst['scored']:.2f} per cent over "
          f"{worst['scored']:,} positions")
    print(f"  unsampled: {sum(r['unsampled_bp'] for r in rows):,} bp of "
          f"{sum(r['mapped_bp'] for r in rows):,} mapped")


if __name__ == "__main__":
    main()
