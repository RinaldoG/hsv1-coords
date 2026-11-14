"""Check every committed block file against a fresh alignment of its genome.

Why this exists. The crosswalk is the artefact the library reads, and until round 10
nothing checked it against an alignment. The round 9 baseline pins each genome's
base-agreement counts, which catches a block whose bases stop matching the reference. It
cannot catch a block moved to a different position where the bases still match, and
inside a tandem array the reference is periodic, so such positions exist: a reviewer
moved 34 blocks in 17 genomes that way, putting 3,868 positions at a wrong reference
coordinate by 1 to 48 bp, with the whole suite green.

The suite re-aligns a rotating sample of twelve genomes plus the four special cases,
because 0.08 s per genome is affordable per run but 183 of them is 14 s. This script does
the whole corpus, and CI runs it against the committed files **before** `make crosswalk`
overwrites them, which is the only moment at which a corrupted committed artefact is
still visible.

It also re-derives `version_crosswalk.tsv` the same way. Round 12 found that artefact
checked by no alignment at all, and its 18/19 junction sits inside a tandem array, so a
boundary could move 800 bp with every published figure and every gate intact.

Run: `make verify-crosswalk`. Exits non-zero and names every block that differs.
"""
from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from build_crosswalk import derive_blocks, genbank_sequence, read_fasta  # noqa: E402
from build_version_crosswalk import derive_blocks as derive_version_blocks  # noqa: E402
from hsv1coords.liftover import Crosswalk  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CROSSWALK = ROOT / "data" / "derived" / "crosswalk"


def main() -> int:
    if not shutil.which("minimap2"):
        print("minimap2 not found; it is a declared dependency of this milestone and "
              "without it the blocks can only be checked against themselves",
              file=sys.stderr)
        return 1
    for need in ("ref_NC_001806.2.gb", "ref_NC_001806.1.gb", "genomes.fasta"):
        if not (RAW / need).exists():
            print(f"missing {RAW / need}; run `make inputs` and `make genomes`",
                  file=sys.stderr)
            return 1
    reference = genbank_sequence(RAW / "ref_NC_001806.2.gb")
    genomes = read_fasta(RAW / "genomes.fasta")
    tmp = pathlib.Path(tempfile.mkdtemp())
    ref_fa = tmp / "ref.fa"
    ref_fa.write_text(">ref\n" + reference + "\n")

    problems: list[str] = []
    checked = 0
    for acc in Crosswalk.available(CROSSWALK):
        sequence = genomes.get(acc)
        if not sequence:
            problems.append(f"{acc}: a block file with no deposited sequence")
            continue
        blocks, orientation, reverse_bp, _, _ = derive_blocks(
            acc, sequence, tmp / "query.fa", ref_fa)
        cw = Crosswalk.load(acc, CROSSWALK)
        stored = sorted((b.query_start, b.target_start, b.length) for b in cw.blocks)
        derived = sorted(blocks)
        checked += 1
        if derived != stored:
            only_alignment = sorted(set(derived) - set(stored))
            only_file = sorted(set(stored) - set(derived))
            problems.append(
                f"{acc}: {len(derived)} blocks from the alignment against "
                f"{len(stored)} in the file; "
                f"{len(only_alignment)} only in the alignment "
                f"(first: {only_alignment[:2]}), "
                f"{len(only_file)} only in the file (first: {only_file[:2]})")
        if orientation != cw.orientation:
            problems.append(f"{acc}: file says {cw.orientation}, alignment says "
                            f"{orientation}")
        if reverse_bp != cw.reverse_strand_bp:
            problems.append(f"{acc}: file says {cw.reverse_strand_bp} bp dropped on the "
                            f"other strand, alignment says {reverse_bp}")
        if cw.deposited_length != len(sequence):
            problems.append(f"{acc}: file says deposited_length "
                            f"{cw.deposited_length}, the sequence is {len(sequence)}")

    print(f"re-aligned {checked} genomes against {len(Crosswalk.available(CROSSWALK))} "
          f"committed block files")
    # Round 11: with no floor, an empty crosswalk directory exited 0 and reported
    # success over nothing. Round 12: the floor was typed as 150 against a corpus of
    # 183, so deleting 30 block files printed "re-aligned 153 genomes against 153
    # committed block files" and exited 0. The count is now compared against the corpus
    # itself, so there is no headroom to delete into and no number to keep in step.
    missing = sorted(set(genomes) - set(Crosswalk.available(CROSSWALK)))
    if checked != len(genomes):
        print(f"{checked} genomes were re-aligned; the FASTA holds {len(genomes)}. "
              f"Every deposited sequence must have a committed block file that an "
              f"alignment reproduces. Without a block file for: "
              f"{missing[:5]}{' and more' if len(missing) > 5 else ''}", file=sys.stderr)
        return 1
    if problems:
        print(f"\n{len(problems)} problems:", file=sys.stderr)
        for line in problems:
            print(f"  {line}", file=sys.stderr)
        return 1
    if version_problems := verify_version_crosswalk(reference, tmp):
        print(f"\n{len(version_problems)} problems in the version crosswalk:",
              file=sys.stderr)
        for line in version_problems:
            print(f"  {line}", file=sys.stderr)
        return 1

    print("every committed block is the block an alignment gives")
    return 0


def verify_version_crosswalk(new_seq: str, tmp: pathlib.Path) -> list[str]:
    """The same check for `version_crosswalk.tsv`, which had none.

    Round 12: the strain crosswalk got a re-alignment check in round 10, and the version
    crosswalk, built by the same code path from the same two sequences, did not. Its
    18/19 junction sits inside a 22 bp tandem array, so a boundary could be moved 800 bp
    with coverage, block count, offsets and the 99.0 boundary-identity floor all intact
    and the suite green, putting 844 positions wrong. Identity cannot see this; the
    alignment can.
    """
    old_path = RAW / "ref_NC_001806.1.gb"
    if not old_path.exists():
        return [f"missing {old_path}; run `make inputs`"]
    committed = ROOT / "data" / "derived" / "version_crosswalk.tsv"
    if not committed.exists():
        return [f"missing {committed.relative_to(ROOT)}; run `make versions`"]
    old_seq = genbank_sequence(old_path)
    derived, _, _ = derive_version_blocks(old_seq, new_seq, tmp)
    stored = []
    for line in committed.read_text().splitlines()[1:]:
        if line.strip():
            qs, ts, ln = line.split("\t")[:3]
            stored.append((int(qs), int(ts), int(ln)))
    if sorted(derived) != sorted(stored):
        only_alignment = sorted(set(derived) - set(stored))
        only_file = sorted(set(stored) - set(derived))
        return [f"version_crosswalk.tsv: {len(derived)} blocks from the alignment "
                f"against {len(stored)} in the file; {len(only_alignment)} only in the "
                f"alignment (first: {only_alignment[:2]}), {len(only_file)} only in the "
                f"file (first: {only_file[:2]})"]
    print(f"the version crosswalk's {len(stored)} blocks are the blocks an alignment "
          f"of {len(old_seq):,} bp against {len(new_seq):,} bp gives")
    return []


if __name__ == "__main__":
    raise SystemExit(main())
