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

Run: `make verify-crosswalk`. Exits non-zero and names every block that differs.
"""
from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from build_crosswalk import derive_blocks, genbank_sequence, read_fasta  # noqa: E402
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
    for need in ("ref_NC_001806.2.gb", "genomes.fasta"):
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
    if problems:
        print(f"\n{len(problems)} problems:", file=sys.stderr)
        for line in problems:
            print(f"  {line}", file=sys.stderr)
        return 1
    print("every committed block is the block an alignment gives")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
