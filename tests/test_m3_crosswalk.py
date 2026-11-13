"""The M3 gate: the crosswalk round-trips, and declines rather than guessing."""
from __future__ import annotations

import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import Crosswalk, RepeatMap, lift  # noqa: E402
from hsv1coords.liftover import REFERENCE  # noqa: E402

sys.path.insert(0, str(ROOT / "src"))
from build_crosswalk import genbank_sequence, read_fasta  # noqa: E402

CROSSWALK = ROOT / "data" / "derived" / "crosswalk"
SEED = 20261113        # fixed, so a failure is reproducible


class CrosswalkExists(unittest.TestCase):
    def test_every_genome_has_blocks(self):
        available = Crosswalk.available(CROSSWALK)
        self.assertGreater(len(available), 150,
                           f"only {len(available)} crosswalks built")

    def test_blocks_are_well_formed(self):
        for acc in Crosswalk.available(CROSSWALK)[:20]:
            cw = Crosswalk.load(acc, CROSSWALK)
            with self.subTest(accession=acc):
                self.assertGreater(len(cw.blocks), 0)
                for b in cw.blocks:
                    self.assertGreater(b.length, 0)
                    self.assertGreaterEqual(b.query_start, 0)
                    self.assertGreaterEqual(b.target_start, 0)

    def test_blocks_do_not_overlap_in_either_sequence(self):
        """Overlapping blocks would let one position map to two places silently."""
        for acc in Crosswalk.available(CROSSWALK)[:20]:
            cw = Crosswalk.load(acc, CROSSWALK)
            with self.subTest(accession=acc):
                q = sorted((b.query_start, b.query_start + b.length) for b in cw.blocks)
                for (_, end), (start, _) in zip(q, q[1:]):
                    self.assertLessEqual(end, start, "query blocks overlap")
                t = sorted((b.target_start, b.target_start + b.length) for b in cw.blocks)
                for (_, end), (start, _) in zip(t, t[1:]):
                    self.assertLessEqual(end, start, "target blocks overlap")


class RoundTrip(unittest.TestCase):
    """Criterion: round-trip returns the original coordinate for at least 99 per cent
    of positions outside the unmappable regions, with failures enumerated."""

    def test_round_trip_through_the_reference(self):
        rng = random.Random(SEED)
        failures: list[tuple[str, int, int | None]] = []
        tested = 0
        for acc in Crosswalk.available(CROSSWALK)[:40]:
            cw = Crosswalk.load(acc, CROSSWALK)
            for _ in range(250):
                b = rng.choice(cw.blocks)
                pos = b.query_start + rng.randrange(b.length) + 1
                tested += 1
                ref = cw.to_reference(pos)
                if ref is None:
                    failures.append((acc, pos, None))
                    continue
                back = cw.from_reference(ref)
                if back != pos:
                    failures.append((acc, pos, back))
        rate = 100 * (tested - len(failures)) / tested
        self.assertGreaterEqual(
            rate, 99.0,
            f"round-trip succeeded for {rate:.2f}% of {tested} positions; "
            f"first failures: {failures[:10]}")

    def test_positions_in_blocks_always_lift(self):
        rng = random.Random(SEED + 1)
        for acc in Crosswalk.available(CROSSWALK)[:15]:
            cw = Crosswalk.load(acc, CROSSWALK)
            for _ in range(100):
                b = rng.choice(cw.blocks)
                pos = b.query_start + rng.randrange(b.length) + 1
                with self.subTest(accession=acc, position=pos):
                    self.assertIsNotNone(cw.to_reference(pos))


class DeclinesRatherThanGuessing(unittest.TestCase):
    """Criterion: returns no confident coordinate where it has none."""

    def test_position_past_the_last_block_returns_none(self):
        for acc in Crosswalk.available(CROSSWALK)[:15]:
            cw = Crosswalk.load(acc, CROSSWALK)
            last = max(b.query_start + b.length for b in cw.blocks)
            with self.subTest(accession=acc):
                self.assertIsNone(cw.to_reference(last + 1000))

    def test_position_before_the_first_block_returns_none(self):
        for acc in Crosswalk.available(CROSSWALK)[:15]:
            cw = Crosswalk.load(acc, CROSSWALK)
            first = min(b.query_start for b in cw.blocks)
            if first == 0:
                continue
            with self.subTest(accession=acc):
                self.assertIsNone(cw.to_reference(0))

    def test_gaps_between_blocks_return_none(self):
        """A position in an insertion relative to the reference has no coordinate."""
        checked = 0
        for acc in Crosswalk.available(CROSSWALK)[:25]:
            cw = Crosswalk.load(acc, CROSSWALK)
            blocks = cw.blocks
            for a, b in zip(blocks, blocks[1:]):
                gap_start = a.query_start + a.length
                if b.query_start > gap_start:
                    checked += 1
                    with self.subTest(accession=acc, position=gap_start + 1):
                        self.assertIsNone(cw.to_reference(gap_start + 1))
                    break
        self.assertGreater(checked, 0, "no insertion gaps found to test")


class Monotonicity(unittest.TestCase):
    """Within one strand the mapping must not jump backwards.

    A block whose target order contradicts its query order has been placed on the other
    copy of a repeat. Review found exactly this: a position 900 bp into one genome
    resolved to 151,879, the far copy of the terminal repeat, which is the error this
    whole layer exists to prevent.
    """

    def test_no_crosswalk_jumps_backwards(self):
        offenders = []
        for acc in Crosswalk.available(CROSSWALK):
            blocks = Crosswalk.load(acc, CROSSWALK).blocks
            for a, b in zip(blocks, blocks[1:]):
                if b.target_start < a.target_start + a.length:
                    offenders.append((acc, a.target_start, b.target_start))
                    break
        self.assertEqual(offenders, [],
                         f"{len(offenders)} crosswalks are non-monotonic: "
                         f"{offenders[:5]}")

    def test_the_known_cross_repeat_case_now_declines(self):
        if "OP297884.1" not in Crosswalk.available(CROSSWALK):
            self.skipTest("OP297884.1 not in the pinned genome set")
        self.assertIsNone(
            lift(900, "OP297884.1", REFERENCE, CROSSWALK),
            "position 900 must not resolve to the far copy of the terminal repeat")


class SequenceAgreement(unittest.TestCase):
    """The round-trip test alone is a mathematical identity and cannot fail.

    Blocks are non-overlapping in both sequences by construction, so lifting out and
    back is guaranteed to return the original coordinate whatever the blocks say. That
    gate has no power on its own. This one does: if a lift is correct, the base at the
    source position and the base at the lifted position are usually the same, and the
    surrounding window is nearly identical. A frame offset of even one base destroys
    both, so this is what actually tests the coordinate mapping.
    """

    @classmethod
    def setUpClass(cls):
        import unittest as _ut
        for _need in ("ref_NC_001806.2.gb", "genomes.fasta"):
            if not (ROOT / "data" / "raw" / _need).exists():
                raise _ut.SkipTest(
                    f"{_need} is fetched, not committed. Run `make inputs` "
                    f"and `make genomes` first.")
        cls.ref = genbank_sequence(ROOT / "data" / "raw" / "ref_NC_001806.2.gb")
        cls.genomes = read_fasta(ROOT / "data" / "raw" / "genomes.fasta")

    COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

    def _base_agreement(self, acc: str, step: int,
                        min_block: int = 200) -> tuple[int, int]:
        """Matching bases over every `step`-th mapped position of one genome.

        A reverse-complemented genome's bases are complemented as well as reordered,
        so the deposited base is compared after complementing. Getting this wrong makes
        a correct crosswalk look like noise, which it did on the first attempt at this
        check.
        """
        seq = self.genomes.get(acc)
        if not seq:
            return 0, 0
        seq = seq.upper()
        cw = Crosswalk.load(acc, CROSSWALK)
        rc = cw.orientation == "reverse_complemented"
        same = total = 0
        for b in cw.blocks:
            # Blocks shorter than `min_block` are skipped. In the two most divergent
            # genomes the aligner fragments into hundreds of 3 to 17 bp blocks, where a
            # single real substitution scores 67 per cent and drags a per-genome rate
            # below any sensible threshold: 2,187 bp of 143,459 in 159 such blocks for
            # ON783214.1. A frame error destroys agreement in long blocks too, so
            # restricting to them keeps the sensitivity and drops the noise.
            if b.length < min_block:
                continue
            for off in range(0, b.length, step):
                dep = cw._to_deposited(b.query_start + off + 1)
                ref_pos = cw.to_reference(dep)
                if ref_pos is None:
                    continue
                a = seq[dep - 1]
                if rc:
                    a = a.translate(self.COMPLEMENT)
                b_ = self.ref[ref_pos - 1].upper()
                if a not in "ACGT" or b_ not in "ACGT":
                    continue   # ambiguous base in a deposited assembly, not a lift error
                total += 1
                same += (a == b_)
        return same, total

    def test_lifted_bases_match_in_every_genome_not_just_the_first_thirty(self):
        """Per genome, over all 183, and that is how round 7 found the frame bug.

        This gate used to score the alphabetically-first 30 genomes and pool the
        result. Both halves hid the defect: the three reverse-complemented genomes sort
        outside the first 30, and pooling 183 genomes means three at 28 per cent cannot
        pull the total below a 95 per cent threshold. A pooled average cannot locate a
        broken genome, so this enumerates them.
        """
        poor = []
        checked = grand_total = 0
        for acc in Crosswalk.available(CROSSWALK):
            same, total = self._base_agreement(acc, step=37)
            if total < 300:
                continue
            checked += 1
            grand_total += total
            rate = 100 * same / total
            if rate < 97.0:
                poor.append((acc, round(rate, 2), total))
        self.assertGreaterEqual(checked, 180,
                                f"only {checked} genomes were scored; the corpus or "
                                f"the sampling changed")
        self.assertEqual(poor, [],
                         f"{len(poor)} genomes disagree with the reference at more "
                         f"than 3 per cent of mapped positions: {poor}")
        print(f"\n    base agreement: {checked} genomes, all at or above 97%, "
              f"{grand_total:,} positions scored")

    def test_mutated_crosswalks_are_rejected(self):
        """Real mutation testing, replacing a control that had no power.

        The control this replaces shifted the LOOKUP by one base and checked that
        agreement collapsed. That measures only that DNA is not self-similar under a
        shift, which is a property of the sequence: it passes on a completely random
        crosswalk. FINDINGS.md claimed for a time that it had already been replaced
        here, while it sat unchanged. Closing a defect in the retraction document
        instead of in the code is its own failure mode, recorded in PITFALLS.md.

        These mutants perturb the crosswalk itself.
        """
        from hsv1coords.liftover import Block, Crosswalk as CW
        rng = random.Random(SEED + 11)
        acc = next(a for a in Crosswalk.available(CROSSWALK)
                   if len(Crosswalk.load(a, CROSSWALK).blocks) > 10
                   and self.genomes.get(a))
        real = Crosswalk.load(acc, CROSSWALK)
        seq = self.genomes[acc]

        def worst_block_identity(cw) -> float:
            """The weakest block, not the average.

            A genome-wide average hides a single wrong block: shifting one block of
            one genome moves the overall figure from about 99 to 97.6 per cent, which
            any sane threshold tolerates. The error is local, so the measurement has
            to be local too. Blocks too short to measure are skipped rather than
            allowed to dominate.
            """
            worst = 100.0
            for b in cw.blocks:
                if b.length < 200:
                    continue
                same = total = 0
                for off in range(0, b.length, max(1, b.length // 40)):
                    pos = b.query_start + off + 1
                    ref_pos = cw.to_reference(pos)
                    if ref_pos is None or ref_pos > len(self.ref):
                        continue
                    a, c = seq[pos - 1].upper(), self.ref[ref_pos - 1].upper()
                    if a not in "ACGT" or c not in "ACGT":
                        continue
                    total += 1
                    same += a == c
                if total >= 20:
                    worst = min(worst, 100 * same / total)
            return worst

        blocks = real.blocks
        mid = len(blocks) // 2
        mutants = {
            "one block shifted by 1": [
                Block(b.query_start, b.target_start + (1 if i == mid else 0), b.length)
                for i, b in enumerate(blocks)],
            "one block shifted by 500": [
                Block(b.query_start, b.target_start + (500 if i == mid else 0), b.length)
                for i, b in enumerate(blocks)],
            "every block shifted by 1": [
                Block(b.query_start, b.target_start + 1, b.length) for b in blocks],
            "targets randomised": [
                Block(b.query_start, rng.randrange(0, 140000), b.length)
                for b in blocks],
        }
        baseline = worst_block_identity(real)
        self.assertGreaterEqual(
            baseline, 90.0,
            f"the real crosswalk's weakest block scores {baseline:.2f}%")
        survivors = []
        for name, mb in mutants.items():
            score = worst_block_identity(CW(acc, mb))
            if score >= 90.0:
                survivors.append((name, round(score, 2)))
        self.assertEqual(survivors, [],
                         f"these mutated crosswalks passed the identity gate: "
                         f"{survivors}")

    def test_flanking_windows_agree(self):
        rng = random.Random(SEED + 9)
        poor = []
        for acc in Crosswalk.available(CROSSWALK)[:20]:
            seq = self.genomes.get(acc)
            if not seq:
                continue
            cw = Crosswalk.load(acc, CROSSWALK)
            if not cw.blocks:
                continue
            for _ in range(40):
                b = rng.choice(cw.blocks)
                if b.length < 120:
                    continue
                off = rng.randrange(50, b.length - 50)
                pos, ref_pos = b.query_start + off, b.target_start + off
                a = seq[pos - 50:pos + 50].upper()
                c = self.ref[ref_pos - 50:ref_pos + 50].upper()
                if len(a) != len(c) or not a:
                    continue
                # Compare only where both assemblies actually called a base. Several
                # deposited genomes carry runs of N inside the short repeats, and
                # scoring those as disagreement measures assembly completeness rather
                # than whether the coordinate is right.
                pairs = [(x, y) for x, y in zip(a, c) if x in "ACGT" and y in "ACGT"]
                if len(pairs) < 50:
                    continue
                ident = sum(1 for x, y in pairs if x == y) / len(pairs)
                if ident < 0.9:
                    poor.append((acc, pos, round(ident, 3), f"{len(pairs)}/100 called"))
        self.assertEqual(poor, [],
                         f"{len(poor)} windows below 90% identity: {poor[:5]}")


    def test_ambiguous_bases_are_reported_not_hidden(self):
        """How much of each genome is N is a property worth knowing, so state it."""
        worst = []
        for acc in Crosswalk.available(CROSSWALK):
            seq = self.genomes.get(acc)
            if not seq:
                continue
            n = sum(1 for c in seq.upper() if c not in "ACGT")
            if n:
                worst.append((acc, n, round(100 * n / len(seq), 3)))
        worst.sort(key=lambda r: -r[1])
        print(f"\n    {len(worst)} of {len(self.genomes)} genomes contain ambiguous "
              f"bases; worst: {worst[:3]}")
        self.assertLess(worst[0][2] if worst else 0, 50.0,
                        "a genome that is mostly N should not be in the crosswalk")


class ZeroBlockGenomes(unittest.TestCase):
    """A genome with no blocts must degrade cleanly, whether or not any exist today.

    Three genomes had none until the build learned to reverse complement deposits made
    on the opposite strand. The handling is still required: a future deposit could be
    unmappable for another reason, and the code must decline rather than raise.
    """

    def test_no_genome_is_left_unmapped_without_explanation(self):
        import csv as _csv
        with (ROOT / "data" / "derived" / "crosswalk_coverage.tsv").open() as fh:
            rows = list(_csv.DictReader(fh, delimiter="\t"))
        unmapped = [r["accession"] for r in rows if float(r["mapped_percent"]) == 0]
        self.assertEqual(unmapped, [],
                         f"{len(unmapped)} genomes map at zero per cent: {unmapped}")

    def test_an_empty_crosswalk_declines_rather_than_raising(self):
        from hsv1coords.liftover import Crosswalk as _CW
        cw = _CW("SYNTHETIC.0", [])
        self.assertIsNone(cw.to_reference(1000))
        self.assertIsNone(cw.from_reference(1000))
        self.assertEqual(cw.mapped_bases, 0)

    def test_real_zero_block_genomes_if_any_are_handled(self):
        empty = [a for a in Crosswalk.available(CROSSWALK)
                 if not Crosswalk.load(a, CROSSWALK).blocks]
        for acc in empty:
            cw = Crosswalk.load(acc, CROSSWALK)
            with self.subTest(accession=acc):
                self.assertIsNone(cw.to_reference(1000))
                self.assertIsNone(cw.from_reference(1000))
                self.assertEqual(cw.mapped_bases, 0)
                self.assertIsNone(lift(1000, acc, REFERENCE, CROSSWALK))


class GenomeToGenome(unittest.TestCase):
    def test_lift_between_two_genomes_is_consistent(self):
        accs = Crosswalk.available(CROSSWALK)
        a, b = accs[0], accs[1]
        cw = Crosswalk.load(a, CROSSWALK)
        rng = random.Random(SEED + 2)
        agreed = total = 0
        for _ in range(300):
            blk = rng.choice(cw.blocks)
            pos = blk.query_start + rng.randrange(blk.length) + 1
            there = lift(pos, a, b, CROSSWALK)
            if there is None:
                continue
            total += 1
            if lift(there, b, a, CROSSWALK) == pos:
                agreed += 1
        self.assertGreater(total, 50, "too few positions lifted to judge")
        rate = 100 * agreed / total
        self.assertGreaterEqual(rate, 99.0,
                                f"genome-to-genome round-trip agreed {rate:.2f}%")

    @staticmethod
    def _unique_position(cw):
        """A deposited position whose reference counterpart names one locus.

        Round 7 made `lift` decline inside the repeats, so a test that wants a number
        back has to ask for a position outside them. Picking `blocks[0]` no longer does:
        on most genomes the first block is inside TRL.
        """
        rm = RepeatMap.load()
        for b in cw.blocks:
            for off in (0, b.length // 2, b.length - 1):
                dep = cw._to_deposited(b.query_start + off + 1)
                ref = cw.to_reference(dep)
                if ref is not None and rm.classify(ref).well_defined:
                    return dep
        raise AssertionError(f"no uniquely defined position found in {cw.accession}")

    def test_identity_lift_returns_the_position_only_where_it_is_mappable(self):
        """The identity case is still a claim, so it declines where it has no basis."""
        acc = Crosswalk.available(CROSSWALK)[0]
        cw = Crosswalk.load(acc, CROSSWALK)
        inside = self._unique_position(cw)
        self.assertEqual(lift(inside, acc, acc, CROSSWALK), inside)
        beyond = (cw.deposited_length or
                  max(b.query_start + b.length for b in cw.blocks)) + 5000
        self.assertIsNone(lift(beyond, acc, acc, CROSSWALK),
                          "identity lift must not confirm an unmappable position")

    def test_reference_is_a_valid_endpoint(self):
        acc = Crosswalk.available(CROSSWALK)[0]
        cw = Crosswalk.load(acc, CROSSWALK)
        pos = self._unique_position(cw)
        self.assertEqual(lift(pos, acc, REFERENCE, CROSSWALK), cw.to_reference(pos))

    def test_lift_declines_inside_a_repeat_while_the_alignment_does_not(self):
        """The difference between the alignment's answer and an honest one.

        `--secondary=no` means minimap2 reports one placement inside the inverted
        repeats and the ambiguity filter never sees the other, so `Crosswalk` returns a
        confident number there. Round 7 measured roughly 92 per cent of repeat
        positions receiving one. `lift` consults the repeat map and declines, and this
        is the test that fails if it stops.
        """
        rm = RepeatMap.load()
        declined = alignment_answered = 0
        for acc in Crosswalk.available(CROSSWALK)[:40]:
            cw = Crosswalk.load(acc, CROSSWALK)
            for b in cw.blocks:
                dep = cw._to_deposited(b.query_start + b.length // 2 + 1)
                ref = cw.to_reference(dep)
                if ref is None or rm.classify(ref).well_defined:
                    continue
                alignment_answered += 1
                if lift(dep, acc, REFERENCE, CROSSWALK) is None:
                    declined += 1
        self.assertGreater(alignment_answered, 200,
                           "too few in-repeat positions to judge; the corpus changed")
        self.assertEqual(declined, alignment_answered,
                         f"lift answered at {alignment_answered - declined} of "
                         f"{alignment_answered} positions the repeat map calls ambiguous")


if __name__ == "__main__":
    unittest.main(verbosity=2)
