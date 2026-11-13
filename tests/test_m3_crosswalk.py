"""The M3 gate: the crosswalk round-trips, and declines rather than guessing."""
from __future__ import annotations

import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import Crosswalk, lift  # noqa: E402
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
        cls.ref = genbank_sequence(ROOT / "data" / "raw" / "ref_NC_001806.2.gb")
        cls.genomes = read_fasta(ROOT / "data" / "raw" / "genomes.fasta")

    def test_lifted_bases_match_far_more_often_than_chance(self):
        rng = random.Random(SEED + 7)
        same = total = 0
        for acc in Crosswalk.available(CROSSWALK)[:30]:
            seq = self.genomes.get(acc)
            if not seq:
                continue
            cw = Crosswalk.load(acc, CROSSWALK)
            if not cw.blocks:
                continue
            for _ in range(100):
                b = rng.choice(cw.blocks)
                pos = b.query_start + rng.randrange(b.length) + 1
                ref_pos = cw.to_reference(pos)
                if ref_pos is None:
                    continue
                a, b_ = seq[pos - 1].upper(), self.ref[ref_pos - 1].upper()
                if a not in "ACGT" or b_ not in "ACGT":
                    continue      # ambiguous base in a deposited assembly, not a lift error
                total += 1
                if a == b_:
                    same += 1
        self.assertGreater(total, 1000, "too few positions lifted to judge")
        rate = 100 * same / total
        self.assertGreaterEqual(
            rate, 95.0,
            f"only {rate:.2f}% of {total} lifted positions have a matching base; "
            f"a frame offset would show here")

    def test_a_deliberate_offset_is_rejected(self):
        """Proof the test above has power: shift by one and it must fail the bar."""
        rng = random.Random(SEED + 8)
        same = total = 0
        for acc in Crosswalk.available(CROSSWALK)[:15]:
            seq = self.genomes.get(acc)
            if not seq:
                continue
            cw = Crosswalk.load(acc, CROSSWALK)
            if not cw.blocks:
                continue
            for _ in range(100):
                b = rng.choice(cw.blocks)
                pos = b.query_start + rng.randrange(b.length) + 1
                ref_pos = cw.to_reference(pos)
                if ref_pos is None or ref_pos + 1 > len(self.ref):
                    continue
                a, b_ = seq[pos - 1].upper(), self.ref[ref_pos].upper()
                if a not in "ACGT" or b_ not in "ACGT":
                    continue
                total += 1
                if a == b_:                                   # deliberate +1 offset
                    same += 1
        rate = 100 * same / total
        self.assertLess(rate, 60.0,
                        f"a one-base offset still agreed {rate:.2f}% of the time, so "
                        f"the sequence check cannot detect a frame error")

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

    def test_identity_lift_returns_the_position_only_where_it_is_mappable(self):
        """The identity case is still a claim, so it declines where it has no basis."""
        acc = Crosswalk.available(CROSSWALK)[0]
        cw = Crosswalk.load(acc, CROSSWALK)
        inside = cw.blocks[0].query_start + 1
        self.assertEqual(lift(inside, acc, acc, CROSSWALK), inside)
        beyond = max(b.query_start + b.length for b in cw.blocks) + 5000
        self.assertIsNone(lift(beyond, acc, acc, CROSSWALK),
                          "identity lift must not confirm an unmappable position")

    def test_reference_is_a_valid_endpoint(self):
        acc = Crosswalk.available(CROSSWALK)[0]
        cw = Crosswalk.load(acc, CROSSWALK)
        pos = cw.blocks[0].query_start + 1
        self.assertEqual(lift(pos, acc, REFERENCE, CROSSWALK), cw.to_reference(pos))


if __name__ == "__main__":
    unittest.main(verbosity=2)
