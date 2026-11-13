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

    def test_identity_lift_is_the_identity(self):
        accs = Crosswalk.available(CROSSWALK)
        self.assertEqual(lift(1234, accs[0], accs[0], CROSSWALK), 1234)

    def test_reference_is_a_valid_endpoint(self):
        acc = Crosswalk.available(CROSSWALK)[0]
        cw = Crosswalk.load(acc, CROSSWALK)
        pos = cw.blocks[0].query_start + 1
        self.assertEqual(lift(pos, acc, REFERENCE, CROSSWALK), cw.to_reference(pos))


if __name__ == "__main__":
    unittest.main(verbosity=2)
