"""The M4 gate, written before the code, per ROADMAP.md.

The last criterion is deliberate: this milestone must demonstrably cover the three
defects that review found in M3, rather than merely claiming to.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import Crosswalk, RepeatMap, lift  # noqa: E402
from hsv1coords.liftover import REFERENCE  # noqa: E402
from hsv1coords.repeats import INVERTED_PAIRS  # noqa: E402

DERIVED = ROOT / "data" / "derived"
CROSSWALK = DERIVED / "crosswalk"


class ReproducesTheReference(unittest.TestCase):
    """Criterion: reproduces all four inverted repeat boundaries exactly."""

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)

    def test_the_four_inverted_repeats(self):
        expected = {"TRL": (1, 9213), "IRL": (117161, 126373),
                    "IRS": (125975, 132607), "TRS": (145590, 152222)}
        for label, (start, end) in expected.items():
            with self.subTest(repeat=label):
                r = next((x for x in self.m.regions
                          if x.name.startswith(label) and x.kind == "repeat_region"), None)
                self.assertIsNotNone(r, f"{label} missing from the map")
                self.assertEqual((r.start, r.end), (start, end))

    def test_the_pairs_are_equal_length(self):
        for a, b in INVERTED_PAIRS:
            ra = next(x for x in self.m.regions
                      if x.name.startswith(a) and x.kind == "repeat_region")
            rb = next(x for x in self.m.regions
                      if x.name.startswith(b) and x.kind == "repeat_region")
            with self.subTest(pair=f"{a}/{b}"):
                self.assertEqual(ra.length, rb.length)

    def test_coverage_is_a_fifth_of_the_genome(self):
        self.assertGreater(self.m.covered_bases, 30000)
        self.assertLess(self.m.covered_bases, 40000)


class DiploidPairs(unittest.TestCase):
    """Criterion: pairs are derived, and each pair's members have the same length."""

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)

    def test_pairs_exist(self):
        self.assertGreater(len(self.m.pairs), 5, "too few diploid pairs derived")

    def test_exact_pairs_have_equal_length_copies(self):
        for p in self.m.pairs:
            if not p.exact:
                continue
            with self.subTest(locus=p.locus):
                self.assertEqual(p.copy_1[1] - p.copy_1[0], p.copy_2[1] - p.copy_2[0])

    def test_partial_pairs_declare_themselves(self):
        """A pair whose copies differ in length must say so, not pass as exact."""
        for p in self.m.pairs:
            if p.copy_1[1] - p.copy_1[0] != p.copy_2[1] - p.copy_2[0]:
                with self.subTest(locus=p.locus):
                    self.assertFalse(p.exact,
                                     f"{p.locus} copies differ in length but the pair "
                                     f"is not marked partial")
                    self.assertIn("partial", p.note)

    def test_the_two_copies_do_not_overlap(self):
        for p in self.m.pairs:
            with self.subTest(locus=p.locus):
                self.assertTrue(p.copy_1[1] < p.copy_2[0] or p.copy_2[1] < p.copy_1[0],
                                f"{p.locus} copies overlap, so they are not two copies")

    def test_lat_is_paired(self):
        """The roadmap claimed LAT was paired while only its TATA boxes were.

        Containment-based matching could never find it: LAT's long-repeat copy runs
        778 bp past the end of IRL.
        """
        lat = [p for p in self.m.pairs if p.locus == "LAT"]
        self.assertTrue(lat, "LAT is not paired")
        self.assertTrue(any(p.feature_type == "gene" for p in lat),
                        f"only {[p.feature_type for p in lat]} paired, not the gene")

    def test_known_diploid_genes_are_paired(self):
        """RL1 and RL2 lie in the long repeat and must appear in both copies."""
        loci = {p.locus for p in self.m.pairs}
        for gene in ("RL1", "RL2"):
            self.assertTrue(any(gene in locus for locus in loci),
                            f"{gene} is in the long repeat and should be paired")


class TheQuestionItAnswers(unittest.TestCase):
    """Criterion: a structured answer, and the counterpart named."""

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)

    def test_a_unique_position_is_well_defined(self):
        s = self.m.classify(60000)
        self.assertTrue(s.well_defined)
        self.assertEqual(s.alternatives, [])

    def test_a_position_in_a_repeat_is_not_and_names_its_counterpart(self):
        s = self.m.classify(900)
        self.assertFalse(s.well_defined)
        self.assertEqual(s.alternatives, [self.m.counterparts(900)[0]])
        self.assertIn("TRL", s.reason)

    def test_counterparts_round_trip(self):
        for pos in (900, 5000, 9000, 118000, 121000, 130000, 148000):
            with self.subTest(position=pos):
                alts = self.m.counterparts(pos)
                self.assertTrue(alts, f"no counterpart for {pos}")
                for alt in alts:
                    self.assertIn(pos, self.m.counterparts(alt),
                                  f"{pos} -> {alt} does not come back")

    def test_counterparts_land_inside_the_other_copy(self):
        for a, b in INVERTED_PAIRS:
            ra = next(x for x in self.m.regions
                      if x.name.startswith(a) and x.kind == "repeat_region")
            rb = next(x for x in self.m.regions
                      if x.name.startswith(b) and x.kind == "repeat_region")
            for pos in (ra.start, ra.start + ra.length // 2, ra.end):
                with self.subTest(pair=f"{a}/{b}", position=pos):
                    alts = self.m.counterparts(pos)
                    self.assertTrue(
                        any(rb.start <= alt <= rb.end for alt in alts),
                        f"{pos} in {a} gives {alts}, none inside {b}")

    def test_the_repeat_junction_has_a_counterpart_in_each_repeat(self):
        """IRL and IRS overlap, so a position in the junction belongs to both."""
        junction = 126000
        alts = self.m.counterparts(junction)
        self.assertGreaterEqual(len(alts), 2,
                                f"{junction} is in both IRL and IRS and should have a "
                                f"counterpart in each, got {alts}")
        s = self.m.classify(junction)
        self.assertFalse(s.well_defined)
        self.assertIn("overlap", s.reason)

    def test_the_lat_tata_box_is_flagged(self):
        """118,777 is the LAT TATA box, inside IRL. It must not look unambiguous."""
        s = self.m.classify(118777)
        self.assertFalse(s.well_defined)
        self.assertTrue(any("LAT" in locus for locus in s.loci),
                        f"LAT not named among {s.loci}")

    def test_the_status_explains_itself(self):
        s = self.m.classify(900)
        self.assertIn("two copies", str(s))
        self.assertIn(str(s.alternatives[0]), str(s))


class CoversTheDefectsReviewFound(unittest.TestCase):
    """Criterion: the three M3 defects are demonstrably covered, not just claimed."""

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)

    def test_the_cross_repeat_coordinate_is_inside_the_map(self):
        """A position 900 bp into OP297884.1 resolved to 151,879, in TRS.

        The lift-over now declines it. This asserts the complementary fact: that
        151,879 is a position the map knows is not well defined, so the mask would have
        caught it even if the chaining had not.
        """
        s = self.m.classify(151879)
        self.assertFalse(s.well_defined,
                         "151,879 is inside TRS and must not be reported as unique")
        self.assertIn("TRS", s.reason)

    def test_the_lift_over_still_declines_it(self):
        if "OP297884.1" in Crosswalk.available(CROSSWALK):
            self.assertIsNone(lift(900, "OP297884.1", REFERENCE, CROSSWALK))

    def test_the_segmental_isomer_genome_is_mapped(self):
        """ON783215.1 has its long segment forward and its short segment reversed."""
        if "ON783215.1" not in Crosswalk.available(CROSSWALK):
            self.skipTest("ON783215.1 not in the pinned genome set")
        cw = Crosswalk.load("ON783215.1", CROSSWALK)
        self.assertGreater(cw.mapped_bases, 50000,
                           "the segmental isomer should still map substantially")

    def test_every_crosswalk_block_boundary_near_a_repeat_is_declarable(self):
        """Any reference position the crosswalk emits can be classified, not crash."""
        checked = 0
        for acc in Crosswalk.available(CROSSWALK)[:25]:
            for b in Crosswalk.load(acc, CROSSWALK).blocks[:5]:
                pos = b.target_start + 1
                s = self.m.classify(pos)
                self.assertIsInstance(s.well_defined, bool)
                checked += 1
        self.assertGreater(checked, 50)


class WorkedExample(unittest.TestCase):
    """Criterion: a tool declines and names the alternatives it declined between."""

    def test_a_tool_can_decline_with_reasons(self):
        m = RepeatMap.load(DERIVED)

        def report_position(position: int) -> str:
            s = m.classify(position)
            if s.well_defined:
                return f"position {position} is unique"
            alts = ", ".join(str(a) for a in s.alternatives)
            loci = ", ".join(s.loci) if s.loci else "no annotated locus"
            return (f"declining position {position}: {s.reason}. "
                    f"It could equally be {alts} ({loci})")

        unique = report_position(60000)
        declined = report_position(900)
        self.assertIn("unique", unique)
        self.assertIn("declining", declined)
        self.assertIn("could equally be", declined)
        self.assertIn("RL1", declined)




class CorrectionsFromReview(unittest.TestCase):
    """Guards the seven defects the M4 correctness review found."""

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)

    def test_repeats_carrying_a_gene_qualifier_are_not_dropped(self):
        """Three repeat_region features were skipped because M2 reclassed them."""
        for start in (988, 2466, 132391):
            with self.subTest(start=start):
                self.assertTrue(any(r.start == start for r in self.m.regions),
                                f"the repeat at {start} is missing from the map")

    def test_the_terminal_a_sequence_is_modelled_as_a_direct_repeat(self):
        """1..399 is identical to 151,824..152,222 in the SAME orientation.

        Only the inverted relation was modelled, so 798 positions had an incomplete
        list of alternatives.
        """
        s = self.m.classify(200)
        self.assertFalse(s.well_defined)
        self.assertIn(152023, s.alternatives,
                      f"the direct 'a' sequence counterpart is missing from "
                      f"{s.alternatives}")
        back = self.m.classify(152023)
        self.assertIn(200, back.alternatives)

    def test_out_of_range_positions_are_not_reported_as_well_defined(self):
        for pos in (0, -5, 152223, 1000000):
            with self.subTest(position=pos):
                s = self.m.classify(pos)
                self.assertFalse(s.well_defined)
                self.assertIn("outside the reference", s.reason)

    def test_no_status_claims_genome_wide_uniqueness(self):
        """The map describes annotated repeats and cannot certify uniqueness.

        The reason string used to assert "the position is still unique in the genome",
        which is false for every tandem reiteration set.
        """
        for pos in (117160, 117159, 71700, 60000, 1000):
            with self.subTest(position=pos):
                self.assertNotIn("unique in the genome",
                                 self.m.classify(pos).reason)

    def test_pairing_includes_loci_that_are_themselves_repeats(self):
        """Excluding every repeat_region from pairing lost nine genuine pairs."""
        self.assertGreaterEqual(len(self.m.pairs), 20,
                                f"only {len(self.m.pairs)} pairs derived")

    def test_the_lift_over_actually_consults_the_map(self):
        """A stated M4 criterion that was marked met while being false.

        `liftover.py` did not import RepeatMap at all. `lift_detailed` now reports
        whether the lifted coordinate identifies a unique locus.
        """
        from hsv1coords import lift_detailed
        accs = Crosswalk.available(CROSSWALK)
        if "JQ673480.1" not in accs or "JQ780693.1" not in accs:
            self.skipTest("the KOS pair is not in the pinned genome set")
        inside = lift_detailed(118777, "JQ673480.1", "JQ780693.1", CROSSWALK)
        self.assertIsNotNone(inside.position)
        self.assertFalse(inside.well_defined, "a position in IRL must be flagged")
        self.assertTrue(inside.alternatives)
        self.assertIn("NOT well defined", str(inside))

        outside = lift_detailed(60000, "JQ673480.1", "JQ780693.1", CROSSWALK)
        self.assertTrue(outside.well_defined)
        self.assertEqual(outside.alternatives, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
