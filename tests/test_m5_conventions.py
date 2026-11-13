"""The M5 gate, written before the code, per ROADMAP.md."""
from __future__ import annotations

import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import ConventionSet, RepeatMap, UnverifiedConvention  # noqa: E402
from hsv1coords.conventions import GENOME_BP  # noqa: E402

DERIVED = ROOT / "data" / "derived"


class EveryDefinitionIsSourced(unittest.TestCase):
    """Criterion: quoted from its own primary source, or marked unverified."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_every_convention_names_a_source(self):
        for c in self.cs.conventions:
            with self.subTest(convention=c.name):
                self.assertTrue(c.source.strip(), f"{c.name} names no source")

    def test_every_verified_convention_carries_a_quote(self):
        for c in self.cs.verified:
            with self.subTest(convention=c.name):
                self.assertTrue(c.quote.strip())
                self.assertNotIn("NOT VERIFIED", c.quote)

    def test_every_verified_source_states_when_it_was_read(self):
        for c in self.cs.verified:
            with self.subTest(convention=c.name):
                self.assertIn("2026", c.source,
                              f"{c.name} does not say when its source was read")

    def test_unverified_conventions_say_why(self):
        for c in self.cs.conventions:
            if c.verified:
                continue
            with self.subTest(convention=c.name):
                self.assertIn("NOT VERIFIED", c.quote)
                self.assertTrue(len(c.source) > 40,
                                "an unverified convention must say what was tried")

    def test_every_convention_is_now_verified(self):
        """All four were read from their own sources in the end.

        szpara2014_trimmed was marked unverified on the strength of one failed API
        call. Review found the paper free on the publisher's own site. The lesson kept
        here is that an unverified row must name what was tried, so the next reader can
        try something else, and the refusal machinery is tested separately below
        against a synthetic row rather than by keeping a real one unverified.
        """
        unverified = [c.name for c in self.cs.conventions if not c.verified]
        self.assertEqual(unverified, [],
                         f"still unverified: {unverified}. If a definition genuinely "
                         f"cannot be read, that is acceptable, but the source field "
                         f"must say what was attempted")


class RefusesWhatItCannotVerify(unittest.TestCase):
    """Criterion: the resolver declines to apply an unverified convention."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    @staticmethod
    def _with_synthetic_unverified():
        """A set containing one deliberately unverified row, to test the refusal."""
        from hsv1coords.conventions import Convention
        base = ConventionSet.load()
        fake = Convention(
            name="synthetic_unverified", kind="trim", intervals=[(1, 100)],
            reference="NC_001806.2", verified=False,
            source="Synthetic row used only to test the refusal path; no source read",
            quote="NOT VERIFIED")
        return ConventionSet(base.conventions + [fake])

    def test_usable_refuses_an_unverified_convention(self):
        cs = self._with_synthetic_unverified()
        with self.assertRaises(UnverifiedConvention):
            cs.usable(1000, "synthetic_unverified")

    def test_convert_refuses_in_either_direction(self):
        cs = self._with_synthetic_unverified()
        with self.assertRaises(UnverifiedConvention):
            cs.convert(1000, "full_length", "synthetic_unverified")
        with self.assertRaises(UnverifiedConvention):
            cs.convert(1000, "synthetic_unverified", "full_length")

    def test_the_refusal_names_the_source(self):
        cs = self._with_synthetic_unverified()
        try:
            cs.usable(1000, "synthetic_unverified")
        except UnverifiedConvention as e:
            self.assertIn("no source read", str(e))
        else:
            self.fail("expected a refusal")

    def test_report_marks_an_unverified_convention_unanswerable(self):
        cs = self._with_synthetic_unverified()
        self.assertIn("synthetic_unverified", cs.report(60000).refused)


class DroppingDataIsNotMoving(unittest.TestCase):
    """Criterion: a convention that only drops data leaves coordinates unchanged,
    and 'unusable' is a different answer from 'moved'."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_exclusion_does_not_shift_coordinates(self):
        for c in self.cs.verified:
            with self.subTest(convention=c.name):
                if c.kind in ("exclude", "mask_then_restore", "none"):
                    self.assertFalse(c.changes_coordinates)

    def test_a_retained_position_converts_to_itself(self):
        self.assertEqual(
            self.cs.convert(60000, "full_length", "lassalle2020_excluded"), 60000)

    def test_an_excluded_position_converts_to_nothing(self):
        self.assertIsNone(
            self.cs.convert(5000, "full_length", "lassalle2020_excluded"))

    def test_usable_and_convert_disagree_in_the_right_way(self):
        """An excluded position is not usable, and has no equivalent. A masked and
        restored one is usable, and is its own equivalent."""
        ok, _ = self.cs.usable(5000, "lassalle2020_excluded")
        self.assertFalse(ok)
        ok, why = self.cs.usable(5000, "whisnant2020_masked")
        self.assertTrue(ok, "masking then restoring retains the position")
        self.assertIn("restored", why)


class RoundTrip(unittest.TestCase):
    """Criterion: a retained position converts out and back to itself."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_round_trip_for_every_retained_position(self):
        failures = []
        for c in self.cs.verified:
            for pos in range(1, GENOME_BP, 977):
                there = self.cs.convert(pos, "full_length", c.name)
                if there is None:
                    continue
                back = self.cs.convert(there, c.name, "full_length")
                if back != pos:
                    failures.append((c.name, pos, there, back))
        self.assertEqual(failures, [], f"{len(failures)} round trips failed")


class TheReportAnswersTheQuestion(unittest.TestCase):
    """Criterion: states which conventions a position is usable under."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_a_unique_region_position_is_comparable_everywhere_verified(self):
        rep = self.cs.report(60000)
        self.assertIn("lassalle2020_excluded", rep.usable_under)
        self.assertIn("whisnant2020_masked", rep.usable_under)
        self.assertEqual(rep.unusable_under, [])

    def test_a_repeat_position_is_not(self):
        rep = self.cs.report(5000)
        self.assertTrue(any(n == "lassalle2020_excluded"
                            for n, _ in rep.unusable_under))

    def test_nothing_is_refused_now_that_all_four_are_verified(self):
        self.assertEqual(self.cs.report(60000).refused, [])


class WorkedExample(unittest.TestCase):
    """Criterion: two papers' coordinates for one locus reconciled, or shown not to be."""

    def test_the_lat_tata_box_across_two_papers(self):
        cs = ConventionSet.load()
        m = RepeatMap.load(DERIVED)
        lat_tata = 118777

        usable_whisnant, _ = cs.usable(lat_tata, "whisnant2020_masked")
        usable_lassalle, why = cs.usable(lat_tata, "lassalle2020_excluded")

        self.assertTrue(usable_whisnant)
        self.assertFalse(usable_lassalle)
        self.assertIn("excluded", why)

        # And the repeat map adds the second half of the problem: even where the
        # position is retained, it does not name a unique locus.
        status = m.classify(lat_tata)
        self.assertFalse(status.well_defined)
        self.assertTrue(status.alternatives)

    def test_the_example_is_stateable_in_one_sentence(self):
        cs = ConventionSet.load()
        rep = cs.report(118777)
        text = str(rep)
        self.assertIn("usable under whisnant2020_masked", text)
        self.assertIn("NOT usable under lassalle2020_excluded", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)


class CollapsedReferences(unittest.TestCase):
    """Criterion three: a convention that removes sequence shifts coordinates, and the
    shift is computed rather than assumed constant.

    Marked met before it was implemented, which is the same error M4 was corrected for.
    """

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_both_collapsed_references_are_declared_shifting(self):
        for name in ("szpara2014_trimmed", "morse2017_nr"):
            with self.subTest(convention=name):
                self.assertTrue(self.cs[name].changes_coordinates)

    def test_removed_positions_have_no_coordinate(self):
        sz = self.cs["szpara2014_trimmed"]
        self.assertIsNone(sz.to_convention(5000), "TRL is removed by Szpara")
        mo = self.cs["morse2017_nr"]
        self.assertIsNone(mo.to_convention(120000), "IRL is removed by Morse")

    def test_round_trip_through_a_collapsed_reference(self):
        failures = []
        for name in ("szpara2014_trimmed", "morse2017_nr"):
            c = self.cs[name]
            for pos in range(1, GENOME_BP + 1, 311):
                there = c.to_convention(pos)
                if there is None:
                    continue
                if c.to_reference(there) != pos:
                    failures.append((name, pos, there, c.to_reference(there)))
        self.assertEqual(failures, [], f"{len(failures)} round trips failed")

    def test_the_collapsed_length_is_the_retained_length(self):
        for name in ("szpara2014_trimmed", "morse2017_nr"):
            c = self.cs[name]
            with self.subTest(convention=name):
                self.assertEqual(c.length(), GENOME_BP - c.affected_bases)

    def test_the_offset_is_not_constant(self):
        """The whole reason this needs computing rather than quoting one number."""
        sz = self.cs["szpara2014_trimmed"]
        offsets = {sz.to_convention(p) - p
                   for p in (10000, 60000, 140000)
                   if sz.to_convention(p) is not None}
        self.assertEqual(len(offsets), 1,
                         "Szpara removes only a leading block, so its offset is "
                         "constant across what it retains")
        mo = self.cs["morse2017_nr"]
        mo_offsets = {mo.to_convention(p) - p
                      for p in (10000, 60000, 140000)
                      if mo.to_convention(p) is not None}
        self.assertGreater(len(mo_offsets), 1,
                           "Morse removes an internal block, so a single offset cannot "
                           "describe it")

    def test_the_two_collapsed_references_disagree(self):
        """Same length, different content: the sharpest form of the problem."""
        sz, mo = self.cs["szpara2014_trimmed"], self.cs["morse2017_nr"]
        self.assertEqual(sz.length(), mo.length(),
                         "both collapse to the same length, which is what makes a "
                         "bare number look safe")
        disagreements = self.cs.offsets_between("szpara2014_trimmed", "morse2017_nr",
                                                step=5000)
        self.assertTrue(disagreements,
                        "two references removing different blocks must disagree")

    def test_convert_between_two_collapsed_references(self):
        cs = self.cs
        there = cs.convert(60000, "full_length", "szpara2014_trimmed")
        self.assertEqual(there, 50787)
        back = cs.convert(there, "szpara2014_trimmed", "full_length")
        self.assertEqual(back, 60000)


class RepeatCaveats(unittest.TestCase):
    """A convention retaining a position does not make the position meaningful."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_the_gap_in_the_lassalle_list_is_caveated(self):
        """Lassalle's exclusion stops at 150,962 while TRS runs to 152,222.

        Without the caveat the resolver calls a terminal-repeat position comparable.
        """
        rep = self.cs.report(152000)
        self.assertIn("lassalle2020_excluded", rep.usable_under)
        self.assertTrue(rep.caveats, "a position in TRS must carry a caveat")
        self.assertIn("not well defined", rep.caveats[0])

    def test_a_unique_position_carries_no_caveat(self):
        self.assertEqual(self.cs.report(60000).caveats, [])

    def test_comparable_everywhere_ignores_the_refused(self):
        """Counting refused conventions made this False at every position."""
        self.assertTrue(self.cs.report(60000).comparable_everywhere)


class FourAnswersForOneLocus(unittest.TestCase):
    """The worked example, in its strongest form."""

    def test_the_lat_tata_box_differs_under_every_convention(self):
        cs = ConventionSet.load()
        pos = 118777
        self.assertEqual(cs["szpara2014_trimmed"].to_convention(pos), 109564)
        self.assertIsNone(cs["morse2017_nr"].to_convention(pos),
                          "Morse removes IRL, so the LAT TATA box is not in it at all")
        ok_l, _ = cs.usable(pos, "lassalle2020_excluded")
        self.assertFalse(ok_l)
        ok_w, _ = cs.usable(pos, "whisnant2020_masked")
        self.assertTrue(ok_w)
