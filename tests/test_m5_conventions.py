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

    def test_at_least_one_convention_is_unverified_and_refused(self):
        """Guards the rule itself: silently dropping the unverified one would pass
        every other test in this file while losing the honesty the rule exists for."""
        unverified = [c for c in self.cs.conventions if not c.verified]
        self.assertTrue(unverified,
                        "szpara2014_trimmed is paywalled and must remain declared "
                        "and refused, not quietly removed")


class RefusesWhatItCannotVerify(unittest.TestCase):
    """Criterion: the resolver declines to apply an unverified convention."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_usable_refuses(self):
        for c in self.cs.conventions:
            if c.verified:
                continue
            with self.subTest(convention=c.name):
                with self.assertRaises(UnverifiedConvention):
                    self.cs.usable(1000, c.name)

    def test_convert_refuses_in_either_direction(self):
        for c in self.cs.conventions:
            if c.verified:
                continue
            with self.subTest(convention=c.name):
                with self.assertRaises(UnverifiedConvention):
                    self.cs.convert(1000, "full_length", c.name)
                with self.assertRaises(UnverifiedConvention):
                    self.cs.convert(1000, c.name, "full_length")

    def test_the_refusal_names_the_source(self):
        try:
            self.cs.usable(1000, "szpara2014_trimmed")
        except UnverifiedConvention as e:
            self.assertIn("Szpara", str(e))
        else:
            self.fail("expected a refusal")


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

    def test_the_unverified_convention_is_reported_as_unanswerable(self):
        self.assertIn("szpara2014_trimmed", self.cs.report(60000).refused)

    def test_comparable_everywhere_accounts_for_the_refused(self):
        """A position cannot be called comparable everywhere while one convention
        cannot be evaluated at all."""
        self.assertFalse(self.cs.report(60000).comparable_everywhere)


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
