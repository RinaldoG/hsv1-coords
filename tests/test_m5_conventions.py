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

    def test_whisnant_masks_the_two_terminal_repeats(self):
        """TRL and TRS, not TRL and IRL.

        The first encoding masked both copies of the long repeat and neither of the
        short, which is impossible: the paper restores read mappings from the unmasked
        copy, and that encoding leaves none.
        """
        c = self.cs["whisnant2020_masked"]
        self.assertEqual(sorted(c.intervals), [(1, 9213), (145590, 152222)])

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
        self.assertIn("synthetic_unverified",
                      [n for n, _ in cs.report(60000).refused])


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

    def test_nothing_is_refused_once_the_version_crosswalk_exists(self):
        """Morse was refused until M6. Now it is translated, not guessed at."""
        self.assertEqual(self.cs.report(60000).refused, [])
        self.assertEqual([c.name for c in self.cs.other_version], ["morse2017_nr"])


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

    def test_szpara_and_whisnant_act_on_the_same_regions_differently(self):
        """Both touch TRL and TRS. One removes them, one masks and restores them.

        Same intervals, opposite consequences for a coordinate, which is the clearest
        illustration of why the intervals alone are not the answer.
        """
        sz, wh = self.cs["szpara2014_trimmed"], self.cs["whisnant2020_masked"]
        self.assertEqual(sorted(sz.intervals), sorted(wh.intervals))
        self.assertTrue(sz.changes_coordinates)
        self.assertFalse(wh.changes_coordinates)
        self.assertIsNone(sz.to_convention(5000))
        ok, _ = self.cs.usable(5000, "whisnant2020_masked")
        self.assertTrue(ok)

    def test_removed_positions_have_no_coordinate(self):
        sz = self.cs["szpara2014_trimmed"]
        self.assertIsNone(sz.to_convention(5000), "TRL is removed by Szpara")
        self.assertIsNone(sz.to_convention(150000), "TRS is removed by Szpara")

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
                           "describe it, even before the reference-version problem")

    def test_comparing_the_two_collapsed_references_needs_the_version_crosswalk(self):
        """Szpara is defined on NC_001806.2 and Morse on NC_001806.1.

        Before M6 this comparison was refused outright. An earlier version of this file
        instead computed a -9,213 difference and reported both as 136,376 bp, which
        treated Morse's intervals as though they were .2 coordinates: the resolver
        committing the error it exists to prevent. It is now translated through the
        version crosswalk, and the translation is not a constant shift.
        """
        from hsv1coords import VersionCrosswalk
        v = VersionCrosswalk.load()
        offsets = {v.to_old(p) - p for p in (1000, 60000, 120000, 150000)
                   if v.to_old(p) is not None}
        self.assertGreater(len(offsets), 1,
                           "if the version offset were constant, the earlier figures "
                           "would have been merely imprecise rather than wrong")

    def test_morses_own_stated_length_is_not_reproducible_here(self):
        """Morse states 136,770 bp. Applying their intervals in .2 gives a different
        number, which is the symptom of the version mismatch, recorded rather than
        reconciled."""
        mo = self.cs["morse2017_nr"]
        self.assertNotEqual(mo.length(), 136770,
                            "if this ever matches, a version crosswalk has appeared "
                            "and this row should be made applicable")

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

    def test_the_lat_tata_box_differs_between_conventions(self):
        cs = ConventionSet.load()
        pos = 118777
        self.assertEqual(cs["szpara2014_trimmed"].to_convention(pos), 109564)
        ok_l, _ = cs.usable(pos, "lassalle2020_excluded")
        self.assertFalse(ok_l)
        ok_w, _ = cs.usable(pos, "whisnant2020_masked")
        self.assertTrue(ok_w)

    def test_morse_removes_the_lat_tata_box(self):
        """The fourth answer, available only once the version crosswalk existed."""
        cs = ConventionSet.load()
        ok, why = cs.usable(118777, "morse2017_nr")
        self.assertFalse(ok)
        self.assertIn("117,161-125,974", why)
        self.assertIn("NC_001806.1", why)


class VersionCrosswalkGate(unittest.TestCase):
    """M6: the two reference versions, and the control proving the gate can fail.

    PITFALLS.md rule 2: every gate ships with a control. Here the control is a
    deliberately shifted lookup that must disagree with the sequence, proving the
    sequence check has power.
    """

    @classmethod
    def setUpClass(cls):
        from hsv1coords import VersionCrosswalk
        cls.v = VersionCrosswalk.load()
        import re as _re
        def seq(p):
            t = (ROOT / "data" / "raw" / p).read_text()
            return "".join(_re.findall(r"[acgtnACGTN]",
                           t.split("ORIGIN")[1].split("//")[0].replace("\n", ""))).upper()
        cls.old = seq("ref_NC_001806.1.gb")
        cls.new = seq("ref_NC_001806.2.gb")

    def test_the_versions_differ_as_published(self):
        self.assertEqual(len(self.old), 152261)
        self.assertEqual(len(self.new), 152222)

    def test_the_offset_is_not_constant(self):
        """The entire reason this crosswalk has to exist."""
        self.assertGreater(len(self.v.distinct_offsets), 1,
                           "a single offset would make the crosswalk unnecessary")

    def test_round_trip(self):
        failures = []
        for pos in range(1, len(self.old), 97):
            there = self.v.to_new(pos)
            if there is None:
                continue
            if self.v.to_old(there) != pos:
                failures.append((pos, there, self.v.to_old(there)))
        self.assertEqual(failures, [], f"{len(failures)} round trips failed")

    def test_mapped_bases_agree_in_sequence(self):
        """A coordinate map is only right if the bases match."""
        same = total = 0
        for pos in range(1, len(self.old), 29):
            there = self.v.to_new(pos)
            if there is None:
                continue
            a, b = self.old[pos - 1], self.new[there - 1]
            if a not in "ACGT" or b not in "ACGT":
                continue
            total += 1
            same += a == b
        self.assertGreater(total, 3000)
        rate = 100 * same / total
        self.assertGreaterEqual(rate, 99.0,
                                f"only {rate:.2f}% of {total} mapped bases agree")

    def test_a_deliberate_offset_fails_the_same_check(self):
        """The control. Shift by one and agreement must collapse."""
        same = total = 0
        for pos in range(1, len(self.old), 29):
            there = self.v.to_new(pos)
            if there is None or there >= len(self.new):
                continue
            a, b = self.old[pos - 1], self.new[there]       # deliberate +1
            if a not in "ACGT" or b not in "ACGT":
                continue
            total += 1
            same += a == b
        rate = 100 * same / total
        self.assertLess(rate, 60.0,
                        f"a one-base offset still agreed {rate:.2f}% of the time, so "
                        f"the sequence check cannot detect a frame error")

    def test_morse_is_now_usable_through_the_crosswalk(self):
        cs = ConventionSet.load()
        refused = [n for n, _ in cs.report(118777).refused]
        self.assertNotIn("morse2017_nr", refused,
                         "the version crosswalk should unblock it")
        ok, why = cs.usable(118777, "morse2017_nr")
        self.assertFalse(ok, "Morse removes IRL, so the LAT TATA box is not in it")
        self.assertIn("NC_001806.1", why,
                      "the answer must state which version it was evaluated in")
        self.assertIn("offset", why, "and the offset applied at that position")
