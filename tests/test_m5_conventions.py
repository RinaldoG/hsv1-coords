"""The M5 gate, written before the code, per ROADMAP.md."""
from __future__ import annotations

import ast
import csv
import fnmatch
import pathlib
import re
import shutil
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import ConventionSet, RepeatMap, UnverifiedConvention  # noqa: E402
from hsv1coords.conventions import GENOME_BP  # noqa: E402

DERIVED = ROOT / "data" / "derived"
CROSSWALK_DIR = DERIVED / "crosswalk"


class EveryDefinitionIsSourced(unittest.TestCase):
    """Criterion: quoted from its own primary source, or marked unverified."""

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_every_convention_names_a_source(self):
        for c in self.cs.conventions:
            with self.subTest(convention=c.name):
                self.assertTrue(c.source.strip(), f"{c.name} names no source")

    def test_every_interval_is_tied_to_its_own_source(self):
        """Round 7: no test tied any interval to the quote it came from.

        Lassalle's three intervals were anchored nowhere at all, and Szpara's were
        anchored to Whisnant's, which were anchored to a literal typed into this file.
        A literal in a test is not a source: it travels with the mistake.

        So each interval is checked against something outside this repository's own
        prose. Where the quote states coordinates, they are parsed from the quote.
        Where it does not, the interval must match a repeat_region in the reference
        record's feature table, which is the primary source for those extents.
        """
        import re
        ref = ROOT / "data" / "raw" / "ref_NC_001806.2.gb"
        if not ref.exists():
            self.skipTest("data/raw/ref_NC_001806.2.gb is fetched, not committed; "
                          "run `make inputs` first")
        body = ref.read_text()
        annotated = set()
        for loc, quals in re.findall(
                r"^     repeat_region\s+(\S+)\n((?:\s+/.*\n)+)", body, re.M):
            m = re.match(r"(\d+)\.\.(\d+)$", loc)
            if m and "inverted repeat flanking" in quals:
                annotated.add((int(m.group(1)), int(m.group(2))))
        self.assertTrue(annotated, "no inverted repeats found in the record")

        checked = 0
        for c in self.cs.verified:
            if not c.intervals:
                continue
            pairs = re.findall(r"(?<![\d,])([\d,]+)\s*-\s*([\d,]+)(?![\d,])",
                               c.quote)
            quoted = {(lo, hi) for lo, hi in
                      ((int(a.replace(",", "")), int(b.replace(",", "")))
                       for a, b in pairs)
                      if 0 < lo < hi <= GENOME_BP}
            for interval in c.intervals:
                with self.subTest(convention=c.name, interval=interval):
                    if quoted:
                        self.assertIn(
                            interval, quoted,
                            f"{c.name} states {interval}, which is not among the "
                            f"intervals its own quote gives ({sorted(quoted)})")
                    else:
                        self.assertIn(
                            interval, annotated,
                            f"{c.name}'s quote gives no coordinates, so {interval} "
                            f"must be an annotated inverted repeat in the reference "
                            f"record, and it is not: {sorted(annotated)}")
                    checked += 1
        self.assertGreaterEqual(checked, 5,
                                f"only {checked} intervals were tied to a source")

    def test_the_morse_edges_come_from_the_previous_records_feature_table(self):
        """Round 12: the one row that has been wrong three times was tied to nothing.

        The test above iterates `cs.verified`, which by construction excludes
        `morse2017_nr`, so its two intervals rested on a sum (136,770 bp, the length the
        paper states) plus one string assertion. This project's own note on that row
        says a sum cannot locate a boundary, which is exactly how the second wrong
        version of it passed: boundaries wrong at three of four edges, with the errors
        cancelling.

        Its intervals are in `NC_001806.1`, so the `NC_001806.2` extents the test above
        uses cannot check them. Each of the four edges is checked here against a number
        in `NC_001806.1`'s own feature table, which is the frame the row declares.
        """
        import re
        record = ROOT / "data" / "raw" / "ref_NC_001806.1.gb"
        if not record.exists():
            self.skipTest("ref_NC_001806.1.gb is fetched, not committed; "
                          "run `make inputs` first")
        body = record.read_text()
        extents = {}
        for loc, quals in re.findall(
                r"^     repeat_region\s+(\S+)\n((?:\s+/.*\n)+)", body, re.M):
            m = re.match(r"(\d+)\.\.(\d+)$", loc)
            note = re.search(r'/note="([^";]+)', quals)
            if m and note and "inverted repeat flanking" in quals:
                extents[note.group(1).strip()] = (int(m.group(1)), int(m.group(2)))
        for need in ("IRL", "IRS", "TRS"):
            self.assertIn(need, extents,
                          f"NC_001806.1 does not annotate {need} as an inverted repeat")
        length = int(re.search(r"^LOCUS\s+\S+\s+(\d+) bp", body, re.M).group(1))

        morse = self.cs["morse2017_nr"]
        self.assertEqual(morse.reference_version, "NC_001806.1",
                         "this test checks the row against NC_001806.1 because that is "
                         "the frame the row declares")
        expected = {
            # The paper's rule: cut at the RL-UL and RS-US junctions, so the retained
            # blocks run from the start of IRL to the base before IRS, and from the
            # start of TRS to the end of the record.
            (extents["IRL"][0], extents["IRS"][0] - 1): "IRL start to IRS start - 1",
            (extents["TRS"][0], length): "TRS start to the end of the record",
        }
        self.assertEqual(
            set(morse.intervals), set(expected),
            f"morse2017_nr states {sorted(morse.intervals)}; the record's own feature "
            f"table gives {sorted(expected)} "
            f"({', '.join(expected.values())})")
        # The intervals are the blocks REMOVED (IRL and TRS), so what the paper states
        # is the record's length less their extent.
        removed = sum(b - a + 1 for a, b in morse.intervals)
        self.assertEqual(
            length - removed, 136770,
            f"the record is {length:,} bp and the row removes {removed:,}, leaving "
            f"{length - removed:,}, where the paper states 136,770")

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
        return ConventionSet((*base.conventions, fake))

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
        self.assertFalse(self.cs["szpara2014_trimmed"].is_split)
        self.assertTrue(self.cs["morse2017_nr"].is_split)

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

    def test_morse_is_a_four_record_reference_not_one_sequence(self):
        """Its coordinates are pairs like UL:50788, confirmed from its deposited VCFs.

        An earlier version modelled it as a single collapsed sequence, copied from
        Szpara, and returned integers that appear in no Morse file.
        """
        mo = self.cs["morse2017_nr"]
        self.assertTrue(mo.is_split)
        self.assertEqual([n for n, _, _ in mo.records], ["RL", "UL", "RS", "US"])
        with self.assertRaises(ValueError):
            mo.to_convention(60000)
        self.assertEqual(self.cs.coordinate_in(60000, "morse2017_nr"), "UL:50787")

    def test_morse_records_sum_to_the_published_length(self):
        """136,770 bp, from the paper's stated cutting rule with nothing fitted.

        The previous row reached this total with boundaries wrong at three of four
        edges, because the errors cancelled. A sum cannot locate a boundary, so this
        is a consistency check and the boundaries stand on the rule, not on the total.
        """
        mo = self.cs["morse2017_nr"]
        self.assertEqual(sum(b - a + 1 for _, a, b in mo.records), 136770)

    def test_convert_translates_versions_like_usable_does(self):
        """Root cause 6b survived its own fix here: usable translated, convert did not."""
        for target in ("lassalle2020_excluded", "whisnant2020_masked",
                       "szpara2014_trimmed"):
            for pos in (1000, 60000, 130000, 152000):
                with self.subTest(target=target, position=pos):
                    ok, _ = self.cs.usable(pos, target)
                    got = self.cs.convert(pos, "full_length", target)
                    self.assertEqual(got is not None, ok,
                                     "usable and convert must not contradict")

    def test_convert_refuses_a_multi_record_convention(self):
        with self.assertRaises(ValueError):
            self.cs.convert(60000, "full_length", "morse2017_nr")

    def test_removed_positions_have_no_coordinate(self):
        sz = self.cs["szpara2014_trimmed"]
        self.assertIsNone(sz.to_convention(5000), "TRL is removed by Szpara")
        self.assertIsNone(sz.to_convention(150000), "TRS is removed by Szpara")

    def test_round_trip_through_a_collapsed_reference(self):
        failures = []
        for name in ("szpara2014_trimmed",):      # Morse is multi-record, tested above
            c = self.cs[name]
            for pos in range(1, GENOME_BP + 1, 311):
                there = c.to_convention(pos)
                if there is None:
                    continue
                if c.to_reference(there) != pos:
                    failures.append((name, pos, there, c.to_reference(there)))
        self.assertEqual(failures, [], f"{len(failures)} round trips failed")

    def test_the_collapsed_length_is_the_retained_length(self):
        """Measured against the length of the version each convention declares."""
        for name in ("szpara2014_trimmed",):
            c = self.cs[name]
            with self.subTest(convention=name):
                self.assertEqual(c.length(), c.genome_length - c.affected_bases)

    def test_each_convention_uses_its_own_versions_length(self):
        """Guards the frame error: the row declared .1 and was measured against .2."""
        self.assertEqual(self.cs["morse2017_nr"].genome_length, 152261)
        self.assertEqual(self.cs["szpara2014_trimmed"].genome_length, 152222)

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
        # Morse has no scalar offset at all: its coordinates are record-and-position
        # pairs, so asking for one is the category error this row used to commit.
        starts = {self.cs.coordinate_in(p, "morse2017_nr").split(":")[0]
                  for p in (10000, 60000, 130000, 140000)
                  if self.cs.coordinate_in(p, "morse2017_nr")}
        self.assertGreater(len(starts), 1,
                           "positions in different records cannot share one offset")

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

    def test_morses_stated_length_is_now_reproduced(self):
        """136,770 bp, reproduced once the intervals were put in the right frame.

        It did not reproduce while the row stored NC_001806.2 coordinates under a
        declared NC_001806.1 frame, and the 5 bp shortfall was exactly the version gap.
        Matching the stated length is what the reconstructed boundaries were chosen to
        do, so it is a consistency check and not independent evidence.
        """
        self.assertEqual(self.cs["morse2017_nr"].length(), 136770)

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
        self.assertIn("117,160-125,973", why,
                      "the interval must be quoted in NC_001806.1, the frame the row "
                      "declares, and at the boundaries the paper's own rule gives")
        self.assertIn("NC_001806.1", why)


class VersionCrosswalkGate(unittest.TestCase):
    """M6: the two reference versions, and the control proving the gate can fail.

    PITFALLS.md rule 2: every gate ships with a control. Here the control is a
    deliberately shifted lookup that must disagree with the sequence, proving the
    sequence check has power.
    """

    @classmethod
    def setUpClass(cls):
        import unittest as _ut
        for _need in ("ref_NC_001806.1.gb", "ref_NC_001806.2.gb"):
            if not (ROOT / "data" / "raw" / _need).exists():
                raise _ut.SkipTest(
                    f"{_need} is fetched, not committed. Run `make inputs` "
                    f"and `make genomes` first.")
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

    # --- the gates, and the mutants that prove they can fail -------------------
    #
    # PITFALLS.md root cause 3 records two gates in this repository that could not
    # fail. Review found this file contained a third: a round trip through the same
    # blocks is exact whatever the blocks say, and the "control" measured only that DNA
    # is not self-similar under a one-base shift, which is a property of the sequence
    # and passes on a random crosswalk. Deleting a whole block passed every gate.
    #
    # The gates below are therefore tested against deliberately broken crosswalks. A
    # gate that passes on a mutant is not a gate.

    # Five gates, because no one of them is sufficient. Review built 15 broken
    # crosswalks against an earlier version of this file and 12 passed everything.
    # Each gate below exists because a specific mutant defeats all the others.

    EXPECTED_BLOCKS = 22
    EXPECTED_COVERAGE = 152209

    @staticmethod
    def _mutants(blocks):
        """Deliberately broken crosswalks, each a realistic failure mode."""
        import random
        from hsv1coords.versions import VersionBlock, VersionCrosswalk
        rng = random.Random(7)
        m = {}

        def shift(i, d, label):
            b2 = list(blocks)
            b2[i] = VersionBlock(blocks[i].old_start, blocks[i].new_start + d,
                                 blocks[i].length)
            m[label] = VersionCrosswalk(b2)

        shift(len(blocks) // 2, 1, "mid block shifted 1")
        shift(len(blocks) // 2, 1000, "mid block shifted 1000")
        shift(0, 5, "first block shifted 5")
        shift(len(blocks) - 1, 40000, "last block pushed out of range")
        for i, b in enumerate(blocks):
            if b.length < 100:
                b2 = list(blocks)
                b2.pop(i)
                m["smallest block deleted"] = VersionCrosswalk(b2)
                break
        big = max(range(len(blocks)), key=lambda i: blocks[i].length)
        b2 = list(blocks)
        b2.pop(big)
        m["largest block deleted"] = VersionCrosswalk(b2)
        i = len(blocks) // 3
        b2 = list(blocks)
        b2[i] = VersionBlock(blocks[i].old_start, blocks[i].new_start,
                             blocks[i].length + 1500)
        m["block extended 1500"] = VersionCrosswalk(b2)
        b2 = list(blocks)
        b2[i] = VersionBlock(blocks[i].old_start, blocks[i + 1].new_start,
                             blocks[i].length)
        m["two source blocks onto one target span"] = VersionCrosswalk(b2)
        m["targets randomised"] = VersionCrosswalk(
            [VersionBlock(b.old_start, rng.randrange(0, 150000), b.length)
             for b in blocks])
        b2 = list(blocks)
        b2[1] = VersionBlock(blocks[1].old_start - 50, blocks[1].new_start,
                             blocks[1].length)
        m["boundary moved 50 in the query"] = VersionCrosswalk(b2)

        # Two mutants that preserve coverage, bounds and non-overlap, so only a gate
        # that looks AT the block edges can see them. Without these, the boundary
        # identity gate could be deleted with the suite still green.
        pairs = [(i, j) for i in range(len(blocks)) for j in range(i + 1, len(blocks))
                 if blocks[i].length == blocks[j].length]
        if pairs:
            i, j = pairs[0]
            b2 = list(blocks)
            b2[i] = VersionBlock(blocks[i].old_start, blocks[j].new_start,
                                 blocks[i].length)
            b2[j] = VersionBlock(blocks[j].old_start, blocks[i].new_start,
                                 blocks[j].length)
            m["two equal-length blocks swapped"] = VersionCrosswalk(b2)
        for i in range(len(blocks) - 1):
            a, b = blocks[i], blocks[i + 1]
            if a.length > 400 and b.length > 400:
                b2 = list(blocks)
                b2[i] = VersionBlock(a.old_start, a.new_start, a.length - 200)
                b2[i + 1] = VersionBlock(b.old_start - 200, b.new_start - 200,
                                         b.length + 200)
                m["boundary moved 200, coverage preserved"] = VersionCrosswalk(b2)
                break
        return m

    def _identity(self, cw, step=29) -> float:
        same = total = 0
        for pos in range(1, len(self.old), step):
            there = cw.to_new(pos)
            if there is None or there > len(self.new):
                continue
            a, b = self.old[pos - 1], self.new[there - 1]
            if a not in "ACGT" or b not in "ACGT":
                continue
            total += 1
            same += a == b
        return 100 * same / total if total else 0.0

    def _boundary_identity(self, cw) -> float:
        """Identity within 30 bases of a block edge, where a moved boundary shows.

        Genome-wide identity misses a boundary move, because the moved region is a
        small fraction of the whole. This is the gate that catches it.
        """
        same = total = 0
        for b in cw.blocks:
            edges = (list(range(max(1, b.old_start - 30), b.old_start + 31))
                     + list(range(max(1, b.old_start + b.length - 30),
                                  min(len(self.old), b.old_start + b.length + 31))))
            for pos in edges:
                there = cw.to_new(pos)
                if there is None or there > len(self.new):
                    continue
                a, c = self.old[pos - 1], self.new[there - 1]
                if a not in "ACGT" or c not in "ACGT":
                    continue
                total += 1
                same += a == c
        return 100 * same / total if total else 0.0

    @staticmethod
    def _in_bounds(cw, old_len, new_len) -> bool:
        return all(b.new_start + b.length <= new_len
                   and b.old_start + b.length <= old_len for b in cw.blocks)

    @staticmethod
    def _non_overlapping(cw) -> bool:
        for key in ("new_start", "old_start"):
            spans = sorted((getattr(b, key), getattr(b, key) + b.length)
                           for b in cw.blocks)
            if any(end > nxt for (_, end), (nxt, _) in zip(spans, spans[1:])):
                return False
        return True

    def _invariants(self, cw) -> dict[str, bool]:
        """Structural properties the builder must guarantee, whatever the data says.

        These are not gates in the mutation-testing sense, and trying to treat them
        as such is what made this file claim power it did not have. Overlapping blocks
        make the lookup itself ill-defined: a bisect would return an arbitrary one of
        two answers. Out-of-range blocks index past the sequence. Neither is a question
        about whether the alignment is biologically right, so neither needs a mutant to
        justify it. They are asserted directly on the real crosswalk instead.
        """
        return {
            "in bounds": self._in_bounds(cw, len(self.old), len(self.new)),
            "non-overlapping": self._non_overlapping(cw),
        }

    def _gates(self, cw) -> dict[str, bool]:
        """The empirical gates: claims about whether the mapping is correct.

        Genome-wide identity is deliberately NOT among them.

        With ungapped blocks, correct edges imply correct interiors, so a genome-wide
        identity check catches nothing the boundary check misses: measured across
        every mutant here, its catches are a strict subset. Keeping it as a fifth gate
        would have meant one gate that could be deleted with the suite still green,
        which is how review found three of the previous five. It survives below as a
        documented sanity check, not as a gate.
        """
        return {
            "boundary identity": self._boundary_identity(cw) >= 99.0,
            "exact coverage": sum(b.length for b in cw.blocks) == self.EXPECTED_COVERAGE,
        }

    def test_round_trip(self):
        failures = []
        for pos in range(1, len(self.old), 97):
            there = self.v.to_new(pos)
            if there is None:
                continue
            if self.v.to_old(there) != pos:
                failures.append((pos, there, self.v.to_old(there)))
        self.assertEqual(failures, [], f"{len(failures)} round trips failed")

    def test_the_round_trip_gate_is_nearly_powerless(self):
        """Measured, not asserted. It sees only mutants that make blocks overlap."""
        missed = []
        for name, mutant in self._mutants(self.v.blocks).items():
            ok = all(mutant.to_old(mutant.to_new(p)) == p
                     for p in range(1, len(self.old), 997)
                     if mutant.to_new(p) is not None)
            if ok:
                missed.append(name)
        for blind_spot in ("mid block shifted 1", "smallest block deleted"):
            self.assertIn(blind_spot, missed,
                          "the round trip cannot see this, which is why the gates "
                          "below exist")

    def test_the_true_crosswalk_passes_every_gate(self):
        failed = [g for g, ok in self._gates(self.v).items() if not ok]
        self.assertEqual(failed, [], f"the real crosswalk fails: {failed}")
        self.assertEqual(len(self.v.blocks), self.EXPECTED_BLOCKS)

    def test_the_structural_invariants_hold(self):
        """Asserted on the real crosswalk, because they are properties of the build."""
        broken = [k for k, ok in self._invariants(self.v).items() if not ok]
        self.assertEqual(broken, [], f"the builder produced blocks that are {broken}")

    def test_no_mutant_survives_the_gates_and_invariants_together(self):
        """The control, over both kinds of check."""
        survivors = []
        for name, mutant in self._mutants(self.v.blocks).items():
            checks = {**self._gates(mutant), **self._invariants(mutant)}
            if all(checks.values()):
                survivors.append(name)
        self.assertEqual(survivors, [],
                         f"these broken crosswalks passed everything: {survivors}")

    def test_each_gate_is_the_sole_catcher_of_something(self):
        """A gate no mutant needs is a gate that is not pulling its weight.

        Recorded so that removing one later shows up as a failure here rather than as
        a silent loss of power.
        """
        mutants = self._mutants(self.v.blocks)
        sole = {}
        for name, mutant in mutants.items():
            caught = [g for g, ok in self._gates(mutant).items() if not ok]
            if len(caught) == 1:
                sole.setdefault(caught[0], []).append(name)
        missing = [g for g in self._gates(self.v) if g not in sole]
        # Invariants are exempt by design: see `_invariants`.
        self.assertEqual(missing, [],
                         f"these gates catch nothing on their own, so each could be "
                         f"deleted with the suite still green: {missing}. Either add "
                         f"a mutant only that gate catches, or remove the gate.")

    def test_genome_wide_identity_is_subsumed_and_says_so(self):
        """Kept as a sanity check, and measured rather than assumed redundant."""
        self.assertGreaterEqual(self._identity(self.v), 99.5)
        for name, mutant in self._mutants(self.v.blocks).items():
            if self._identity(mutant) < 99.5:
                with self.subTest(mutant=name):
                    self.assertFalse(
                        all(self._gates(mutant).values()),
                        f"{name} is caught by genome-wide identity and by no gate, so "
                        f"identity is not subsumed after all and should be a gate")

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


class TheFixesAreReachable(unittest.TestCase):
    """Tests that fail if a previously fixed defect is reintroduced.

    Round 3 found that reverting `convert`'s version translation, re-breaking
    `report`'s warning handling, and deleting the Morse `quote` each left the whole
    suite green. A fix nothing can reach is not covered, however correct the code is.
    PITFALLS.md root cause: verifying a fix with inputs that cannot exercise it.
    """

    @staticmethod
    def _with_other_version_row():
        """A non-split convention on NC_001806.1, so `convert` must translate.

        The real conventions cannot exercise this: the only NC_001806.1 row is
        multi-record, and `convert` refuses those, so every other pair is .2 to .2
        where translation is a no-op.
        """
        from hsv1coords.conventions import Convention
        base = ConventionSet.load()
        row = Convention(
            name="synthetic_v1_exclude", kind="exclude",
            intervals=[(1, 1000)], reference="NC_001806.1", verified=True,
            source="Synthetic row, used only to exercise version translation",
            quote="synthetic", reference_version="NC_001806.1",
            applicable=False, boundary_uncertainty_bp=3)
        return ConventionSet((*base.conventions, row))

    def test_convert_actually_translates_between_versions(self):
        """Fails if `convert` stops calling the version crosswalk.

        The expected value is computed here from the crosswalk directly, NOT from
        `convert`. The first version of this test compared `convert` with itself,
        which is the tautology pattern in PITFALLS.md root cause 3, written into the
        fix for a different root cause. Re-breaking `convert` left it green.
        """
        from hsv1coords import VersionCrosswalk
        cs = self._with_other_version_row()
        v = VersionCrosswalk.load()
        checked = 0
        for v1_pos in (60000, 130000, 150000):
            expected = v.to_new(v1_pos)          # NC_001806.1 position into .2
            if expected is None or expected == v1_pos:
                continue                          # no detectable difference here
            with self.subTest(position=v1_pos, expected=expected):
                got = cs.convert(v1_pos, "synthetic_v1_exclude", "full_length")
                self.assertEqual(got, expected,
                                 f"a NC_001806.1 position {v1_pos} must arrive in "
                                 f"NC_001806.2 as {expected}, not {got}")
                checked += 1
        self.assertGreater(checked, 0,
                           "no position was checked where the two versions differ, "
                           "so this test had no power")

    def test_convert_and_usable_agree_on_an_other_version_row(self):
        cs = self._with_other_version_row()
        for pos in range(500, 152000, 4999):
            with self.subTest(position=pos):
                ok, _ = cs.usable(pos, "synthetic_v1_exclude")
                got = cs.convert(pos, "full_length", "synthetic_v1_exclude")
                self.assertEqual(got is not None, ok,
                                 "usable and convert must not contradict, which is "
                                 "exactly what they did before round 2")

    def test_report_keeps_a_warning_on_a_usable_answer(self):
        """Fails if `report` goes back to dropping warnings when the answer is usable.

        That was the original behaviour, and it discarded every warning on precisely
        the positions a reader would act on.
        """
        cs = self._with_other_version_row()
        found = None
        for pos in range(1, 4000):
            rep = cs.report(pos)
            if any("synthetic_v1_exclude" in c for c in rep.caveats):
                found = (pos, rep)
                break
        self.assertIsNotNone(found, "no position produced a warning to check")
        pos, rep = found
        self.assertTrue(any("unreliable" in c.lower() for c in rep.caveats))

    def test_morse_records_are_pinned_to_the_feature_table_not_to_a_sum(self):
        """Derives the four records from NC_001806.1's own annotation and compares.

        The length check alone cannot hold these boundaries: moving the RS/US edge
        5,000 bp with the total preserved leaves the whole suite green and makes
        `coordinate_in` return a position past the end of a record, which is the error
        this row was corrected for in round 2. PITFALLS.md says a sum cannot locate a
        boundary; that lesson was applied to the derivation and not to the test that
        pins it.

        The paper's rule: cut at the RL-UL and RS-US junctions.
        """
        import re as _re
        gb = ROOT / "data" / "raw" / "ref_NC_001806.1.gb"
        if not gb.exists():
            self.skipTest("ref_NC_001806.1.gb is fetched, not committed: "
                          "run `make inputs` first.")
        feats = {}
        for loc, quals in _re.findall(r"^     repeat_region\s+(\S+)\n((?:\s+/.*\n)+)",
                                      gb.read_text(), _re.M):
            note = _re.search(r'/note="([^"]+)"', quals)
            if not note or "reiteration" in note.group(1):
                continue
            m = _re.match(r"(\d+)\.\.(\d+)$", loc)
            if m and "inverted repeat flanking" in note.group(1):
                feats[note.group(1).split(";")[0]] = (int(m.group(1)), int(m.group(2)))
        for need in ("TRL", "IRL", "IRS", "TRS"):
            self.assertIn(need, feats, f"{need} not found in the .1 feature table")

        expected = {
            "RL": (feats["TRL"][0], feats["TRL"][1]),
            "UL": (feats["TRL"][1] + 1, feats["IRL"][0] - 1),
            "RS": (feats["IRS"][0], feats["IRS"][1]),
            "US": (feats["IRS"][1] + 1, feats["TRS"][0] - 1),
        }
        stored = {n: (a, b) for n, a, b in ConventionSet.load()["morse2017_nr"].records}
        self.assertEqual(stored, expected,
                         "the stored records do not match the paper's rule applied to "
                         "the .1 feature table")

    def test_a_coordinate_never_runs_past_the_end_of_its_record(self):
        """The symptom a sum-only check misses."""
        cs = ConventionSet.load()
        mo = cs["morse2017_nr"]
        lengths = {n: b - a + 1 for n, a, b in mo.records}
        for pos in range(1, 152222, 503):
            coord = cs.coordinate_in(pos, "morse2017_nr")
            if coord is None:
                continue
            rec, off = coord.split(":")
            with self.subTest(position=pos, coordinate=coord):
                self.assertLessEqual(int(off), lengths[rec],
                                     f"{coord} is past the end of {rec}, which is "
                                     f"{lengths[rec]:,} bp")

    def test_every_convention_row_is_individually_checked(self):
        """Fails if any row's provenance is deleted, including ones a filter excludes.

        The provenance tests iterate `cs.verified`, which excludes exactly the row
        that has been wrong three times. A filter that exists to exclude one row means
        that row needs its own named test.
        """
        cs = ConventionSet.load()
        by_name = {c.name: c for c in cs.conventions}
        self.assertIn("morse2017_nr", by_name)
        morse = by_name["morse2017_nr"]
        self.assertTrue(morse.quote.strip(), "the Morse quote is missing")
        self.assertIn("NR-reference", morse.quote)
        self.assertIn("152,261", morse.quote,
                      "the quote must carry the length that fixes its frame")
        self.assertIn("Reliable Detection", morse.source,
                      "the citation title was invented once and must stay checked")
        self.assertIn("2017", morse.source)
        self.assertIn("PMC5580483", morse.source)
        self.assertIn("November 2026", morse.source,
                      "every verified row records when its source was read")
        for c in cs.conventions:
            with self.subTest(convention=c.name):
                self.assertIn("2026", c.source,
                              f"{c.name} does not say when its source was read")
        for c in cs.conventions:
            with self.subTest(convention=c.name):
                self.assertTrue(c.quote.strip())
                self.assertTrue(c.source.strip())


class EveryTestFileRunsAllOfItself(unittest.TestCase):
    """Fails if a `__main__` block sits anywhere but the end of its file.

    In this file it sat at line 227 of 500, so running it directly executed 20 of its
    53 tests and printed OK. The 33 it hid were every Morse test, the whole
    version-gate apparatus and the class holding two previously closed blockers.
    `make test` uses discovery and saw them all, which is why nothing ever failed.

    Checked statically rather than by running the file, because a guard that invokes
    its own file recurses forever.
    """

    def test_no_test_class_is_defined_after_the_main_block(self):
        import re as _re
        for path in sorted((ROOT / "tests").glob("test_*.py")):
            text = path.read_text()
            m = _re.search(r'^if __name__ == "__main__":', text, _re.M)
            with self.subTest(file=path.name):
                # Round 11: this skipped a file with no `__main__` block, and a slice
                # edit of mine had deleted that block from two files along with a whole
                # test class. The absence of the block is the thing worth catching: it
                # is what a truncated file looks like.
                self.assertIsNotNone(
                    m,
                    f"{path.name} has no `__main__` block. Every other test file has "
                    f"one, so this file has probably been truncated by an edit")
            if not m:
                continue
            after = text[m.end():]
            orphans = _re.findall(r"^class (\w+)\(", after, _re.M)
            orphans += _re.findall(r"^    def (test_\w+)", after, _re.M)
            with self.subTest(file=path.name):
                self.assertEqual(
                    orphans, [],
                    f"{path.name} defines {orphans} after its `__main__` block, so "
                    f"running the file directly silently skips them")

    TESTS_PER_FILE = {
        "test_docs_match_artefacts": 16,
        "test_m1_dictionary": 48,
        "test_m2_annotation": 45,
        "test_m3_crosswalk": 53,
        "test_m4_repeat_map": 58,
        "test_m5_conventions": 119,
        "test_placeholder": 1,
        "test_readme_examples": 7,
    }
    """How many tests each file contributes to discovery.

    Round 12: the `__main__` guard above catches a file truncated at the end, which is
    what my round 11 slice edit did, and nothing else. Deleting the class that edit
    destroyed, `WhetherTheRepeatCopySurvives`, while leaving the `__main__` block in
    place, left the suite green at 265 tests instead of 269. Nothing in the repository
    counted tests, classes or methods, or put a floor under any of them, so a whole
    class could go missing in a commit and the only trace was the total scrolling past.

    These numbers rise whenever tests are added, which is ordinary and is why the
    failure says so. A number that falls is a test that stopped running, and that needs
    a reason in the commit message rather than a quiet edit here.
    """

    def test_a_checkout_holds_every_artefact_the_gates_skip_without(self):
        """Round 13: the counting gate counts DISCOVERED tests, and a skip discovers.

        Sixteen `setUpClass` methods raise `SkipTest` when a fetched input or a built
        artefact is absent, and `data/raw/*` and the derived tables are gitignored, so
        on a fresh clone whole classes do not run while the suite reports OK. The three
        tests of `EveryCitationAttestsItsAlias`, which are the only gate tying an alias
        to the record it cites, are among them: 35 tests ran where 38 were recorded,
        and `TESTS_PER_FILE` saw 38 either way because discovery is unaffected by a
        skip. That is root cause 15's compounding cause, a red suite reported as green,
        reached by a different route.

        Rather than one assertion per skip site, which would be sixteen typed lists,
        this asserts the condition they all share: in a git checkout every input and
        artefact the library and the gates read is present. The set is derived from the
        manifest, from what the loaders request, and from `facts.missing_inputs`, so a
        new input is covered by whichever of those names it.
        """
        if not (ROOT / ".git").exists():
            raise unittest.SkipTest(
                "not a git checkout, so the artefacts cannot be built here. This is "
                "the sdist or installed case; in a checkout this test is the one that "
                "turns a silent skip into a failure.")
        import check_inputs
        from hsv1coords import _paths
        from hsv1coords.facts import missing_inputs
        data = ROOT / "data"
        absent = [f"data/raw/{row[0]}" for row in check_inputs.manifest_rows()
                  if not (data / "raw" / row[0]).exists()]
        absent += [f"data/raw/{name}" for name in sorted(check_inputs.PINNED_ELSEWHERE)
                   if not (data / "raw" / name).exists()]
        # Exercise every loader first: `requested()` reports what has been asked for in
        # this process, so a loader nothing called contributes nothing.
        import hsv1coords
        hsv1coords.Dictionary.load()
        hsv1coords.ConventionSet.load()
        hsv1coords.RepeatMap.load()
        hsv1coords.VersionCrosswalk.load()
        available = hsv1coords.Crosswalk.available()
        if available:
            hsv1coords.Crosswalk.load(available[0])
        absent += [f"data/{rel}" for rel in _paths.requested()
                   if not data.joinpath(*rel.split("/")).exists()]
        absent += missing_inputs()
        # `HSV1_Annotation.gtf` is in the manifest, so the first loop names it.
        #
        # Round 14: three conditions no file in those three sources names, each the
        # key to a skip the gate could not see. A missing `crosswalk_report.md` let
        # `test_the_alignment_preset_is_not_contradicted` skip with the gate green;
        # minimap2 off PATH let four tests skip, among them the only class that
        # compares committed blocks against a fresh alignment; and an accession with
        # no block file let three more skip.
        if not shutil.which("minimap2"):
            absent.append("minimap2 on PATH (a declared dependency of M3; without it "
                          "the blocks can only be checked against themselves)")
        for name in sorted(p.name for p in (data / "derived").glob("*.md")):
            if not (data / "derived" / name).exists():       # pragma: no cover
                absent.append(f"data/derived/{name}")
        for report in ("crosswalk_report.md", "repeat_map.md", "strain_inventory.md",
                       "version_crosswalk.md", "annotation_summary.md",
                       "mechanical_report.md"):
            if not (data / "derived" / report).exists():
                absent.append(f"data/derived/{report}")
        # `accessions.txt` lists the GI numbers the fetch asked for, one per genome,
        # so its line count is how many block files a complete build leaves behind.
        # Round 14's remaining skip: `tests/test_readme_examples.py` skips when the
        # documents state no Python block at all, which is a condition about the
        # documents rather than about an artefact, and it silences the gate that runs
        # the quickstart a new reader runs first.
        import re as _re
        blocks = sum(len(_re.findall(r"```python\n", (ROOT / name).read_text()))
                     for name in ("README.md", "ROADMAP.md", "FINDINGS.md",
                                  "CHANGELOG.md")
                     if (ROOT / name).exists())
        if not blocks:
            absent.append("any ```python block in the documents; the README examples "
                          "gate skips without one")
        accessions = data / "accessions.txt"
        if accessions.exists():
            wanted = len([line for line in accessions.read_text().splitlines()
                          if line.strip()])
            built = len(list((data / "derived" / "crosswalk").glob("*.blocks.tsv")))
            if built < wanted:
                absent.append(f"{wanted - built} of {wanted} crosswalk block files; "
                              f"the tests keyed on an accession skip without theirs")
        self.assertEqual(
            sorted(set(absent)), [],
            f"{len(set(absent))} inputs or artefacts are absent in a git checkout, so "
            f"the classes that read them are skipping and the suite is smaller than it "
            f"reports. Run `make all`. Missing: {sorted(set(absent))}")

    SKIP_SITES = {
        "test_docs_match_artefacts": 5, "test_m1_dictionary": 2,
        "test_m2_annotation": 11, "test_m3_crosswalk": 7, "test_m4_repeat_map": 3,
        "test_m5_conventions": 9, "test_placeholder": 0, "test_readme_examples": 2,
    }
    r"""Every way each test file can decline to run something.

    Round 14: a new skip is a new way for the suite to report OK over less than it
    says, so adding one fails here until its condition is added to the checkout gate
    below. Round 15: the count was `re.findall(r"SkipTest\(")`, which sees one of the
    three spellings. It missed all 11 `self.skipTest(...)` calls and every
    `@unittest.skipUnless`/`skipIf` decorator, so putting `SequenceAgreement` behind
    an environment variable left the suite green at `OK (skipped=10)`, and it counted
    a `SkipTest` mentioned in a docstring, so a real skip could be added for free by
    rewording that sentence. It is counted from the syntax tree now: 29 sites against
    the 18 the regex saw.
    """

    @staticmethod
    def _names_of(node) -> str:
        return (node.attr if isinstance(node, ast.Attribute)
                else getattr(node, "id", ""))

    @classmethod
    def _skip_sites(cls, path) -> int:
        """Skips of every spelling this file can hold.

        Round 16: `ast.Call` matched `raise unittest.SkipTest(...)` and
        `self.skipTest(...)` and missed `raise _ut.SkipTest` with no parentheses,
        which is legal and raises the class. One of those as the first statement of
        `SequenceAgreement.setUpClass` removed ten tests, among them the per-genome
        base-agreement equality gate that closed cause 27, with the suite reporting
        `OK (skipped=1)` and this count unchanged.
        """
        sites = 0
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call):
                if cls._names_of(node.func) in ("SkipTest", "skipTest"):
                    sites += 1
            elif isinstance(node, ast.Raise) and node.exc is not None:
                # A bare `raise SkipTest` raises the class itself.
                if cls._names_of(node.exc) == "SkipTest":
                    sites += 1
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                for decorator in node.decorator_list:
                    target = (decorator.func if isinstance(decorator, ast.Call)
                              else decorator)
                    if cls._names_of(target).startswith("skip"):
                        sites += 1
        return sites

    @classmethod
    def _expected_failures(cls, path) -> int:
        """`@unittest.expectedFailure` sites, which are a skip under another name.

        Round 16: the name does not begin with "skip", the test stays discovered so
        `TESTS_PER_FILE` is unmoved, and the run reports OK. Decorating the GTF
        multiset gate and then relabelling twelve rows reported
        `OK (expected failures=1)`.
        """
        sites = 0
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                 ast.ClassDef)):
                for decorator in node.decorator_list:
                    target = (decorator.func if isinstance(decorator, ast.Call)
                              else decorator)
                    if cls._names_of(target) == "expectedFailure":
                        sites += 1
        return sites

    def test_no_test_is_marked_as_an_expected_failure(self):
        """There are none, and one would be a passing test that does not pass."""
        marked = {path.stem: self._expected_failures(path)
                  for path in sorted((ROOT / "tests").glob("test_*.py"))}
        self.assertEqual(
            {k: v for k, v in marked.items() if v}, {},
            f"a test is decorated `@unittest.expectedFailure`, which keeps it "
            f"discovered, keeps the suite green and asserts nothing: {marked}")

    def test_no_test_file_gains_a_skip_the_checkout_gate_does_not_cover(self):
        counted = {path.stem: self._skip_sites(path)
                   for path in sorted((ROOT / "tests").glob("test_*.py"))}
        self.assertEqual(
            counted, self.SKIP_SITES,
            "the number of skip sites has changed. Every skip is a test that does "
            "not run, so each one's condition belongs in "
            "test_a_checkout_holds_every_artefact_the_gates_skip_without, which fails "
            "in a checkout rather than skipping. Update that gate, then this count.")

    def test_no_test_class_defines_the_same_method_twice(self):
        """Round 15, twice in one session, both times mine.

        A slice rewrite of `tests/test_m4_repeat_map.py` deleted a whole test and I
        hid it by lowering `TESTS_PER_FILE`, which that table's own docstring says not
        to do; restoring it from git then duplicated a neighbouring test, and the
        second definition shadowed the first silently, so the file reported 58 tests
        while running 57 distinct ones. Python allows it, discovery counts one, and
        the diff is the only place it shows: root cause 40's shape in a different
        edit.
        """
        import ast
        duplicates = []
        for path in sorted((ROOT / "tests").glob("test_*.py")):
            for node in ast.parse(path.read_text()).body:
                if not isinstance(node, ast.ClassDef):
                    continue
                names = [n.name for n in node.body
                         if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
                for name in sorted(set(names)):
                    if names.count(name) > 1:
                        duplicates.append(f"{path.name}:{node.name}.{name}")
        self.assertEqual(
            duplicates, [],
            f"{len(duplicates)} methods are defined twice in one class, so the "
            f"second shadows the first and only one of them runs: {duplicates}")

    def test_every_test_file_parses_without_a_syntax_warning(self):
        r"""Round 15: a docstring in test_m4_repeat_map.py held `[\d,]*` unescaped.

        Python warns and will one day refuse. Found while counting skip sites with
        `ast`, which is the sort of thing a regex over source text cannot see.
        """
        import ast
        import warnings
        offenders = []
        for path in sorted((ROOT / "tests").glob("test_*.py")):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                ast.parse(path.read_text())
            for warning in caught:
                offenders.append((path.name, str(warning.message)[:60]))
        self.assertEqual(offenders, [], f"{len(offenders)} test files parse with a "
                                        f"warning: {offenders}")

    def test_no_file_contributes_fewer_tests_than_it_did(self):
        import unittest as _unittest
        counted: dict[str, int] = {}

        def walk(suite):
            for item in suite:
                if isinstance(item, _unittest.TestSuite):
                    walk(item)
                else:
                    name = type(item).__module__
                    counted[name] = counted.get(name, 0) + 1

        walk(_unittest.TestLoader().discover(str(ROOT / "tests")))
        lost, gained = [], []
        for name, was in sorted(self.TESTS_PER_FILE.items()):
            now = counted.get(name, 0)
            if now < was:
                lost.append(f"{name}: {now} tests discovered, {was} recorded")
            elif now > was:
                gained.append(f"{name}: {was} -> {now}")
        self.assertEqual(
            lost, [],
            f"{len(lost)} test files run fewer tests than TESTS_PER_FILE records. A "
            f"class or a method has stopped being discovered, which is what round 11's "
            f"slice edit did to a whole class with the suite green:\n  "
            + "\n  ".join(lost))
        # Round 16: a rise was printed, and round 15's own last commit shipped the
        # pin one behind, so the single test holding its cause-68 fix together could
        # be deleted with the suite green and the notice silenced. Raising the entry
        # is one line in the commit that added the test (cause 61), and the slack is
        # worth more than the line.
        self.assertEqual(
            gained, [],
            f"{len(gained)} test files run MORE tests than TESTS_PER_FILE records, so "
            f"the pin has slack and a new test could be deleted inside it. Raise each "
            f"entry in the commit that added the tests:\n  " + "\n  ".join(gained))

class VersionNamesAreValidatedBeforeAnswering(unittest.TestCase):
    """Round 7: `convert` short-circuited on source == target before validating.

    `convert(100, "hg38", "hg38")` returned 100, so an entirely unknown reference was
    accepted as long as both names matched. Two further mutations survived here: the
    identity path returning None, and the unknown-pair branch returning None instead of
    raising, each changing 314 live answers.
    """

    @classmethod
    def setUpClass(cls):
        from hsv1coords.versions import VersionCrosswalk
        cls.vc = VersionCrosswalk.load()

    def test_an_unknown_reference_is_refused_even_when_both_names_match(self):
        with self.assertRaises(ValueError):
            self.vc.convert(100, "hg38", "hg38")
        with self.assertRaises(ValueError):
            self.vc.convert(100, "NC_001806.2", "GRCh38")

    def test_the_identity_case_still_returns_the_position(self):
        from hsv1coords.versions import NEW, OLD
        self.assertEqual(self.vc.convert(100, NEW, NEW), 100)
        self.assertEqual(self.vc.convert(100, OLD, OLD), 100)

    def test_the_two_known_directions_still_work(self):
        from hsv1coords.versions import NEW, OLD
        self.assertEqual(self.vc.convert(152261, OLD, NEW), 152222)
        self.assertEqual(self.vc.convert(152222, NEW, OLD), 152261)

class ASplitReferenceHasARouteBackOut(unittest.TestCase):
    """Round 8: a user holding `UL:50788` from Morse's own VCFs had nowhere to go.

    `locate` goes reference to pair, `coordinate_in` goes reference to pair, and
    `convert` refused and then named `coordinate_in`, which is the same direction
    again. `to_reference` did answer, with an integer that appears in no Morse file.
    """

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()
        cls.m = cls.cs["morse2017_nr"]

    def test_from_record_is_the_inverse_of_locate(self):
        # Native throughout: `locate` takes a position in the convention's own
        # version and these are NC_001806.1 positions, which round 16 made explicit.
        for ref in (1, 5000, 9212, 9213, 60000, 117159, 125974, 132604, 132605,
                    145584):
            hit = self.m.locate(ref, native=True)
            if hit is None:
                continue
            with self.subTest(reference=ref, pair=hit):
                self.assertEqual(self.m.from_record(*hit, native=True), ref)

    def test_every_record_boundary_round_trips(self):
        """In the convention's own version, which is what `locate` speaks.

        Round 16: `from_record` answered in the convention's own NC_001806.1 frame
        under a docstring promising "back to the reference", so it now translates to
        the base version by default and `native=True` is the untranslated route. The
        round trip against `locate` is a native-to-native question, so it asks for
        the native answer explicitly.
        """
        for name, a, b in self.m.records:
            for offset, expected in ((1, a), (b - a + 1, b)):
                with self.subTest(record=name, offset=offset):
                    self.assertEqual(self.m.from_record(name, offset, native=True),
                                     expected)
                    self.assertEqual(self.m.locate(expected, native=True),
                                     (name, offset))
            with self.subTest(record=name, offset="past the end"):
                self.assertIsNone(self.m.from_record(name, b - a + 2, native=True))

    def test_from_record_answers_in_the_base_version_by_default(self):
        """The frame travels with the value, and the two routes agree.

        A user holding `UL:50788` from Morse's deposited VCFs asks `from_record`, and
        the number they get must be the one `coordinate_in` would turn back into
        `UL:50788`. Round 16 found 109,639 of 152,222 positions where the
        `Convention` methods and the set disagreed, because the version translation
        lived only in the set.
        """
        from hsv1coords.versions import VersionCrosswalk
        for record, offset in (("UL", 50788), ("RL", 1), ("US", 12391)):
            with self.subTest(record=record, offset=offset):
                base = self.m.from_record(record, offset)
                own = self.m.from_record(record, offset, native=True)
                self.assertEqual(
                    base,
                    VersionCrosswalk.load().convert(own, "NC_001806.1",
                                                    "NC_001806.2"),
                    "from_record's answer is not the base-version coordinate")
                self.assertEqual(self.cs.coordinate_in(base, "morse2017_nr"),
                                 f"{record}:{offset}",
                                 "the two routes disagree about the same position")

    def test_a_split_reference_refuses_a_bare_integer_in_both_directions(self):
        with self.assertRaises(ValueError):
            self.m.to_convention(50788)
        with self.assertRaises(ValueError):
            self.m.to_reference(50788)

    def test_an_unknown_record_or_offset_is_refused(self):
        with self.assertRaises(KeyError):
            self.m.from_record("ZZ", 1)
        for bad in (0, -1, 1.5, "1", True):
            with self.subTest(offset=bad):
                with self.assertRaises(ValueError):
                    self.m.from_record("UL", bad)

    def test_a_non_split_convention_refuses_from_record(self):
        with self.assertRaises(ValueError):
            self.cs["szpara2014_trimmed"].from_record("UL", 1)


class DocstringFiguresMatchTheArtefact(unittest.TestCase):
    """The documentation gate reads four top-level documents, so `src/` is uncovered.

    Round 8 found every sign in `ReferenceVersionMismatch`'s docstring inverted
    against `version_crosswalk.tsv`, in the text a user reads to decide whether to
    trust a conversion. This gates the figures in that docstring specifically, which is
    narrower than gating every docstring and is the one that describes an artefact.
    """

    def test_the_version_mismatch_docstring_states_the_real_offsets(self):
        import re

        from hsv1coords.conventions import ReferenceVersionMismatch
        from hsv1coords.facts import SPELLED
        from hsv1coords.versions import VersionCrosswalk
        if not (DERIVED / "version_crosswalk.tsv").exists():
            self.skipTest("version_crosswalk.tsv is built; run `make versions`")
        vc = VersionCrosswalk.load()
        # Flatten newlines too: "nine distinct\n    offsets" matches nothing
        # otherwise, which is the line-break trap in standing rule 16, hit again here.
        doc = re.sub(r"\s+", " ", ReferenceVersionMismatch.__doc__ or "")

        stated = re.findall(r"offset at ([\d,]+) is ([+-]\d+)", doc)
        self.assertGreaterEqual(
            len(stated), 4,
            "the docstring states fewer than four measured offsets, so it cannot be "
            "checked positionally")
        for position, offset in stated:
            pos = int(position.replace(",", ""))
            with self.subTest(position=pos):
                self.assertEqual(
                    vc.offset_at(pos), int(offset),
                    f"the docstring says the offset at {position} is {offset}; the "
                    f"crosswalk says {vc.offset_at(pos)}")

        count = len(vc.distinct_offsets)
        word = SPELLED.get(count, str(count))
        self.assertTrue(f"{count} distinct offsets" in doc
                        or f"{word} distinct offsets" in doc,
                        f"the docstring does not state the {count} distinct offsets")
        self.assertIn(f"{len(vc.blocks)} aligned blocks", doc)

class TheGenomeLengthsAreNotJustTypedIn(unittest.TestCase):
    """Round 8: `GENOME_BP` is declared in five files and gated nowhere.

    Setting `repeats.GENOME_BP` to 152200 left the suite green at 173 tests, while
    `facts.py` already derives the same figure from the record. A wrong value makes
    `classify` declare the tail of the genome outside the reference, so `lift` silently
    declines there: a confident refusal, which is the quiet half of the error this
    repository exists to prevent.

    The constants stay, because a module that reads a GenBank file to learn a constant
    is worse, but they are checked against the records here.
    """

    @classmethod
    def setUpClass(cls):
        cls.raw = ROOT / "data" / "raw"
        for name in ("ref_NC_001806.2.gb", "ref_NC_001806.1.gb"):
            if not (cls.raw / name).exists():
                raise unittest.SkipTest(f"{name} is fetched, not committed; "
                                        f"run `make inputs` first")
        sys.path.insert(0, str(ROOT / "src"))
        from build_crosswalk import genbank_sequence
        cls.lengths = {
            "NC_001806.2": len(genbank_sequence(cls.raw / "ref_NC_001806.2.gb")),
            "NC_001806.1": len(genbank_sequence(cls.raw / "ref_NC_001806.1.gb")),
        }

    def test_every_declared_genome_length_matches_the_record(self):
        from hsv1coords import conventions, repeats
        import build_repeat_map
        import check_inputs
        expected = self.lengths["NC_001806.2"]
        for module in (repeats, conventions, build_repeat_map, check_inputs):
            with self.subTest(module=module.__name__):
                self.assertEqual(
                    module.GENOME_BP, expected,
                    f"{module.__name__}.GENOME_BP is {module.GENOME_BP} and "
                    f"NC_001806.2 is {expected} bases")

    def test_the_version_lengths_match_both_records(self):
        from hsv1coords.conventions import VERSION_LENGTHS
        self.assertEqual(VERSION_LENGTHS, self.lengths,
                         f"VERSION_LENGTHS is {VERSION_LENGTHS} and the records are "
                         f"{self.lengths}")

    def test_the_declared_lengths_agree_with_each_other(self):
        """Five declarations of one number, so they are checked against each other too.

        A reader changing one would otherwise see a green suite on four of five files.
        """
        from hsv1coords import conventions, repeats
        import build_repeat_map
        import check_inputs
        values = {m.__name__: m.GENOME_BP for m in
                  (repeats, conventions, build_repeat_map, check_inputs)}
        self.assertEqual(len(set(values.values())), 1,
                         f"the declared genome lengths disagree: {values}")


class TheInputGateIsActuallyRun(unittest.TestCase):
    """Round 12: nothing in the suite executed `src/check_inputs.py`.

    It ran only via `make check` in CI, and the suite imported it for one constant. So
    every assumption in the M0 gate could be made inert without a test noticing:
    changing `ref_cds == 77` to `ref_cds == ref_cds` left the suite green AND left
    `make check` printing "ok   reference carries 77 CDS" and "all inputs and
    assumptions verified", exit 0. The gate that exists to catch an input moving
    reported success while asserting nothing, and the same edit worked on any of the
    eight assumptions.

    The gate is now built as two tables, a measurement per assumption and a literal per
    assumption, with one comparison between them, and this runs it.
    """

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "src"))
        for name in ("ref_NC_001806.2.gb", "alt_BK012101.1.gb", "base_JN555585.1.gb"):
            if not (ROOT / "data" / "raw" / name).exists():
                raise unittest.SkipTest(f"{name} is fetched, not committed; "
                                        f"run `make inputs` first")

    def test_the_gate_passes_on_the_committed_inputs(self):
        import check_inputs
        failures = [a for a in check_inputs.assumptions() if not a.ok]
        self.assertEqual(
            [(a.key, a.got, a.want) for a in failures], [],
            f"{len(failures)} of the M0 assumptions do not hold on the inputs in "
            f"data/raw")

    def test_every_expectation_is_a_literal_and_not_a_measurement(self):
        """The split is the fix, so the split itself is asserted.

        An expectation that is anything but a plain number, string or boolean could be
        a function of the measurement, which is how a tautology gets written.
        """
        import check_inputs
        for key, (want, label) in check_inputs.EXPECTED.items():
            with self.subTest(key=key):
                self.assertIsInstance(want, (int, str, bool),
                                      f"{key}'s expectation is a {type(want).__name__}")
                self.assertTrue(label.strip(), f"{key} states no sentence")
        self.assertEqual(set(check_inputs.EXPECTED), set(check_inputs.MEASUREMENTS))
        # And the other direction: a measurement that reads the expectation table would
        # pass every comparison while reading no record at all.
        import inspect
        for key, measure in check_inputs.MEASUREMENTS.items():
            with self.subTest(key=key):
                source = inspect.getsource(measure)
                self.assertNotIn(
                    "EXPECTED", source,
                    f"{key}'s measurement reads the expectation table, so the gate "
                    f"compares the answer with itself")
        self.assertGreaterEqual(
            len(check_inputs.EXPECTED), 8,
            "the M0 gate asserts fewer assumptions than it did when this was written")

    def test_every_measurement_is_sensitive(self):
        """Each measurement must disagree with a wrong expectation.

        A measurement that cannot disagree with anything is the tautology this class
        exists for. Perturbing the expectation rather than the record keeps this cheap
        enough to run every time.

        Round 13: this leg perturbs the expectation only, so it says nothing about the
        measurement. `test_every_measurement_actually_reads_its_record` below is the
        leg that does, and it is the one that catches a literal.
        """
        import check_inputs
        inert = []
        for a in check_inputs.assumptions():
            wrong = (not a.want) if isinstance(a.want, bool) else (
                a.want + 1 if isinstance(a.want, int) else str(a.want) + "x")
            if check_inputs.Assumption(a.key, a.label, a.got, wrong).ok:
                inert.append((a.key, a.got, wrong))
        self.assertEqual(
            inert, [],
            f"{len(inert)} assumptions hold against a deliberately wrong expectation, "
            f"so they are not comparing the record to anything: {inert}")

    def test_every_measurement_actually_reads_its_record(self):
        """Round 13: round 12's four new tests all perturbed the expectation.

        So the shape they were written to prevent came straight back in the other
        table. Replacing a measurement with the literal it is compared against,
        `"ref_cds": lambda: 77`, passed all 286 tests, and `make check` went on
        printing "ok reference carries 77 CDS" over a reference that no longer carried
        77: the sensitivity test perturbed `want`, and `got` was a constant.

        A measurement is sensitive if its own input moves it. Both legs here point
        `check_inputs.RAW` somewhere else, so neither needs a per-measurement mutation
        table, which is the shape standing rule 13 warns about: with the records absent
        the measurement must fail rather than answer, and with the records replaced by
        stubs every measured value must move. A literal survives neither.
        """
        import check_inputs
        committed = {k: m() for k, m in check_inputs.MEASUREMENTS.items()}
        real = check_inputs.RAW

        def measured_with(raw):
            check_inputs.RAW = raw
            out = {}
            for key, measure in check_inputs.MEASUREMENTS.items():
                try:
                    out[key] = measure()
                except Exception as exc:          # noqa: BLE001 - that is the signal
                    out[key] = f"raised {type(exc).__name__}"
            return out

        try:
            with tempfile.TemporaryDirectory() as tmp:
                absent = measured_with(pathlib.Path(tmp))
            answered = [(k, v) for k, v in absent.items()
                        if not str(v).startswith("raised ")]
            self.assertEqual(
                answered, [],
                f"{len(answered)} M0 measurements answered with their records absent, "
                f"so they are not reading them: {answered}")

            # Round 14: stubbing every record at once asks only whether the value
            # moved, never which file produced it. `base_JN555585.1.gb` is the
            # identical sequence to the reference, so swapping `_ref()` for `_base()`
            # in four measurements left `make check` certifying four assumptions about
            # NC_001806.2 against a different accession, with the whole suite green.
            # Each record is stubbed on its own now, and `check_inputs.READS`, which
            # declares what each measurement reads, is checked against what moves
            # rather than trusted.
            # Round 15: the loop was driven by `READS.values()`, so only the three
            # declared files were ever stubbed and a dependency on a FOURTH input was
            # invisible. A defensive one passed: `"ref_cds": lambda: count_feature(
            # _ref(), "CDS") if not (RAW / "HSV1_Annotation.gtf").exists() or "UL19"
            # in (RAW / "HSV1_Annotation.gtf").read_text() else 0` left 307 tests
            # green with `make check` printing "ok the reference carries 77 CDS",
            # and an upstream change to the GTF deciding a reference assumption. The
            # loop is driven by what `data/raw` HOLDS now, so a measurement that
            # reads anything it does not declare is caught whatever it reads.
            present = sorted(path.name for path in real.iterdir()
                             if path.is_file() and path.name not in (".gitkeep",))
            declared = sorted({name for names in check_inputs.READS.values()
                               for name in names})
            self.assertTrue(declared, "no measurement declares a record")
            for name in declared:
                self.assertIn(name, present,
                              f"{name} is declared and is not in data/raw")
            stubs_for = {"ref_NC_001806.2.gb": "acgt",
                         "alt_BK012101.1.gb": "ggca",
                         "base_JN555585.1.gb": "tttt"}
            records = present
            wrong = []
            for target in records:
                with tempfile.TemporaryDirectory() as tmp:
                    stubs = pathlib.Path(tmp)
                    for name in records:
                        if name != target:
                            (stubs / name).write_bytes((real / name).read_bytes())
                        elif name in stubs_for:
                            (stubs / name).write_text(
                                f"LOCUS       STUB   4 bp    DNA     linear   VRL\n"
                                f"FEATURES             Location/Qualifiers\n"
                                f"ORIGIN\n        1 {stubs_for[name]}\n//\n")
                        else:
                            # An input no measurement claims to read: replaced with a
                            # file that is not what it was, so a measurement that
                            # reads it undeclared moves and is named below.
                            (stubs / name).write_text("stub\n")
                    stubbed = measured_with(stubs)
                for key, declared in check_inputs.READS.items():
                    moved = stubbed[key] != committed[key]
                    if target in declared and not moved:
                        wrong.append((key, f"declares {target} and did not move when "
                                           f"it was replaced by a stub"))
                    elif target not in declared and moved:
                        wrong.append((key, f"moved when {target} was stubbed and does "
                                           f"not declare reading it"))
            self.assertEqual(
                wrong, [],
                f"{len(wrong)} M0 measurements disagree with check_inputs.READS about "
                f"which record they read, so a measurement is either a constant or is "
                f"reading a different accession from the one its sentence names: "
                f"{wrong}")
        finally:
            check_inputs.RAW = real
        # And the committed values are unchanged by the detour.
        self.assertEqual({k: m() for k, m in check_inputs.MEASUREMENTS.items()},
                         committed)

    def test_a_file_in_data_raw_that_nothing_pins_is_named(self):
        """Round 12: the manifest loop iterated the manifest's own rows.

        So a file sitting in `data/raw` with no row was never checksummed and never
        missed, which is how `genomes.fasta` came to be outside this gate entirely
        until round 11 pinned it elsewhere. The other direction is checked now, and
        what pins each exception is stated.
        """
        import check_inputs
        self.assertEqual(check_inputs.check_manifest(), [],
                         "the M0 manifest check does not pass on the committed inputs")
        self.assertIn("genomes.fasta", check_inputs.PINNED_ELSEWHERE,
                      "genomes.fasta is in data/raw and in no manifest row, so it must "
                      "state what pins it")
        for name, pin in check_inputs.PINNED_ELSEWHERE.items():
            with self.subTest(name=name):
                self.assertTrue((ROOT / pin).exists(),
                                f"{name} says it is pinned by {pin}, which is missing")


class TheLibraryRefusesWhatItCannotAnswer(unittest.TestCase):
    """Round 16, from a review that used the library as a stranger would.

    Each of these was a wrong answer or an undiagnosable crash at a public entry
    point, and each is now a refusal that names the fix.
    """

    def test_an_empty_or_truncated_repeat_map_is_refused(self):
        """A header-only artefact made every position well defined.

        `repeats.py`'s structural check was `if rows and ...`, so an empty table
        skipped it: 0 regions, 0 pairs, no error, `classify(118777)` well defined,
        `lift(118777, ref, 'JQ780693.1')` answering 118,185 where the README
        documents None, and the report's ambiguity caveat gone. Reachable without
        tampering, because the builder writes the header and then loops.
        """
        from hsv1coords import RepeatMap
        with tempfile.TemporaryDirectory() as tmp:
            empty = pathlib.Path(tmp)
            for name in ("repeat_map.tsv", "diploid_pairs.tsv"):
                source = (ROOT / "data" / "derived" / name).read_text()
                (empty / name).write_text(source.splitlines()[0] + "\n")
            with self.assertRaises(ValueError) as refusal:
                RepeatMap.load(empty)
            self.assertIn("make repeatmap", str(refusal.exception))
        # And the pair table is required, not optional: deleting it emptied every
        # locus list in silence.
        with tempfile.TemporaryDirectory() as tmp:
            partial = pathlib.Path(tmp)
            (partial / "repeat_map.tsv").write_text(
                (ROOT / "data" / "derived" / "repeat_map.tsv").read_text())
            with self.assertRaises(FileNotFoundError):
                RepeatMap.load(partial)

    def test_an_empty_convention_table_is_refused(self):
        """It made `comparable_everywhere` True over zero conventions."""
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "conventions.tsv"
            source = (ROOT / "data" / "conventions.tsv").read_text()
            path.write_text(source.splitlines()[0] + "\n")
            with self.assertRaises(ValueError):
                ConventionSet.load(path)

    def test_a_coordinate_must_be_an_int_at_every_entry_point(self):
        """Round 16: five public methods leaked a raw comparison TypeError.

        `RepeatMap.counterparts('900')` gave "'<=' not supported between instances of
        'int' and 'str'", which names nothing a caller can act on, and
        `lift(None, REFERENCE, target)` answered None and printed a decline about a
        position that does not exist.
        """
        from hsv1coords import Crosswalk, Dictionary, RepeatMap, lift
        from hsv1coords.liftover import lift_detailed
        m, cs = RepeatMap.load(), ConventionSet.load()
        for label, call in (
                ("counterparts", lambda: m.counterparts("900")),
                ("counterparts_detailed", lambda: m.counterparts_detailed("900")),
                ("loci_at", lambda: m.loci_at("900")),
                ("convert", lambda: cs.convert("60000", "full_length",
                                               "szpara2014_trimmed")),
                ("coordinate_in", lambda: cs.coordinate_in("60000", "full_length")),
                ("to_reference", lambda: cs["szpara2014_trimmed"].to_reference("1")),
                ("locate", lambda: cs["morse2017_nr"].locate("60000")),
                ("lift(None) from the reference",
                 lambda: lift(None, "NC_001806.2", "JQ780693.1")),
                ("lift_detailed(None) from the reference",
                 lambda: lift_detailed(None, "NC_001806.2", "JQ780693.1")),
                ("resolve(None)", lambda: Dictionary.load().resolve(None))):
            with self.subTest(entry_point=label):
                with self.assertRaises(TypeError) as refusal:
                    call()
                self.assertRegex(str(refusal.exception),
                                 r"must be an int|must be a string")

    def test_a_convention_in_another_version_refuses_a_base_position(self):
        """109,639 of 152,222 positions disagreed with the set, four of them naming a
        different record."""
        from hsv1coords.conventions import ReferenceVersionMismatch
        cs = ConventionSet.load()
        morse = cs["morse2017_nr"]
        for call in (lambda: morse.coordinate(60000), lambda: morse.locate(60000)):
            with self.assertRaises(ReferenceVersionMismatch) as refusal:
                call()
            self.assertIn("coordinate_in", str(refusal.exception))
        self.assertEqual(morse.coordinate(60000, native=True), "UL:50788")
        self.assertEqual(cs.coordinate_in(60000, "morse2017_nr"), "UL:50787")

    def test_a_position_outside_a_convention_has_no_reference_coordinate(self):
        """`to_reference(0)` gave 9,213, inside the interval the convention removes."""
        cs = ConventionSet.load()
        trimmed = cs["szpara2014_trimmed"]
        for bad in (0, -20, trimmed.length() + 1):
            with self.subTest(position=bad):
                self.assertIsNone(trimmed.to_reference(bad))
        self.assertIsNone(cs["full_length"].to_convention(1000000))

    def test_convert_answers_the_same_way_in_both_directions(self):
        """The exclusion check was applied to the target only."""
        cs = ConventionSet.load()
        excluded = 118777
        self.assertIsNone(cs.convert(excluded, "full_length",
                                     "lassalle2020_excluded"))
        self.assertIsNone(cs.convert(excluded, "lassalle2020_excluded",
                                     "full_length"))
        retained = 60000
        self.assertEqual(cs.convert(retained, "full_length",
                                    "lassalle2020_excluded"), retained)
        self.assertEqual(cs.convert(retained, "lassalle2020_excluded",
                                    "full_length"), retained)

    def test_the_interval_bases_table_cannot_be_rewritten_by_a_caller(self):
        """One assignment turned the inference caveat off for every later caller."""
        from hsv1coords.conventions import Convention
        with self.assertRaises(TypeError):
            Convention.INTERVAL_BASES["inferred_from_record"] = (True, "nothing")

    def test_the_cached_loaders_hand_one_object_to_every_caller(self):
        """Round 16: `load(x)` and `load(x, None)` were two cache entries.

        The docstrings' "hands the same object to every caller" is the premise of the
        immutability argument in three modules, and it was not true as written: a
        mixed workload held two objects per accession against a cache sized for one.
        """
        from hsv1coords import Crosswalk, Dictionary, RepeatMap
        accession = Crosswalk.available()[0]
        self.assertIs(Crosswalk.load(accession), Crosswalk.load(accession, None))
        self.assertIs(RepeatMap.load(), RepeatMap.load(None))
        self.assertIs(Dictionary.load(), Dictionary.load(None))

    def test_an_accession_is_matched_against_the_built_set(self):
        """Round 16: the lookup was a filesystem probe, so it was platform-dependent.

        On a case-insensitive filesystem `Crosswalk.load('jq673480.1')` answered from
        `JQ673480.1.blocks.tsv` and returned a crosswalk whose accession
        `available()` does not list; on Linux the same call raised.
        """
        from hsv1coords import Crosswalk
        with self.assertRaises(FileNotFoundError) as refusal:
            Crosswalk.load(Crosswalk.available()[0].lower())
        self.assertIn("available()", str(refusal.exception))

    def test_an_impossible_position_has_no_unmapped_stretch(self):
        """None there means "that position lifts fine", which it does not."""
        from hsv1coords import Crosswalk
        cw = Crosswalk.load(Crosswalk.available()[0])
        for bad in (0, -1, (cw.deposited_length or 0) + 5000):
            with self.subTest(position=bad):
                with self.assertRaises(ValueError):
                    cw.unmapped_stretch(bad)


class TheCachedLoadersCannotBeCorruptedByACaller(unittest.TestCase):
    """Round 9: `load` is cached, so every caller holds the same object.

    Nothing in the repository tested the caching at all. `Crosswalk.blocks` and
    `RepeatMap.regions` were plain lists, so one caller clearing or sorting one changed
    what every later `load` returned: `to_reference(60000)` went from 61474 to None for
    a caller that did nothing wrong. The docstrings said the objects were "treated as
    immutable within a run", which is a request rather than a guarantee, and this is the
    guarantee.
    """

    def test_a_cached_crosswalk_cannot_be_mutated_through_its_blocks(self):
        from hsv1coords import Crosswalk
        cw = Crosswalk.load("JQ673480.1", CROSSWALK_DIR)
        before = cw.to_reference(60000)
        self.assertIsNotNone(before)
        with self.assertRaises(AttributeError):
            cw.blocks.clear()
        self.assertEqual(Crosswalk.load("JQ673480.1", CROSSWALK_DIR)
                         .to_reference(60000), before)

    def test_a_cached_repeat_map_cannot_be_mutated_through_its_regions(self):
        from hsv1coords import RepeatMap
        m = RepeatMap.load()
        before = len(m.regions)
        with self.assertRaises(AttributeError):
            m.regions.clear()
        with self.assertRaises(AttributeError):
            m.pairs.clear()
        self.assertEqual(len(RepeatMap.load().regions), before)

    def test_a_convention_cannot_be_mutated_through_its_intervals(self):
        """Round 12: the two siblings round 9's fix missed.

        `Convention` is `frozen=True`, which stops a field being reassigned and says
        nothing about what a field points at. `intervals` was a list, so
        `cs["szpara2014_trimmed"].intervals.append((50000, 70000))` changed that
        convention's `to_convention(60000)` from 50787 to None for every later holder
        of that set.

        Round 13: the second assertion used to re-load the set and compare, and
        `ConventionSet.load` is not cached, unlike `Crosswalk`, `RepeatMap` and
        `VersionCrosswalk`. So it passed because the file was read again and asserted
        nothing about the object that was attacked. It asserts the same object now,
        which is where the guarantee has to hold.
        """
        cs = ConventionSet.load()
        before = cs["szpara2014_trimmed"].to_convention(60000)
        self.assertIsNotNone(before)
        with self.assertRaises(AttributeError):
            cs["szpara2014_trimmed"].intervals.append((50000, 70000))
        self.assertEqual(cs["szpara2014_trimmed"].to_convention(60000), before)

    def test_a_convention_cannot_be_mutated_through_its_records(self):
        """Round 13: the third mutable field in the dataclass round 12 hardened.

        `records` carries Morse's four record extents, and `usable` reads them to
        place a position in a named record, so `cs["morse2017_nr"].records.clear()`
        turned `usable(140000, "morse2017_nr")` from "present as US:7,392" into
        "inside a removed block, which is not in the NR-reference" for every later
        holder of that set. Observed on HEAD before the fix.
        """
        cs = ConventionSet.load()
        before = cs.usable(140000, "morse2017_nr")
        self.assertTrue(before[0], f"the fixture position is unusable: {before}")
        with self.assertRaises(AttributeError):
            cs["morse2017_nr"].records.clear()
        self.assertEqual(cs.usable(140000, "morse2017_nr"), before)

    def test_a_set_cannot_be_mutated_through_its_list_of_conventions(self):
        """Round 13: the fourth container holding the same data, and the widest.

        `ConventionSet.conventions` was a plain list, and it is what `facts.py`,
        `report` and every gate in this suite walk, so a caller sorting or clearing it
        changed what the next caller saw. A caller constructing a set from a list still
        works: the constructor copies it into a tuple, so a later change to that list
        cannot reach the set either.
        """
        from hsv1coords.conventions import Convention
        rows = list(ConventionSet.load().conventions)
        cs = ConventionSet(rows)
        before = len(cs.conventions)
        with self.assertRaises(AttributeError):
            cs.conventions.append(rows[0])
        rows.clear()
        self.assertEqual(len(cs.conventions), before,
                         "the set followed the list it was constructed from")
        self.assertIsInstance(cs.conventions, tuple)
        # And a row constructed with lists, as the synthetic fixtures in this file do,
        # is frozen just the same: round 13 found the guarantee belonged to the loader
        # rather than to the type.
        row = Convention(name="synthetic_lists", kind="exclude",
                         intervals=[[1, 100]], reference="NC_001806.2", verified=True,
                         source="Synthetic row, round 13, for the coercion only",
                         quote="synthetic", records=[["RL", 1, 100]])
        self.assertEqual(row.intervals, ((1, 100),))
        self.assertEqual(row.records, (("RL", 1, 100),))
        with self.assertRaises(AttributeError):
            row.intervals.append((2, 3))
        with self.assertRaises(AttributeError):
            row.records.append(("US", 2, 3))

    def test_the_version_crosswalk_cannot_be_mutated_through_its_blocks(self):
        """`blocks.clear()` made `to_new` raise IndexError rather than answer None."""
        from hsv1coords.versions import VersionCrosswalk
        vc = VersionCrosswalk.load()
        before = vc.to_new(1000)
        with self.assertRaises(AttributeError):
            vc.blocks.clear()
        self.assertEqual(VersionCrosswalk.load().to_new(1000), before)

    def test_the_version_crosswalk_is_read_once_per_process(self):
        """Round 12: it was re-read and re-parsed once per Morse query.

        `_to_convention_version` called `VersionCrosswalk.load()` on every invocation,
        so a sweep over positions under the Morse convention re-parsed the TSV for each
        one. `Crosswalk` and `RepeatMap` are cached for exactly this reason.
        """
        cs = ConventionSet.load()
        cs.usable(118777, "morse2017_nr")
        info = type(cs)._version_crosswalk.cache_info()
        before = info.misses
        for _ in range(50):
            cs.usable(118777, "morse2017_nr")
        after = type(cs)._version_crosswalk.cache_info().misses
        self.assertEqual(after, before,
                         f"50 further queries re-read the version crosswalk "
                         f"{after - before} times")

    def test_the_crosswalk_cache_holds_the_whole_corpus(self):
        """At `maxsize=64` a sweep over 183 genomes gave 6 hits and 185 misses."""
        from hsv1coords import Crosswalk
        from hsv1coords.liftover import DEFAULT_DIR
        available = Crosswalk.available(CROSSWALK_DIR)
        # Round 16 split `load` into a normalising wrapper and the cached `_load`,
        # so the cache lives on the latter.
        self.assertGreaterEqual(
            Crosswalk.__dict__["_load"].__func__.cache_parameters()["maxsize"],
            len(available),
            f"the cache holds fewer entries than the {len(available)} genomes in the "
            f"corpus, so a corpus-wide sweep evicts what it is about to need")
        self.assertEqual(DEFAULT_DIR.resolve(), CROSSWALK_DIR.resolve(),
                         "the default crosswalk directory is not the one the tests "
                         "pass, so `load(acc)` and `load(acc, CROSSWALK)` are two cache "
                         "entries reading two different sets of files")


class EveryInferredBoundarySaysSo(unittest.TestCase):
    """Round 9: `boundary_uncertainty_bp` was 0 on all five conventions.

    Three of them say in their own source field that the source prints no coordinates,
    so `near_a_boundary` returned False everywhere and the UNRELIABLE path was
    unreachable on the shipped data. The field was the wrong instrument: the intervals
    taken for Whisnant and Szpara are the record's own TRL and TRS extents, exact to the
    base, and what is uncertain is which repeat they are rather than where the edge sits.
    So each convention states whether its intervals were read or inferred, and the
    inference reaches the user.

    This gate ties the flag to the source text, which is the evidence actually held, so
    the flag cannot be set to `read` on a row whose source says the coordinates were not
    printed.
    """

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_every_inferred_row_quotes_the_words_that_justify_it(self):
        """The flag is tied to the row's own source text, by a quotation.

        Round 10: this gate matched four phrases of this repository's own prose, so
        rewording a source field honestly ("The paper gives none of these figures
        itself") turned the suite red on correct data, and the gate could not see a
        mislabelled row whose prose had been reworded either. The row quotes its
        evidence now and the quote has to be in the source.
        """
        wrong = []
        for c in self.cs.conventions:
            source = (c.source or "") + " " + (c.quote or "")
            if c.boundaries_are_inferred:
                if not c.boundaries_evidence:
                    wrong.append((c.name, "marked inferred and quotes no evidence"))
                elif c.boundaries_evidence not in source:
                    wrong.append((c.name, f"quotes {c.boundaries_evidence!r}, which is "
                                          f"not in its own source field"))
                elif len(c.boundaries_evidence) < 20:
                    wrong.append((c.name, f"quotes {c.boundaries_evidence!r}, which is "
                                          f"too short to be evidence of anything"))
                elif not re.search(r"\b(no|not|never|none)\b[^.;:]*"
                                   r"(coordinate|print|retriev|publish|state)|"
                                   r"(coordinate|print|retriev|publish|state)[^.;:]*"
                                   r"\b(no|not|never|none)\b",
                                   c.boundaries_evidence, re.I):
                    # Round 14: the alternation here was a list of words any of which
                    # would do, so seven clauses of the Morse source passed it and two
                    # of them argue the opposite of the flag. `the paper does name its
                    # coordinates` was legal evidence for "the intervals were not read
                    # from the source", and round 13 had just made this string
                    # user-facing in the report. An inference is justified by the
                    # source NOT giving the numbers, so the clause has to say that:
                    # a negation and a word about printing or retrieving them, in
                    # either order, within the clause.
                    # Round 11: the check was that the row quotes itself, so setting
                    # whisnant2020_masked's evidence to "Whisnant" passed. The quote has
                    # to be about where the intervals came from, not merely present.
                    wrong.append((c.name, f"quotes {c.boundaries_evidence!r}, which "
                                          f"says nothing about where the intervals "
                                          f"came from"))
                elif not _is_whole_clauses(c.boundaries_evidence, source):
                    # Round 12, the fourth version of this check to be gamed. "At least
                    # 20 characters" and a keyword alternation whose weakest member is
                    # `no\b` were both satisfied by 'rative functional geno', a 22
                    # character window cut out of the middle of the paper's own title,
                    # with "no" matching at the truncation. The suite stayed green and
                    # `usable(118777, 'whisnant2020_masked')` went on reporting the
                    # inference caveat as if it were justified.
                    #
                    # A quotation is now required to be whole clauses of the source, so
                    # a window cut across word and clause boundaries fails whatever
                    # letters it happens to contain. That is a property of the text
                    # rather than a list of words to avoid, which is what the three
                    # previous versions were.
                    wrong.append((c.name, f"quotes {c.boundaries_evidence!r}, which is "
                                          f"a window cut out of the source rather than "
                                          f"whole clauses of it"))
            elif c.boundaries_evidence:
                wrong.append((c.name, "quotes evidence for an inference it does not "
                                      "claim"))
        self.assertEqual(wrong, [], f"{len(wrong)} conventions disagree with their own "
                                    f"source about where their intervals came from: "
                                    f"{wrong}")

    def test_the_flag_and_the_basis_cannot_be_edited_apart(self):
        """Round 15: the `read` flag rested on columns this repository writes.

        A reviewer rewrote one clause of `source` to dodge round 14's contradiction
        check, appended a fabricated sentence to `quote`, flipped `boundaries` to
        `read` and cleared the evidence: `usable(118777, 'whisnant2020_masked')`
        stopped reporting the inference, and the only failures were two README
        quotations, which the reviewer then updated as a developer would. Cause 49's
        own lesson is that a burden checked against our own prose is not checked.

        So the basis is a controlled value, the flag is derived from it, and both are
        in the artefact: editing either alone fails here, and what a reader is owed is
        stated per basis and printed by `report`, for a transcription as well as for
        an inference.
        """
        from hsv1coords.conventions import Convention
        rows = list(csv.DictReader(
            (ROOT / "data" / "conventions.tsv").open(), delimiter="\t"))
        self.assertTrue(rows)
        wrong = []
        for row in rows:
            basis = (row.get("interval_basis") or "").strip()
            if basis not in Convention.INTERVAL_BASES:
                wrong.append((row["convention"], f"basis {basis!r} is not one of "
                                                 f"{sorted(Convention.INTERVAL_BASES)}"))
                continue
            reads, caveat = Convention.INTERVAL_BASES[basis]
            stated = (row.get("boundaries") or "").strip()
            if stated != ("read" if reads else "inferred"):
                wrong.append((row["convention"],
                              f"basis {basis!r} means {'read' if reads else 'inferred'} "
                              f"and the boundaries column says {stated!r}"))
            if not caveat.strip():
                wrong.append((row["convention"], f"basis {basis!r} owes a reader "
                                                 f"nothing, so it says nothing"))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} rows disagree with their own interval basis: {wrong}")
        # Both halves present, so neither column is a constant.
        bases = {(r.get("interval_basis") or "").strip() for r in rows}
        self.assertGreater(len(bases), 2,
                           f"the basis column holds only {bases}, so it carries "
                           f"little more than the flag it replaced")
        # Round 16: the obligations were a loop over the rows that HAVE a basis, so
        # moving the one `paper_transcribed` row to `deposited_record` left nothing
        # asserting anything about it: one cell deleted round 15's whole disclosure,
        # and `deposited_record` carried no burden at all. That is the
        # filtered-iteration shape of cause 6b inside the fix for cause 68.
        #
        # So the row set per basis is recorded, both directions, and each basis owes
        # something checkable:
        by_basis: dict[str, list[str]] = {}
        for row in rows:
            by_basis.setdefault(row["interval_basis"].strip(), []).append(
                row["convention"])
        self.assertEqual(
            {k: sorted(v) for k, v in by_basis.items()},
            {"deposited_record": ["full_length"],
             "paper_transcribed": ["lassalle2020_excluded"],
             "inferred_from_record": ["szpara2014_trimmed", "whisnant2020_masked"],
             "derived_from_stated_rule": ["morse2017_nr"]},
            "the interval basis of a convention has changed. That is the flag that "
            "decides what a reader is told about where its numbers came from, so the "
            "change is deliberate and belongs here as well as in the artefact.")
        report = str(self.cs.report(118777, RepeatMap.load()))
        for name in by_basis.get("deposited_record", []):
            with self.subTest(convention=name, basis="deposited_record"):
                # Its own caveat says "the convention is the deposited record
                # itself", which is only true of a convention that states no
                # intervals of its own.
                self.assertEqual(self.cs[name].intervals, (),
                                 "a convention whose basis is the deposited record "
                                 "states intervals of its own, so the record is not "
                                 "what its coordinates rest on")
                self.assertEqual(self.cs[name].records, ())
        for name in by_basis.get("paper_transcribed", []):
            with self.subTest(convention=name, basis="paper_transcribed"):
                self.assertTrue(self.cs[name].intervals,
                                "a transcription with nothing transcribed")
                self.assertIn("transcription:", report,
                              "the report does not disclose the transcription")
                self.assertIn(name, report)
        for name in by_basis.get("inferred_from_record", []) + by_basis.get(
                "derived_from_stated_rule", []):
            with self.subTest(convention=name, basis="inferred"):
                self.assertTrue(self.cs[name].boundaries_evidence,
                                "an inferred interval quotes no evidence")
                self.assertIn(name, report)

    def test_every_read_row_shows_its_source_printing_the_coordinates(self):
        """Round 13: the gate above guards only the `inferred` branch.

        Everything that checked this flag asked what an `inferred` row owed, so the
        flag could be answered by writing `read`. Flipping `whisnant2020_masked` to
        `read` and clearing its evidence left every test in this class green, and
        `usable(118777, 'whisnant2020_masked')` stopped reporting the inference
        caveat: the only tests that went red were the README output quotations, which
        catch it by coincidence because the README happens to quote a provenance line
        that names the row. Reword the README and the mislabelling ships.

        So `read` carries its own burden, and it is the burden the word means: the
        source prints these numbers. Every edge of a read row's intervals and records
        has to appear as a numeral in that row's verbatim quote of its source, with or
        without thousands separators. The row's `source` field is this repository's own
        prose about the source and is deliberately not consulted here, because quoting
        our own commentary is how the three previous versions of the sibling check were
        gamed.
        """
        wrong = []
        for c in self.cs.conventions:
            if c.boundaries_are_inferred:
                continue
            edges = [e for a, b in c.intervals for e in (a, b)]
            edges += [e for _, a, b in c.records for e in (a, b)]
            if not edges:
                # `full_length` names no intervals: the convention is the deposited
                # record as it stands, so there are no coordinates to have been read
                # from anywhere, and the record itself is the source.
                continue
            # Round 14: the burden rests on `quote`, which is a column nothing
            # verifies against the paper, so appending one fabricated sentence to it
            # ("The masked terminal repeats are 1-9,213 and 145,590-152,222.") bought
            # a `read` flag and dropped the inference caveat out of `usable()`. The
            # row's own `source` field is deliberately not the evidence here, because
            # quoting our own prose is how the sibling check was gamed three times,
            # but it is the one field that CONTRADICTS the lie, and discarding it
            # discarded that too. A row cannot say its coordinates were read from a
            # source its own prose says prints none.
            contradiction = re.search(
                r"\b(no|not|never|none)\b[^.;:]*(coordinate|print|retriev)|"
                r"(coordinate|print|retriev)[^.;:]*\b(no|not|never|none)\b",
                c.source or "", re.I)
            if contradiction:
                wrong.append((c.name, f"says its intervals were read from the source "
                                      f"and its own source field says "
                                      f"{contradiction.group(0)!r}"))
                continue
            printed = {int(n.replace(",", ""))
                       for n in re.findall(r"\d[\d,]*", c.quote or "")}
            missing = sorted(e for e in edges if e not in printed)
            if missing:
                wrong.append((c.name, f"says its intervals were read from the source "
                                      f"and its quote of that source does not print "
                                      f"{missing}"))
        self.assertEqual(wrong, [], f"{len(wrong)} conventions claim coordinates read "
                                    f"from a source that does not print them: {wrong}")

    def test_the_inference_reaches_a_caller_of_usable(self):
        """Round 10: the docstring said it reaches `usable`. It reached `report` only."""
        for c in self.cs.conventions:
            if not c.boundaries_are_inferred:
                continue
            with self.subTest(convention=c.name):
                _, why = self.cs.usable(118777, c.name)
                self.assertIn("interval not read from the source", why,
                              f"{c.name} takes its intervals from somewhere other than "
                              f"its source's own coordinates and usable() says nothing "
                              f"about it")

    def test_a_read_convention_says_nothing_of_the_kind(self):
        _, why = self.cs.usable(118777, "full_length")
        self.assertNotIn("interval not read from the source", why)

    def test_both_kinds_are_present_so_the_flag_is_not_a_constant(self):
        kinds = {c.boundaries for c in self.cs.conventions}
        self.assertEqual(kinds, {"read", "inferred"},
                         f"the boundaries column holds only {kinds}, so it is a "
                         f"constant and carries no information")

    def test_an_inferred_convention_is_named_in_the_report(self):
        """Once per report, not once per convention per position.

        The first version of this attached a caveat per convention, which put three
        paragraphs of provenance on every position in the genome and buried the caveat
        that was about the position. Two existing tests caught that, which is what they
        are for.
        """
        inferred = sorted(c.name for c in self.cs.conventions
                          if c.boundaries_are_inferred)
        self.assertTrue(inferred, "no convention is marked inferred")
        report = self.cs.report(9213)
        self.assertEqual(sorted(report.inferred_intervals), inferred)
        provenance = [line for line in str(report).splitlines()
                      if line.strip().startswith("provenance:")]
        self.assertEqual(len(provenance), 1,
                         f"the report states its interval provenance {len(provenance)} "
                         f"times; it should state it once")
        for name in inferred:
            with self.subTest(convention=name):
                self.assertIn(name, provenance[0])

    def test_a_read_convention_is_not_named_as_inferred(self):
        """The control: the provenance line must not list everything."""
        report = self.cs.report(9213)
        for c in self.cs.conventions:
            if c.boundaries_are_inferred:
                continue
            with self.subTest(convention=c.name):
                self.assertNotIn(c.name, report.inferred_intervals,
                                 f"{c.name} prints its own coordinates and is listed as "
                                 f"inferred anyway")

    def test_the_provenance_line_does_not_crowd_out_the_position_caveats(self):
        """A unique position gets no caveat; an ambiguous one still gets its own."""
        unique = self.cs.report(60000, RepeatMap.load())
        self.assertEqual(unique.caveats, [],
                         f"a well defined position carries caveats: {unique.caveats}")
        ambiguous = self.cs.report(9213, RepeatMap.load())
        self.assertTrue(any("not well defined" in c for c in ambiguous.caveats),
                        f"the ambiguity caveat is missing: {ambiguous.caveats}")


class TheIdentityConversionIsBounded(unittest.TestCase):
    """Round 9: `convert(999999, NEW, NEW)` returned 999999, and `convert(-3, ...)` -3.

    The version names were validated before the identity short-circuit, which is round
    7's fix for the same function, and the position was not validated at all. Fixing the
    lines a blocker quoted rather than the whole function is PITFALLS root cause 14.

    Written after the fix, because reverting the fix left the suite green: a fix with no
    test is a fix that will be reverted by the next person who tidies the branch.
    """

    @classmethod
    def setUpClass(cls):
        from hsv1coords.versions import VersionCrosswalk
        cls.v = VersionCrosswalk.load()

    def test_a_position_past_the_end_does_not_survive_the_identity_case(self):
        from hsv1coords.conventions import VERSION_LENGTHS
        for version, length in sorted(VERSION_LENGTHS.items()):
            with self.subTest(version=version):
                self.assertIsNone(self.v.convert(length + 1, version, version))
                self.assertIsNone(self.v.convert(999999, version, version))
                self.assertIsNone(self.v.convert(0, version, version))
                self.assertIsNone(self.v.convert(-3, version, version))
                self.assertEqual(self.v.convert(length, version, version), length,
                                 "the last base of the version must still convert")

    def test_a_fractional_or_boolean_position_is_refused(self):
        with self.assertRaises(TypeError):
            self.v.convert(1.5, "NC_001806.2", "NC_001806.2")
        with self.assertRaises(TypeError):
            self.v.convert(True, "NC_001806.2", "NC_001806.2")


class TheOutOfRangeGuardIsTested(unittest.TestCase):
    """Round 9: removing the bounds check in `_usable_native` left the suite green.

    The mutant answered `usable(999999, 'full_length')` with "not in any affected
    region", which is True, and reported a position 847,777 bp past the end of the
    genome as usable under four conventions. The shipped code was right; only the guard
    was uncovered.
    """

    @classmethod
    def setUpClass(cls):
        cls.cs = ConventionSet.load()

    def test_a_position_past_the_end_is_refused_under_every_convention(self):
        for name in sorted(c.name for c in self.cs.conventions):
            with self.subTest(convention=name):
                ok, why = self.cs.usable(999999, name)
                self.assertFalse(ok, f"{name} called 999999 usable: {why}")
                # `morse2017_nr` is defined against NC_001806.1, so the version
                # conversion declines before the bounds check is reached. Both are
                # refusals; only the wording differs.
                self.assertRegex(why, r"outside|no coordinate in")

    def test_zero_and_negative_positions_are_refused(self):
        for position in (0, -1, -152222):
            with self.subTest(position=position):
                ok, why = self.cs.usable(position, "full_length")
                self.assertFalse(ok)
                self.assertIn("outside", why)

    def test_the_last_base_of_each_version_is_still_accepted(self):
        """The control: the guard must not refuse the genome's own last base."""
        ok, why = self.cs.usable(152222, "full_length")
        self.assertTrue(ok, f"152222 is the last base of NC_001806.2: {why}")

    def test_the_report_does_not_list_an_impossible_position_as_usable(self):
        report = str(self.cs.report(999999))
        offending = [line.strip() for line in report.splitlines()
                     if "usable under" in line and "NOT usable" not in line]
        self.assertEqual(offending, [],
                         f"the report calls an impossible position usable: {offending}")



def _clauses(text: str) -> set[str]:
    """The source's clauses, split at sentence and clause punctuation.

    Round 13: the split was on any `[.;:,]`, and these source fields are full of
    accession versions and thousands separators, so `NC_001806.1` and `9,338`
    manufactured clause boundaries in the middle of a sentence. That let a 40 character
    window, `2 coordinates under a declared NC_001806`, cut out of "the first stored
    NC_001806.2 coordinates under a declared NC_001806.1 frame", satisfy every one of
    the four filters: present in the source, long enough, mentions `coordinate`, and
    "whole clauses". The fifth version of this check to be walked past.

    Clause punctuation is followed by a space or the end of the text. The punctuation
    inside `NC_001806.1`, `9,338` and `150,962` is not.
    """
    return {part.strip() for part in re.split(r"[.;:,](?=\s|$)", text) if part.strip()}


def _is_whole_clauses(quote: str, source: str) -> bool:
    """Whether a quotation is made of whole clauses of its source.

    Every clause of the quote has to be a clause of the source. A character window cut
    across word boundaries ('rative functional geno') produces a fragment that is a
    clause of nothing, which is what round 12 used to pass the previous three versions
    of this check.
    """
    available = _clauses(source)
    parts = _clauses(quote)
    return bool(parts) and parts <= available


class ThePackageCarriesItsOwnData(unittest.TestCase):
    """Round 11: `pip install .` gave an importable library in which nothing worked.

    Every module resolved its data as `parents[2] / "data"`, the repository layout, so in
    a non-editable install each loader raised FileNotFoundError and
    `Crosswalk.available()` returned an empty list, which made a corpus sweep a silent
    no-op. `pip install -e .` worked, which is what README.md recommends, so this hit
    only the user who installed it as a dependency.
    """

    def test_the_loaders_resolve_through_one_place(self):
        """Six modules computed the same path six times; they share it now."""
        import hsv1coords._paths as paths
        source = pathlib.Path(paths.__file__).parent
        offenders = []
        for module in sorted(source.glob("*.py")):
            if module.name == "_paths.py":
                continue
            text = module.read_text()
            # No exemptions. Round 12: `facts.py` was exempt, on the grounds that it
            # reads data/raw as well and is a build helper, and it was computing
            # `parents[2] / "data"` for BOTH directories, so in an installed copy
            # `derive()` returned an empty list and every figure went ungated instead of
            # failing. An exemption from the rule that found a defect five times is
            # where the sixth one lives.
            # Code, not prose: a comment or docstring recording the defect is not the
            # defect, and this check is textual.
            code = "\n".join(line.split("#")[0] for line in text.splitlines()
                             if not line.strip().startswith("#"))
            if "parents[2]" in code and "data" in code:
                offenders.append(module.name)
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} modules still compute a data path from the repository "
            f"layout, which does not hold in an installed copy: {offenders}")

    def test_what_a_wheel_carries_is_derived_from_what_the_loaders_read(self):
        """Round 12: the staged list was a tuple, and this test was membership-only.

        It asserted that five file names appeared in `staging.WANTED` or
        `staging.DERIVED`, so dropping `annotation_v0.1.0.tsv` or the block files from
        those tuples passed, and it never looked at whether a staged copy existed or was
        current. `accessions.txt` was in the tuple and matched by no `package-data` glob,
        so it was staged and then not shipped while the script printed success.

        The list is derived now: `_paths.requested()` records every data file the
        library asks for, and the staging script stages exactly those. This asserts the
        derivation rather than a typed list, by exercising the loaders in this process
        and requiring each path they read to be a real file or directory under `data/`.
        """
        sys.path.insert(0, str(ROOT / "src"))
        import stage_package_data as staging
        read = staging.what_the_loaders_read()
        self.assertGreaterEqual(
            len(read), 5,
            f"the loaders report reading only {len(read)} data paths: {read}")
        for relative in read:
            with self.subTest(path=relative):
                source = ROOT / "data" / pathlib.Path(relative)
                self.assertTrue(
                    source.exists(),
                    f"a loader reads data/{relative}, which does not exist, so a wheel "
                    f"cannot carry it")
        globs = staging.shipped_globs()
        unshipped = [r for r in read
                     if not any(fnmatch.fnmatch(f"data/{r}/x.tsv", g)
                                or fnmatch.fnmatch(f"data/{r}", g) for g in globs)]
        self.assertEqual(
            unshipped, [],
            f"{len(unshipped)} paths the loaders read match no package-data glob in "
            f"pyproject.toml, so a wheel would be built without them: {unshipped}")

    def test_the_staged_copy_if_present_is_identical_to_the_source(self):
        """Round 12: the staged copy was stale and nothing could see it.

        `src/hsv1coords/data/aliases_literature.tsv` held 220 rows claiming
        `verbatim_in_citation=yes` against 205 in the checkout: the exact 15 rows round
        11 corrected were still asserted as attested in the shipped library. Nothing
        read it, because `_paths.data_dir()` prefers the checkout and the staged tree is
        gitignored, so neither the suite nor CI ever opened it.
        """
        staged = ROOT / "src" / "hsv1coords" / "data"
        if not staged.is_dir():
            # Round 13: this skipped, and the staged tree is gitignored, so the skip
            # was unconditional on every clone and in any run that had not staged. The
            # defect it exists for is a STALE staged copy, which is only visible to a
            # reader who has one, so an absent one is a reason to build it rather than
            # to pass. `make all` stages before it tests.
            if (ROOT / ".git").exists():
                self.fail("no staged copy in a checkout, so nothing compares what a "
                          "wheel would ship with what these tests pass against. Run "
                          "`make package-data`.")
            self.skipTest("no staged copy and no checkout to build one from")
        differing, checked = [], 0
        for path in sorted(staged.rglob("*")):
            if not path.is_file():
                continue
            source = ROOT / "data" / path.relative_to(staged)
            checked += 1
            if not source.exists():
                differing.append(f"{path.relative_to(staged)} (no source)")
            elif source.read_bytes() != path.read_bytes():
                differing.append(str(path.relative_to(staged)))
        self.assertEqual(
            differing, [],
            f"{len(differing)} staged files are not what data/ holds, so a wheel built "
            f"now would ship a different library from the one these tests pass "
            f"against: {differing[:5]}")
        # Round 13: this was `assertGreater(checked, 150)`, a number typed beside the
        # thing it measures, so a third of the crosswalk could go unstaged and
        # uncompared (standing rule 47). The count is derived from the paths the
        # loaders read: every file under each of them has to be in the staged tree.
        # In a subprocess, because `_paths.requested()` is per process and accumulates:
        # anything else this suite has imported (`facts`, which names the whole of
        # `data/derived`) would otherwise be counted as something the staging script
        # stages, which it is not. The script runs alone, so the derivation is measured
        # alone. Round 13.
        import subprocess as _subprocess
        reader = _subprocess.run(
            [sys.executable, "-I", "-c",
             "import sys; sys.path.insert(0, 'src');"
             "import stage_package_data as s;"
             "print('\\n'.join(s.what_the_loaders_read()))"],
            cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(reader.returncode, 0, reader.stderr[-800:])
        expected: set[str] = set()
        for relative in reader.stdout.split():
            source = ROOT / "data" / pathlib.Path(relative)
            if source.is_dir():
                expected.update(f"{relative}/{p.name}" for p in source.iterdir()
                                if p.is_file())
            elif source.is_file():
                expected.add(relative)
        present = {str(p.relative_to(staged)) for p in staged.rglob("*") if p.is_file()}
        self.assertEqual(
            sorted(expected - present), [],
            f"{len(expected - present)} files the loaders read are not in the staged "
            f"tree, so a wheel built now would not carry them: "
            f"{sorted(expected - present)[:5]}")
        self.assertEqual(
            checked, len(present),
            "the comparison above did not read every staged file")

    def test_an_absent_crosswalk_directory_raises_rather_than_answering(self):
        from hsv1coords import Crosswalk
        with self.assertRaises(FileNotFoundError):
            Crosswalk.available(ROOT / "data" / "derived" / "not-a-directory")
        self.assertGreater(len(Crosswalk.available(CROSSWALK_DIR)), 150)


if __name__ == "__main__":
    unittest.main(verbosity=2)
