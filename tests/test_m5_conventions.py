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
        return ConventionSet(base.conventions + [row])

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
        for c in cs.conventions:
            with self.subTest(convention=c.name):
                self.assertTrue(c.quote.strip())
                self.assertTrue(c.source.strip())
