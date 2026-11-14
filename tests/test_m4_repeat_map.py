"""The M4 gate, written before the code, per ROADMAP.md.

The last criterion is deliberate: this milestone must demonstrably cover the three
defects that review found in M3, rather than merely claiming to.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import (Crosswalk, Dictionary, RepeatMap, lift,  # noqa: E402
                        lift_detailed)
from hsv1coords.liftover import REFERENCE  # noqa: E402
from hsv1coords.repeats import GENOME_BP  # noqa: E402
from build_crosswalk import genbank_sequence  # noqa: E402

INVERTED_PAIRS = RepeatMap.load().inverted_pairs

DERIVED = ROOT / "data" / "derived"


def require_reference() -> pathlib.Path:
    """The reference record is fetched, not committed, so skip without it.

    A fresh clone must skip a test that needs it rather than error: round 7 found one
    erroring there, which is the same class of fault as a test that cannot fail.
    """
    ref = ROOT / "data" / "raw" / "ref_NC_001806.2.gb"
    if not ref.exists():
        raise unittest.SkipTest(
            "data/raw/ref_NC_001806.2.gb is fetched, not committed. "
            "Run `make inputs` first.")
    return ref

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

    def test_counterpart_declines_wherever_there_is_more_than_one(self):
        """The library's central contract, which had no test until round 7.

        Changing `len(alts) == 1` to `len(alts) >= 1` made `counterpart` pick the first
        of two or three equally valid answers and changed 32 live answers with the whole
        suite green. Every position with several counterparts must get None.
        """
        several = [p for p in range(1, GENOME_BP + 1, 23)
                   if len(self.m.counterparts(p)) > 1]
        self.assertGreater(len(several), 50,
                           "no position has several counterparts, so this proves "
                           "nothing; the map or the sampling changed")
        for pos in several:
            with self.subTest(position=pos):
                self.assertIsNone(
                    self.m.counterpart(pos),
                    f"{pos} has {self.m.counterparts(pos)} and must not be given one")
        exactly_one = [p for p in range(1, GENOME_BP + 1, 23)
                       if len(self.m.counterparts(p)) == 1]
        self.assertGreater(len(exactly_one), 50, "nothing exercises the single case")
        for pos in exactly_one[:200]:
            with self.subTest(single=pos):
                self.assertEqual(self.m.counterpart(pos),
                                 self.m.counterparts(pos)[0])

    def test_no_repeat_geometry_is_hand_listed_in_the_library(self):
        """The module's docstring says nothing is hand-listed. Round 7: it was.

        The three 'a' sequence copies with their orientations, and the inverted pair
        names, were literals in `hsv1coords/repeats.py`. Both are derivable from the
        artefact, so a literal reappearing there is a defect this catches.
        """
        body = (ROOT / "src" / "hsv1coords" / "repeats.py").read_text()
        for banned in ("A_SEQUENCE_COPIES", "INVERTED_PAIRS", "125975", "151824",
                       '"TRL", "IRL"', '"IRS", "TRS"'):
            with self.subTest(literal=banned):
                self.assertNotIn(banned, body,
                                 f"{banned} is hand-listed in repeats.py; derive it "
                                 f"in src/build_repeat_map.py instead")

    def test_the_derived_orientations_match_the_sequence(self):
        """The orientation column decides every counterpart, so it is checked here.

        Not against the builder that wrote it, which would be circular, but against
        the reference sequence itself.
        """
        ref = genbank_sequence(require_reference()).upper()
        comp = str.maketrans("ACGTN", "TGCAN")
        groups: dict[str, list] = {}
        for r in self.m.regions:
            groups.setdefault(r.name, []).append(r)
        checked = 0
        for name, copies in groups.items():
            if len(copies) < 2 or len({r.length for r in copies}) != 1:
                continue
            first = ref[copies[0].start - 1:copies[0].end]
            for r in copies:
                here = ref[r.start - 1:r.end]
                rc = here.translate(comp)[::-1]
                fwd = sum(1 for a, b in zip(here, first) if a == b)
                rev = sum(1 for a, b in zip(rc, first) if a == b)
                with self.subTest(region=name, start=r.start):
                    expected = "forward" if fwd > rev else "reverse"
                    self.assertEqual(r.orientation, expected,
                                     f"{name} at {r.start} is recorded "
                                     f"{r.orientation} but matches {expected} "
                                     f"({fwd} forward against {rev} reverse)")
                checked += 1
        self.assertGreaterEqual(checked, 3,
                                "no multi-copy region was checked, so this is empty")

    def test_the_inverted_repeats_are_reverse_complementary_base_by_base(self):
        """The premise every counterpart rests on, checked in code at last.

        ROADMAP.md has stated 15,846 of 15,846 reverse-complementary since M4, and the
        claim is true, but nothing executable re-asserted it: the inversion arithmetic
        in `repeats.py` is the whole product and its premise lived only in prose. One
        commit away from being silently wrong.
        """
        ref = genbank_sequence(require_reference()).upper()
        comp = str.maketrans("ACGTN", "TGCAN")
        total = matched = 0
        for a, b in INVERTED_PAIRS:
            ra = next(r for r in self.m.regions
                      if r.name == a and r.kind == "repeat_region")
            rb = next(r for r in self.m.regions
                      if r.name == b and r.kind == "repeat_region")
            first = ref[ra.start - 1:ra.end]
            second = ref[rb.start - 1:rb.end].translate(comp)[::-1]
            self.assertEqual(len(first), len(second),
                             f"{a} and {b} differ in length")
            total += len(first)
            matched += sum(1 for x, y in zip(first, second) if x == y)
        self.assertEqual(matched, total,
                         f"only {matched} of {total} bases are "
                         f"reverse-complementary")
        self.assertEqual(total, 15846,
                         f"the two pairs cover {total} bases, not the 15,846 the "
                         f"documents state")

    def test_a_coordinate_must_be_an_integer(self):
        """A fractional position produced fractional counterparts.

        `classify(1.5)` reported "could equally be 126372.5", which is not a position.
        A string reached a comparison deep inside and raised there instead of at the
        boundary.
        """
        for bad in (1.5, "1", None, True, [1]):
            with self.subTest(value=bad):
                with self.assertRaises(TypeError):
                    self.m.classify(bad)

    def test_out_of_range_is_distinguishable_from_ambiguous(self):
        """Both were `well_defined=False` with the difference only in the prose.

        `LiftResult.well_defined` is fed from this field, so a caller branching on it
        treated a coordinate that does not exist as an in-range ambiguous one.
        """
        missing = self.m.classify(999999)
        ambiguous = self.m.classify(1)
        self.assertFalse(missing.in_range)
        self.assertTrue(ambiguous.in_range)
        self.assertFalse(missing.well_defined)
        self.assertFalse(ambiguous.well_defined)
        self.assertEqual(missing.alternatives, [])
        self.assertTrue(ambiguous.alternatives)
        for pos in (1, GENOME_BP, GENOME_BP // 2):
            with self.subTest(position=pos):
                self.assertTrue(self.m.classify(pos).in_range)
        for pos in (0, -1, GENOME_BP + 1):
            with self.subTest(position=pos):
                self.assertFalse(self.m.classify(pos).in_range)

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

class EveryIntervalComesFromTheRecord(unittest.TestCase):
    """Round 9: 30 of the 34 regions and every diploid pair were tied to nothing.

    `test_the_four_inverted_repeats` hard-codes the four large boundaries, and that was
    the whole of it. A reviewer moved `oriL` 20 kb, moved an ICP0 copy 90 kb, and shortened
    a LAT copy by one base, each with all 180 tests green, while `repeats.py` claimed
    "every region, every pair and every copy's orientation is derived from the reference
    record's own feature table and sequence". The claim was true of the builder and
    untrue of the artefact, which is the file the library actually reads.

    This ties each interval back to a feature in the record. It is not a reimplementation
    of the builder's matching logic, which would drift from it (root cause 13): it asks
    only whether every interval in the artefact is an interval the record states, and
    whether each pair's two copies sit in the two copies of their repeat.
    """

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        from build_repeat_map import spans_of
        record = require_reference()
        cls.spans = spans_of
        cls.features = parse_genbank(record)
        cls.extents = set()
        for f in cls.features:
            spans = spans_of(f.get("_location", ""))
            if not spans:
                continue
            cls.extents.add((min(s for s, _ in spans), max(e for _, e in spans)))
            for span in spans:
                cls.extents.add(span)
        cls.m = RepeatMap.load(DERIVED)

    def test_every_region_interval_is_one_the_record_states(self):
        strays = [(r.name, r.start, r.end) for r in self.m.regions
                  if (r.start, r.end) not in self.extents]
        self.assertEqual(
            strays, [],
            f"{len(strays)} regions in repeat_map.tsv have an interval that is not a "
            f"feature extent in {REFERENCE}: {strays[:8]}. Either the artefact was "
            f"edited or the record changed; rebuild with `make repeatmap`.")
        self.assertGreaterEqual(len(self.m.regions), 30,
                                "the region count collapsed, so this gate is checking "
                                "almost nothing")

    def test_every_diploid_copy_interval_is_one_the_record_states(self):
        strays = []
        for pair in self.m.pairs:
            for label, copy in (("copy_1", pair.copy_1), ("copy_2", pair.copy_2)):
                if tuple(copy) not in self.extents:
                    strays.append((pair.locus, label, tuple(copy)))
        self.assertEqual(
            strays, [],
            f"{len(strays)} diploid-pair copies have an interval that is not a feature "
            f"extent in {REFERENCE}: {strays[:8]}")
        self.assertGreaterEqual(len(self.m.pairs), 25,
                                "the pair count collapsed, so this gate is checking "
                                "almost nothing")

    def test_each_pair_sits_in_the_two_copies_of_its_repeat(self):
        """A copy that moved out of its repeat is the 90 kb mutation the reviewer made."""
        repeats = {label: next(r for r in self.m.regions
                               if r.name.startswith(label)
                               and r.kind == "repeat_region")
                   for label in ("TRL", "IRL", "IRS", "TRS")}
        misplaced = []
        for pair in self.m.pairs:
            first, second = pair.repeat_pair.split("/")
            for copy, label in ((pair.copy_1, first), (pair.copy_2, second)):
                r = repeats.get(label)
                if r is None:
                    continue
                # Overlap, not containment: the LAT copy in the long repeat runs 778 bp
                # past the end of IRL, which is why containment is the wrong test here.
                if not (copy[0] <= r.end and copy[1] >= r.start):
                    misplaced.append((pair.locus, label, tuple(copy),
                                      (r.start, r.end)))
        self.assertEqual(
            misplaced, [],
            f"{len(misplaced)} diploid copies sit outside the repeat they are recorded "
            f"against (locus, repeat, copy, repeat interval): {misplaced[:8]}")

    def test_the_two_copies_of_a_pair_are_equal_length_unless_marked_partial(self):
        wrong = []
        for pair in self.m.pairs:
            a = pair.copy_1[1] - pair.copy_1[0]
            b = pair.copy_2[1] - pair.copy_2[0]
            partial = "partial" in (pair.note or "").lower()
            if a != b and not partial:
                wrong.append((pair.locus, a + 1, b + 1))
            if a == b and partial:
                wrong.append((pair.locus, "marked partial but equal length", a + 1))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} pairs state two copies of different length without being "
            f"marked partial, or the reverse: {wrong[:8]}")


class EveryPublishedPairStatesItsStrand(unittest.TestCase):
    """Round 10: the diploid artefacts gave two copies with no strand.

    `diploid_pairs.tsv` and the table in `repeat_map.md` list Copy 1 and Copy 2 beside an
    explanation saying a position in one maps to the other by the repeat's own
    arithmetic. For a pair in the large inverted repeats that arithmetic reverses, and
    the base at the paired position is the complement. The same file also called TRL and
    IRL "the same 9,213 bp of sequence in two copies", which omits the inversion their
    names denote.
    """

    @classmethod
    def setUpClass(cls):
        cls.rows = list(csv.DictReader(
            (DERIVED / "diploid_pairs.tsv").open(), delimiter="\t"))
        cls.published = (DERIVED / "repeat_map.md").read_text()

    def test_every_pair_states_a_strand_relationship(self):
        missing = [r["locus"] for r in self.rows
                   if r.get("strand_relationship") not in
                   ("inverted", "direct", "undetermined")]
        self.assertEqual(missing, [],
                         f"{len(missing)} pairs state no strand relationship: "
                         f"{missing[:6]}")

    def test_an_equal_length_pair_of_real_size_is_measured_not_undetermined(self):
        """`undetermined` is allowed only where it cannot be measured.

        That is a shared length under 12 bp, where a match could be chance, or two
        copies of different length, where no constant offset relates them at all. Any
        other undetermined row would mean the measurement silently failed.
        """
        unexplained = []
        for r in self.rows:
            if r["strand_relationship"] != "undetermined":
                continue
            length_1 = int(r["copy_1_end"]) - int(r["copy_1_start"]) + 1
            length_2 = int(r["copy_2_end"]) - int(r["copy_2_start"]) + 1
            if length_1 != length_2 or min(length_1, length_2) < 12:
                continue
            unexplained.append((r["locus"], length_1, length_2))
        self.assertEqual(
            unexplained, [],
            f"{len(unexplained)} pairs are of equal, measurable length and still state "
            f"an undetermined strand: {unexplained[:6]}")

    def test_the_measurable_pairs_in_the_large_repeats_are_inverted(self):
        """The biology: the copies in TRL/IRL and IRS/TRS read in opposite directions."""
        measured = [r for r in self.rows
                    if r["strand_relationship"] in ("inverted", "direct")]
        self.assertGreater(len(measured), 10,
                           f"only {len(measured)} pairs could be measured, so this "
                           f"proves little")
        wrong = [(r["locus"], r["repeat_pair"], r["strand_relationship"])
                 for r in measured if r["strand_relationship"] != "inverted"]
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} pairs in the large inverted repeats are not inverted, which "
            f"would contradict the repeat structure: {wrong[:6]}")

    def test_the_published_table_carries_the_strand_column(self):
        self.assertIn("| Strand |", self.published,
                      "the published diploid table states no strand")
        self.assertIn("inverted", self.published)

    def test_the_published_text_does_not_call_inverted_copies_the_same_sequence(self):
        self.assertNotIn("are the same 9,213 bp of sequence in two copies",
                         self.published,
                         "the published text still describes the inverted repeats as "
                         "two copies of the same sequence without saying they are "
                         "inverted")
        self.assertIn("inverted", self.published)


class EveryQuotedIntervalIsOneTheRecordStates(unittest.TestCase):
    """The documentation gate exempts coordinate ranges. This closes the exemption.

    `test_no_occurrence_of_a_large_figure_is_out_of_the_gate_s_reach` skips `a..b`
    ranges, on the stated grounds that the repeat intervals are gated against the record
    in this file. They were gated as artefact rows, not as the text of a document, so
    round 9 changed `117161..126373` in README.md to `126673` with the suite green: a
    wrong boundary for IRL, contradicting both the README's own arithmetic and
    repeat_map.tsv, in the sentence that introduces the inverted repeats.

    Every such range in the four documents must be an interval the reference record
    states, either as a feature's extent or as one span of a joined location, which is
    what makes the LAT spans quoted in FINDINGS.md legitimate.
    """

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        from build_repeat_map import spans_of
        record = require_reference()
        cls.extents = set()
        for f in parse_genbank(record):
            spans = spans_of(f.get("_location", ""))
            if not spans:
                continue
            cls.extents.add((min(s for s, _ in spans), max(e for _, e in spans)))
            cls.extents.update(spans)

    NUMBER = r"(?:\d{1,3}(?:,\d{3})+|\d+)"
    RANGE = re.compile(rf"({NUMBER})\s*\.\.\s*({NUMBER})")
    """Any interval, not only four-digit ones.

    Round 10: requiring four or more digits on both sides left six quoted intervals
    unchecked, including `1..9213`, which is the TRL extent and the first interval the
    documents state. Changing it to `1..8999` passed. The thousands-separator branch has
    to be written out, because `[\d,]*` runs across the comma in a joined location such
    as `join(<1..4953,6910..7569)` and reads it as one 11-digit number.
    """

    NAMED_FEATURES = {
        "TRL": (1, 9213), "IRL": (117161, 126373),
        "IRS": (125975, 132607), "TRS": (145590, 152222),
    }
    """The four large repeats, whose extents a reader is most likely to act on.

    Checked for adjacency as well as membership: round 10 replaced the README's IRL
    interval with the extent of TRS, which is a real interval of the record, and the
    membership test passed. A wrong boundary on one of these four is the error this
    repository exists to prevent.
    """
    DOCS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")

    def test_every_interval_quoted_in_the_documents_is_in_the_record(self):
        wrong, checked = [], 0
        for name in self.DOCS:
            path = ROOT / name
            if not path.exists():
                continue
            for m in self.RANGE.finditer(path.read_text()):
                a = int(m.group(1).replace(",", ""))
                z = int(m.group(2).replace(",", ""))
                checked += 1
                if (a, z) not in self.extents:
                    wrong.append((name, f"{a}..{z}"))
        self.assertGreater(checked, 8,
                           f"only {checked} intervals were checked, so either the "
                           f"documents stopped quoting coordinates or the pattern no "
                           f"longer matches them")
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} intervals quoted in the documents are not intervals "
            f"{REFERENCE} states, as a feature extent or as a span of a joined "
            f"location: {wrong[:8]}")
        self.assertGreater(checked, 18,
                           f"only {checked} intervals were checked; the documents or "
                           f"the pattern changed")

    def test_an_interval_quoted_beside_a_repeat_name_is_that_repeat_s(self):
        """Membership is not enough: any real interval passed for any other.

        Round 10 gave IRL the extent of TRS in README.md and the suite stayed green.
        Both are real intervals of the record, so only adjacency to the name can tell
        them apart.
        """
        mismatched, checked = [], 0
        for name in self.DOCS:
            path = ROOT / name
            if not path.exists():
                continue
            body = re.sub(r"\s+", " ", path.read_text())
            for m in self.RANGE.finditer(body):
                a = int(m.group(1).replace(",", ""))
                z = int(m.group(2).replace(",", ""))
                before = body[max(0, m.start() - 80):m.start()]
                # The nearest repeat name before the interval, if any is close by.
                hits = [(before.rfind(label), label)
                        for label in self.NAMED_FEATURES
                        if label in before]
                if not hits:
                    continue
                _, label = max(hits)
                checked += 1
                if (a, z) != self.NAMED_FEATURES[label]:
                    # Another of the four named closer still, or a sub-interval of this
                    # one, are both legitimate; a DIFFERENT repeat's exact extent is not.
                    if (a, z) in self.NAMED_FEATURES.values():
                        mismatched.append((name, label, f"{a}..{z}",
                                           f"which is {self._which(a, z)}"))
        self.assertGreater(checked, 2,
                           f"only {checked} intervals sit beside a repeat name, so this "
                           f"gate is checking almost nothing")
        self.assertEqual(
            mismatched, [],
            f"{len(mismatched)} intervals are quoted beside the name of a different "
            f"repeat (document, name beside it, interval, what it really is): "
            f"{mismatched[:6]}")

    @classmethod
    def _which(cls, start: int, end: int) -> str:
        for label, extent in cls.NAMED_FEATURES.items():
            if extent == (start, end):
                return label
        return "not one of the four"


class ThePublishedTableMarksThePartialPairs(unittest.TestCase):
    """Round 9: `repeat_map.md` printed the seven partial pairs as if they were exact.

    The published table is the artefact a reader opens without the TSV beside it, and it
    omitted the `offset_delta` column and said nothing about length. So it printed LAT at
    1..7,569 against 118,777..127,151, two copies differing by 806 bp, in the same form
    as an exact pair. A reader transfers a LAT coordinate between them by a constant
    offset, which is the error this repository exists to prevent, on its own worked
    example.
    """

    def test_every_partial_pair_is_marked_partial_in_the_published_table(self):
        published = (DERIVED / "repeat_map.md").read_text()
        rows = [line for line in published.splitlines()
                if line.startswith("| `") and line.count("|") >= 6
                and ".." in line]
        self.assertGreater(len(rows), 20,
                           f"only {len(rows)} diploid rows found in repeat_map.md, so "
                           f"this gate is reading the wrong table")
        wrong = []
        for line in rows:
            cells = [c.strip() for c in line.strip("|").split("|")]
            try:
                a = [int(x.replace(",", "")) for x in cells[3].split("..")]
                b2 = [int(x.replace(",", "")) for x in cells[4].split("..")]
            except (ValueError, IndexError):
                continue
            differs = (a[1] - a[0]) != (b2[1] - b2[0])
            marked = "partial" in line.lower()
            if differs != marked:
                wrong.append((cells[0], a, b2, "marked" if marked else "unmarked"))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} published rows state two copies of unequal length without "
            f"the word partial, or the reverse: {wrong[:6]}")

    def test_the_partial_pairs_are_still_there_to_be_marked(self):
        """The control: if nothing is partial, the test above passes on an empty set."""
        partial = [p for p in RepeatMap.load(DERIVED).pairs if not p.exact]
        self.assertGreaterEqual(
            len(partial), 5,
            f"only {len(partial)} partial pairs, so the gate above has almost nothing "
            f"to check and the LAT case may have been dropped")
        self.assertTrue(any(p.locus == "LAT" for p in partial),
                        "LAT is the partial pair that matters and it is not in the set")


class EveryAlternativeCarriesItsStrand(unittest.TestCase):
    """Round 9's blocker: the alternatives were on the opposite strand and said nothing.

    `counterparts` computes an inverted-repeat counterpart as an offset from the other
    copy's end, so every one of them is on the opposite strand by construction. Both
    printed forms said only "could equally be 7597", and the README quoted that line as
    flagship output. A user who took the alternative got a position whose allele is the
    complement of theirs: the same class of error as root cause 25, one layer inside the
    fix for root cause 25.

    This checks the claim against the sequence rather than against the code that makes
    it: for each alternative, a 20-mer from the position must equal the 20-mer at the
    alternative (same strand) or its reverse complement (opposite strand), and the flag
    must say which.
    """

    @classmethod
    def setUpClass(cls):
        cls.ref = genbank_sequence(require_reference()).upper()
        cls.m = RepeatMap.load()

    @staticmethod
    def _rc(s: str) -> str:
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def _window(self, position: int, alternative: int, inverted: bool,
                width: int = 20):
        """Compare the two positions over a window clipped to their own repeats.

        Round 10: the unclipped 20-mer criterion was false at 38 of 32,376 alternatives,
        all within 20 bp of a repeat edge, where the window runs past the boundary and
        the flanking sequence legitimately differs. The old test sampled on fixed strides
        and never landed in either range, so it passed by not looking. Clipping the
        window to the region makes the criterion true everywhere, and the edges are now
        sampled at every base.
        """
        home_regions = [r for r in self.m.regions if r.contains(position)]
        alt_regions = [r for r in self.m.regions if r.contains(alternative)]
        if not home_regions or not alt_regions:
            return None
        home_r = max(home_regions, key=lambda r: r.length)
        alt_r = max(alt_regions, key=lambda r: r.length)
        if inverted:
            span = min(width, home_r.end - position + 1,
                       alternative - alt_r.start + 1)
            if span < 1:
                return None
            here = self.ref[position - 1:position - 1 + span]
            there = self._rc(self.ref[alternative - span:alternative])
        else:
            span = min(width, home_r.end - position + 1, alt_r.end - alternative + 1)
            if span < 1:
                return None
            here = self.ref[position - 1:position - 1 + span]
            there = self.ref[alternative - 1:alternative - 1 + span]
        if not here or len(here) != len(there):
            return None
        return here == there, span

    def test_the_strand_flag_matches_the_sequence(self):
        edges = (list(range(9150, 9260)) + list(range(117110, 117220))
                 + list(range(125930, 126420)) + list(range(132540, 132660))
                 + list(range(145560, 145660)) + list(range(1, 450))
                 + list(range(151790, GENOME_BP + 1)))
        positions = sorted(set(list(range(1, GENOME_BP + 1, 7)) + edges))
        wrong, checked, flipped = [], 0, 0
        for position in positions:
            status = self.m.classify(position)
            if status.well_defined:
                continue
            for alternative in status.alternatives:
                claims_inverted = alternative in status.inverted_alternatives
                result = self._window(position, alternative, claims_inverted)
                if result is None:
                    continue
                agrees, _span = result
                checked += 1
                flipped += claims_inverted
                if not agrees:
                    wrong.append((position, alternative,
                                  "opposite" if claims_inverted else "same"))
        self.assertGreater(checked, 5000,
                           f"only {checked} alternatives were checked; the sampling or "
                           f"the repeat map changed")
        self.assertGreater(flipped, 100, "the inverted flag is not being exercised")
        self.assertLess(flipped, checked,
                        "every alternative was claimed inverted, so the flag is a "
                        "constant and the two terminal copies of the 'a' sequence, "
                        "which are direct repeats, are being mislabelled")
        self.assertEqual(wrong, [],
                         f"{len(wrong)} alternatives state the wrong strand "
                         f"(position, alternative, claimed): {wrong[:8]}")

    def test_the_two_terminal_a_sequence_copies_are_called_direct(self):
        """The control on the flag: it must not be True everywhere.

        The 'a' sequence sits at both ends of the linear record in the same orientation
        and once inverted between the segments. Position 200's counterpart at 152,023 is
        a direct repeat; its counterpart at 126,174 is inverted. If the flag were a
        constant, both the biology and this test would be wrong.
        """
        status = self.m.classify(200)
        self.assertIn(152023, status.alternatives)
        self.assertNotIn(152023, status.inverted_alternatives,
                         "the two terminal copies of the 'a' sequence are direct "
                         "repeats of each other, not inverted")
        self.assertIn(126174, status.inverted_alternatives,
                      "the internal copy of the 'a' sequence is inverted relative to "
                      "the terminal copies")

    def test_both_printed_forms_name_the_strand(self):
        """The warning has to reach the text a caller actually reads."""
        status = self.m.classify(118777)
        self.assertIn("opposite strand", str(status))
        result = lift_detailed(118777, REFERENCE, "JQ780693.1", CROSSWALK)
        self.assertIn("opposite strand", str(result),
                      "a lift into an inverted repeat prints its alternative without "
                      "saying the allele there is complemented")
        self.assertEqual(tuple(result.inverted_alternatives),
                         tuple(status.inverted_alternatives))

    def test_a_declined_lift_claims_no_strand_relationship(self):
        """`complemented` described a comparison with a coordinate that does not exist.

        Round 9, minor: a declined lift out of a reverse-complemented genome returned
        position None with complemented True, which reads as a fact about an answer that
        was never given.
        """
        declined = lift_detailed(5, "ON783216.1", REFERENCE, CROSSWALK)
        self.assertIsNone(declined.position)
        self.assertFalse(declined.complemented)
        self.assertNotIn("OPPOSITE", str(declined))

    def test_a_lift_across_strands_still_says_so(self):
        """The control for the test above: the flag must still fire when it should."""
        crossed = lift_detailed(50000, REFERENCE, "ON783216.1", CROSSWALK)
        self.assertIsNotNone(crossed.position)
        self.assertTrue(crossed.complemented)
        self.assertIn("OPPOSITE STRAND", str(crossed))
        same = lift_detailed(50000, REFERENCE, "JQ673480.1", CROSSWALK)
        self.assertFalse(same.complemented)
        self.assertNotIn("OPPOSITE", str(same))


class TheNameColumnIsJoinable(unittest.TestCase):
    """Round 8: three reiteration sets were named `RL1`, `RL2` and `US1`.

    Those are ICP34.5, ICP0 and ICP22, genes of up to 3.6 kb, and the rows carrying
    their names are 53 to 159 bp repeats inside them. Joining this column against the
    gene dictionary returned a 53 bp interval for ICP34.5. Not a coordinate error, but
    a joinability error in the file most likely to be joined.
    """

    @classmethod
    def setUpClass(cls):
        cls.m = RepeatMap.load(DERIVED)
        cls.d = Dictionary.load()

    def test_no_repeat_region_is_named_after_a_gene_it_sits_inside(self):
        wrong = []
        for r in self.m.regions:
            try:
                canonical = self.d.resolve(r.name)
            except Exception:
                continue
            # A region named for a gene is only honest if it covers that gene, not a
            # fragment of it. Anything much shorter is the mislabel this catches.
            if r.length < 500:
                wrong.append((r.name, canonical, r.start, r.end, r.length))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} repeat regions are named after a gene while covering a "
            f"fragment of it, so a join against the dictionary returns the wrong "
            f"interval: {wrong}")

    def test_the_gene_qualifier_is_kept_rather_than_discarded(self):
        """Renaming must not lose the record's own `/gene` value."""
        import csv
        with (DERIVED / "repeat_map.tsv").open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertIn("within_gene", rows[0],
                      "the repeat map has no within_gene column, so the record's gene "
                      "qualifier was dropped rather than moved")
        kept = [r for r in rows if r["within_gene"]]
        self.assertGreaterEqual(
            len(kept), 3,
            f"only {len(kept)} rows carry a gene qualifier; the reference annotates "
            f"three reiteration sets with one, so they have been lost")
        for r in kept:
            with self.subTest(region=r["name"]):
                self.assertNotEqual(r["name"], r["within_gene"])

if __name__ == "__main__":
    unittest.main(verbosity=2)
