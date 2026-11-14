"""The M3 gate: the crosswalk round-trips, and declines rather than guessing."""
from __future__ import annotations

import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import Crosswalk, RepeatMap, lift  # noqa: E402
from hsv1coords.liftover import REFERENCE, SEGMENT_BP  # noqa: E402
from hsv1coords import lift_detailed  # noqa: E402

sys.path.insert(0, str(ROOT / "src"))
from build_crosswalk import genbank_sequence, read_fasta  # noqa: E402
import build_base_agreement  # noqa: E402

CROSSWALK = ROOT / "data" / "derived" / "crosswalk"
BASELINE = ROOT / "data" / "derived" / "base_agreement.tsv"
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


class TheBlocksComeFromTheAlignment(unittest.TestCase):
    """Round 10's blocker: a block can be moved to another matching position.

    The round 9 baseline pins each genome's base-agreement counts, which catches a block
    whose bases stop matching the reference. It cannot catch a block moved to a different
    position where the bases still match, and inside a tandem array the reference is
    periodic, so such positions exist. A reviewer moved 34 blocks in 17 genomes that way:
    3,868 positions lifted to a wrong reference coordinate, by 1 to 48 bp, with all 221
    tests green. One of 36 such moves was caught, by a random 800-window sample, which is
    luck rather than a gate.

    Only an alignment can say where a block belongs, so this runs one and requires the
    stored blocks to be identical to the derived ones, byte for byte, including the
    header fields. `minimap2` is already a declared dependency of this milestone and the
    cost is about 0.08 s per genome.

    The sample rotates with the day so that the whole corpus is covered over time, and
    the three reverse-complemented deposits plus the segmental isomer are always
    included, because they are the cases where the derivation is least obvious.
    """

    ALWAYS = ("ON783216.1", "ON783217.1", "PQ569920.1", "ON783215.1")
    SAMPLE = 12

    @classmethod
    def setUpClass(cls):
        import shutil
        for need in ("ref_NC_001806.2.gb", "genomes.fasta"):
            if not (ROOT / "data" / "raw" / need).exists():
                raise unittest.SkipTest(f"{need} is fetched, not committed")
        if not shutil.which("minimap2"):
            raise unittest.SkipTest(
                "minimap2 is not installed. It is a declared dependency of this "
                "milestone, and without it the blocks cannot be checked against an "
                "alignment, only against themselves.")
        from build_crosswalk import genbank_sequence, read_fasta
        cls.ref = genbank_sequence(ROOT / "data" / "raw" / "ref_NC_001806.2.gb")
        cls.genomes = read_fasta(ROOT / "data" / "raw" / "genomes.fasta")

    def _sample(self, seed: int | None = None) -> list[str]:
        """The genomes to re-align. Deterministic for a given seed.

        The seed is the day, so a failure is reproducible for anyone who runs the suite
        that day, and the corpus is covered over a fortnight rather than the same twelve
        genomes forever. The seed is a parameter so the rotation itself can be tested;
        `datetime.date.today` cannot be monkeypatched, which is how the first version of
        that test failed.
        """
        import datetime
        available = [a for a in Crosswalk.available(CROSSWALK) if self.genomes.get(a)]
        rotating = [a for a in available if a not in self.ALWAYS]
        if seed is None:
            seed = datetime.date.today().toordinal()
        rng = random.Random(seed)
        rng.shuffle(rotating)
        return [a for a in self.ALWAYS if a in available] + rotating[:self.SAMPLE]

    def test_the_stored_blocks_are_the_blocks_an_alignment_gives(self):
        import tempfile
        from build_crosswalk import derive_blocks
        tmp = pathlib.Path(tempfile.mkdtemp())
        ref_fa = tmp / "ref.fa"
        ref_fa.write_text(">ref\n" + self.ref + "\n")
        checked, differences = 0, []
        for acc in self._sample():
            blocks, orientation, reverse_bp, _, _ = derive_blocks(
                acc, self.genomes[acc], tmp / "query.fa", ref_fa)
            cw = Crosswalk.load(acc, CROSSWALK)
            stored = sorted((x.query_start, x.target_start, x.length)
                            for x in cw.blocks)
            checked += 1
            with self.subTest(accession=acc):
                if sorted(blocks) != stored:
                    missing = sorted(set(map(tuple, blocks)) - set(stored))[:3]
                    extra = sorted(set(stored) - set(map(tuple, blocks)))[:3]
                    differences.append((acc, f"{len(blocks)} derived against "
                                             f"{len(stored)} stored",
                                        f"in the alignment only: {missing}",
                                        f"in the file only: {extra}"))
                self.assertEqual(orientation, cw.orientation,
                                 f"{acc}: the file states the wrong frame")
                self.assertEqual(reverse_bp, cw.reverse_strand_bp,
                                 f"{acc}: the file states the wrong dropped amount")
                self.assertEqual(cw.deposited_length, len(self.genomes[acc]),
                                 f"{acc}: the file states the wrong length")
        self.assertEqual(
            differences, [],
            f"{len(differences)} genomes' stored blocks are not what an alignment of "
            f"the deposited sequence gives: {differences}")
        self.assertGreaterEqual(checked, self.SAMPLE,
                                f"only {checked} genomes were re-aligned")
        print(f"\n    re-aligned {checked} genomes; stored blocks identical to the "
              f"alignment in all of them")

    def test_the_sample_rotates(self):
        """A fixed sample would protect twelve genomes and leave 171 open."""
        import datetime
        today = datetime.date.today().toordinal()
        self.assertEqual(self._sample(today), self._sample(today),
                         "the sample is not reproducible for a given day, so a failure "
                         "could not be repeated")
        seen = set()
        for offset in range(14):
            seen.update(self._sample(today + offset))
        self.assertGreater(
            len(seen), self.SAMPLE * 3,
            f"a fortnight of samples covers only {len(seen)} genomes, so the rotation "
            f"is not actually rotating")
        self.assertLessEqual(len(self._sample(today)), self.SAMPLE + len(self.ALWAYS))


class ADroppedSegmentSaysSo(unittest.TestCase):
    """Round 9: `ON783215.1` is deposited as a segmental isomer.

    Its long segment is forward and its short segment reversed. The builder keeps only
    the forward alignments, so about 30 kb of that genome, essentially IRS, US and TRS,
    has no reference coordinate at all. A lift of a variant in US6 (gD) or RS1 (ICP4)
    out of it returned None with the same words as a position in a two-base insertion,
    and nothing in the documents said a whole segment was absent.

    Representing the segment needs per-segment orientation, which is a larger change and
    is listed as an accepted limitation. What is fixed here is the silence.
    """

    def test_the_isomer_is_the_only_genome_with_a_segment_dropped(self):
        """If this stops being one genome, the wording of the decline is wrong."""
        big = {acc: Crosswalk.load(acc, CROSSWALK).reverse_strand_bp
               for acc in Crosswalk.available(CROSSWALK)}
        over = sorted((bp, acc) for acc, bp in big.items() if bp >= SEGMENT_BP)
        self.assertEqual([acc for _, acc in over], ["ON783215.1"],
                         f"genomes with a segment's worth of dropped reverse-strand "
                         f"alignment: {over}")
        others = sorted(bp for acc, bp in big.items() if acc != "ON783215.1")
        self.assertLess(
            others[-1], SEGMENT_BP,
            f"the largest repeat artefact is {others[-1]:,} bp, which is not clearly "
            f"separated from a segment, so the threshold is arbitrary")

    def test_a_position_in_the_dropped_segment_declines_with_the_reason(self):
        result = lift_detailed(140000, "ON783215.1", REFERENCE, CROSSWALK)
        self.assertIsNone(result.position)
        self.assertIn("segmental isomer", result.reason)
        self.assertIn("segmental isomer", str(result),
                      "the reason does not reach the printed form, so nobody is told")

    def test_an_ordinary_decline_in_the_same_genome_is_not_dressed_up(self):
        """The control, in the genome that has the segment.

        Round 10: the first version of this used a different genome, so it could not
        exercise the case at all, which is PITFALLS root cause 12 (verifying a fix with
        inputs that cannot reach it). Inside `ON783215.1` every ordinary decline was
        being told it was in a missing segment: position 1,776 sits in a two-base
        insertion and got the whole paragraph about segmental isomers.
        """
        result = lift_detailed(1776, "ON783215.1", REFERENCE, CROSSWALK)
        self.assertIsNone(result.position)
        self.assertNotIn("segmental isomer", result.reason,
                         "a two-base insertion is being described as a missing segment")
        self.assertIn("2 bp", result.reason,
                      "the decline does not say how large the unmapped stretch is")

    def test_a_decline_in_another_genome_stays_plain(self):
        result = lift_detailed(900, "OP297884.1", REFERENCE, CROSSWALK)
        self.assertIsNone(result.position)
        self.assertNotIn("segmental isomer", result.reason)

    def test_the_unmapped_stretch_is_the_one_the_position_is_in(self):
        """The stretch is measured, so a wrong one would be visible here."""
        cw = Crosswalk.load("ON783215.1", CROSSWALK)
        start, end = cw.unmapped_stretch(140000)
        self.assertLessEqual(start, 140000)
        self.assertGreaterEqual(end, 140000)
        for position in (start, end, 140000):
            with self.subTest(position=position):
                self.assertIsNone(cw.to_reference(position),
                                  f"{position} is inside the stretch and lifts anyway")
        self.assertIsNotNone(cw.to_reference(start - 1),
                             "the stretch starts one base too late")
        self.assertIsNotNone(cw.to_reference(end + 1),
                             "the stretch ends one base too early")
        self.assertIsNone(cw.unmapped_stretch(60000),
                          "a position that lifts has no unmapped stretch")

    def test_every_block_file_states_how_much_it_dropped(self):
        missing = [acc for acc in Crosswalk.available(CROSSWALK)
                   if "reverse_strand_bp="
                   not in (CROSSWALK / f"{acc}.blocks.tsv").read_text()
                   .splitlines()[0]]
        self.assertEqual(missing, [],
                         f"{len(missing)} block files do not state how much alignment "
                         f"was dropped: {missing[:5]}. Rebuild with `make crosswalk`.")


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

    SHORT_BLOCK = build_base_agreement.SHORT_BLOCK

    @classmethod
    def _measured(cls) -> dict[str, dict]:
        """Every genome measured once, by the builder's own function.

        The gate calls `build_base_agreement.measure` rather than reimplementing the
        sampling. Round 7 and round 8 both lost a gate to a reimplementation that
        drifted from the thing it was checking (PITFALLS root causes 13 and 21), and
        this gate's reported coverage is only meaningful if the coverage and the
        scoring come from one definition.
        """
        if not hasattr(cls, "_cache"):
            cls._cache = {acc: build_base_agreement.measure(acc, seq, cls.ref,
                                                            CROSSWALK)
                          for acc in Crosswalk.available(CROSSWALK)
                          if (seq := cls.genomes.get(acc))}
        return cls._cache

    @staticmethod
    def _recorded() -> dict[str, dict]:
        rows = {}
        for line in BASELINE.read_text().splitlines():
            if line.startswith("#") or line.startswith("accession\t"):
                continue
            values = line.split("\t")
            row = dict(zip(build_base_agreement.COLUMNS, values))
            rows[row["accession"]] = {
                k: (v if k in ("accession", "orientation") else int(v))
                for k, v in row.items()}
        return rows

    def test_every_genome_agrees_with_its_recorded_measurement_exactly(self):
        """Equality against `data/derived/base_agreement.tsv`, not a threshold.

        The threshold version of this gate asserted no genome below 95 per cent. The
        corpus runs 95.8 to 99.5, so each genome had a three to four point budget, and
        review round 9 spent it: 1,764,625 bp shifted by one base across 136 genomes,
        all 180 tests green, and the same again at up to 500 bp of displacement. A
        floor beneath what the data achieves measures only how badly it would have to
        break before anyone noticed.

        Exact equality has no budget. One shifted block changes its genome's matching
        count and this names that genome. The cost is that a deliberate rebuild has to
        regenerate the artefact, which is a diff someone can see.
        """
        measured, recorded = self._measured(), self._recorded()
        self.assertEqual(
            sorted(measured), sorted(recorded),
            "the set of genomes with a crosswalk and a sequence is not the set "
            "recorded in base_agreement.tsv; run `make agreement`")
        differences = []
        for acc, row in sorted(measured.items()):
            # Every column, not seven of fourteen. Round 10: six of the recorded
            # columns were written and never read, so `orientation`,
            # `sequence_length`, `sampled_blocks`, `unsampled_bp`, `thin_blocks` and
            # `sampled_positions` could be set to anything in a file whose own header
            # says "pinned exactly".
            for field in build_base_agreement.COLUMNS:
                if field == "accession":
                    continue
                if row[field] != recorded[acc][field]:
                    differences.append(
                        (acc, field, recorded[acc][field], row[field]))
        self.assertEqual(
            differences, [],
            f"{len(differences)} measurements differ from the recorded baseline "
            f"(accession, field, recorded, measured): {differences[:8]}. Either the "
            f"crosswalk changed or a coordinate moved. If the change is intended, "
            f"`make agreement` records it.")
        scored = sum(r["scored"] for r in measured.values())
        matching = sum(r["matching"] for r in measured.values())
        short = sum(r["short_scored"] for r in measured.values())
        print(f"\n    base agreement: {len(measured)} genomes, {scored:,} positions "
              f"scored against the recorded baseline, of which {short:,} in blocks "
              f"under {self.SHORT_BLOCK} bp; pooled {100 * matching / scored:.2f} "
              f"per cent")

    def test_the_worst_genome_is_still_named_and_above_the_floor(self):
        """The floor stays, below the equality gate, because it says something else.

        Equality catches change. A floor catches a corpus that was always wrong, which
        equality against a baseline built from it cannot. Both are needed: the baseline
        is only as good as the day it was written.
        """
        measured = self._measured()
        rates = sorted((100 * r["matching"] / r["scored"], acc)
                       for acc, r in measured.items() if r["scored"] >= 300)
        self.assertEqual(len(rates), len(measured),
                         f"{len(measured) - len(rates)} genomes scored fewer than 300 "
                         f"positions and would have been dropped silently: "
                         f"{[a for a, r in measured.items() if r['scored'] < 300]}")
        poor = [(acc, round(rate, 2)) for rate, acc in rates if rate < 95.0]
        self.assertEqual(poor, [],
                         f"{len(poor)} genomes disagree with the reference at more "
                         f"than 5 per cent of scored positions: {poor}")
        short_poor = [(acc, round(100 * r["short_matching"] / r["short_scored"], 2))
                      for acc, r in sorted(measured.items())
                      if r["short_scored"] >= 200
                      and 100 * r["short_matching"] / r["short_scored"] < 90.0]
        self.assertEqual(short_poor, [],
                         f"{len(short_poor)} genomes disagree inside their short "
                         f"blocks at more than 10 per cent: {short_poor}")
        print(f"    worst genome: {rates[0][1]} at {rates[0][0]:.2f} per cent")

    def test_the_gate_samples_every_block_densely_enough_to_say_anything(self):
        """This gate's coverage, counted by the function that does the scoring.

        It has been wrong twice. The first version set `unscorable = 0`, never
        incremented it, asserted `unscorable == 0` and printed "0 are excluded from
        scoring" whatever the gate did; both round 9 reviewers found that independently.
        The replacement asserted every block contributes at least one sample, which a
        stride of 37 makes true of a 3 bp block, so it was a second tautology in the
        same place and I found it by trying to make it fail.

        What it asserts now is density: a block is sampled at `MIN_SAMPLES` positions
        or at every position. A block sampled once scores 0 or 100 per cent and says
        nothing about the rest of itself, which is the defect the 200 bp cutoff had.
        """
        measured = self._measured()
        mapped = sum(r["mapped_bp"] for r in measured.values())
        unsampled = sum(r["unsampled_bp"] for r in measured.values())
        thin = [(acc, r["thin_blocks"]) for acc, r in sorted(measured.items())
                if r["thin_blocks"]]
        blocks = sum(r["blocks"] for r in measured.values())
        sampled_blocks = sum(r["sampled_blocks"] for r in measured.values())
        positions = sum(r["sampled_positions"] for r in measured.values())
        self.assertEqual(
            unsampled, 0,
            f"{unsampled:,} of {mapped:,} mapped bases sit in blocks the agreement "
            f"gate never samples")
        self.assertEqual(sampled_blocks, blocks,
                         f"{blocks - sampled_blocks} of {blocks} blocks contributed no "
                         f"sampled position")
        self.assertEqual(
            thin, [],
            f"{sum(n for _, n in thin)} blocks are sampled at fewer than "
            f"{build_base_agreement.MIN_SAMPLES} positions without being sampled "
            f"exhaustively, so their agreement figure carries no information about the "
            f"bases between the samples: {thin[:8]}")
        # The spread, because the count cannot see clustering. Round 10: `thin_blocks`
        # is inert at the current parameters, so this check was the third version of
        # itself that could not fail, and a sampling that packed six samples into the
        # first six bases of every block left every printed figure identical.
        clustered = [(acc, r["clustered_blocks"]) for acc, r in sorted(measured.items())
                     if r["clustered_blocks"]]
        worst_gap = max(r["worst_sample_gap"] for r in measured.values())
        self.assertEqual(
            clustered, [],
            f"{sum(n for _, n in clustered)} blocks have a gap between consecutive "
            f"sampled positions larger than "
            f"{build_base_agreement.MAX_SAMPLE_GAP} bp, so part of the block is further "
            f"from any scored base than the sampling claims: {clustered[:8]}")
        self.assertLessEqual(
            worst_gap, build_base_agreement.MAX_SAMPLE_GAP,
            f"the largest gap between sampled positions is {worst_gap} bp")
        print(f"    coverage: {positions:,} positions sampled in {blocks:,} blocks "
              f"over {mapped:,} mapped bases ({100 * positions / mapped:.1f} per cent "
              f"of them), 0 blocks unsampled, 0 sampled too thinly")

    def test_the_short_block_cutoff_cannot_be_widened_without_failing(self):
        """A control on the parameter, because the parameter was unprotected.

        Round 9: setting `SHORT_BLOCK` to 0 reverted round 8's fix, dropped scored
        positions from 1,676,839 to 752,745, and the suite stayed green. The gate that
        reports coverage has to be sensitive to the knob that sets it, so this runs the
        real measurement with the cutoff at 0 and requires the density assertion above
        to have something to catch.
        """
        original = build_base_agreement.SHORT_BLOCK
        # Genomes that actually have short blocks. The first three accessions do not
        # all have one, and a control run on a genome the parameter cannot affect is
        # the "verifying a fix with inputs that cannot exercise it" failure (PITFALLS
        # root cause 12).
        accs = [a for a, r in sorted(self._measured().items())
                if r["short_scored"] > 0][:3]
        self.assertEqual(len(accs), 3, "no genome has a short block to test with")
        try:
            build_base_agreement.SHORT_BLOCK = 0
            mutants = [build_base_agreement.measure(a, self.genomes[a], self.ref,
                                                    CROSSWALK) for a in accs]
        finally:
            build_base_agreement.SHORT_BLOCK = original
        self.assertEqual(build_base_agreement.SHORT_BLOCK, original)
        for acc, mutant in zip(accs, mutants):
            real = self._measured()[acc]
            with self.subTest(accession=acc):
                self.assertGreater(
                    real["scored"], mutant["scored"],
                    "setting the short-block cutoff to 0 changed nothing, so the "
                    "cutoff is not what decides this gate's coverage")
                self.assertGreater(
                    mutant["thin_blocks"] + mutant["clustered_blocks"], 0,
                    "with the cutoff at 0, short blocks get one sample each and must "
                    "be reported as sampled too thinly or too unevenly; they were not, "
                    "so the coverage count is not measuring coverage")

    def test_the_last_base_of_every_block_matches(self):
        """A sharper invariant than the pooled rate, and sensitive to one base.

        Round 9 flagged an asymmetry: pooled over the corpus, the last base of a block
        matches the reference 99.97 per cent of the time and the first base 95.04 per
        cent. That is alignment behaviour, not an offset. A substitution is kept as an
        aligned base, so block boundaries fall only at indels, and the aligner extends
        matches up to a gap while often placing a mismatch just after one.

        It is worth asserting because it is also the signature a one-base frame error
        would produce, and the block-tail figure is far tighter than the 98.58 per cent
        pooled rate: a shift of one base moves it to near zero.
        """
        comp = str.maketrans("ACGT", "TGCA")
        ok = scored = 0
        for acc in Crosswalk.available(CROSSWALK):
            seq = self.genomes.get(acc)
            if not seq:
                continue
            seq = seq.upper()
            cw = Crosswalk.load(acc, CROSSWALK)
            rc = cw.orientation == "reverse_complemented"
            for block in cw.blocks:
                if block.length < 5:
                    continue
                dep = cw._to_deposited(block.query_start + block.length)
                ref_pos = cw.to_reference(dep)
                if ref_pos is None:
                    continue
                if ref_pos > len(self.ref):
                    # Counted as a mismatch, not skipped. Under a one-base shift this
                    # raised IndexError instead of reporting a rate, which still fails
                    # the suite but fails it as a crash rather than as a measurement.
                    scored += 1
                    continue
                a = seq[dep - 1]
                if rc:
                    a = a.translate(comp)
                c = self.ref[ref_pos - 1].upper()
                if a not in "ACGT" or c not in "ACGT":
                    continue
                scored += 1
                ok += a == c
        rate = 100 * ok / scored
        self.assertGreater(scored, 20000,
                           f"only {scored} block tails were scored")
        self.assertGreaterEqual(
            rate, 99.5,
            f"the last base of a block matches at {rate:.3f} per cent over {scored:,} "
            f"blocks. This figure is near 100 by construction, so anything lower means "
            f"the block offsets have moved, not that the assemblies differ.")
        print(f"    block tails: {rate:.3f} per cent match over {scored:,} blocks")

    def test_a_position_before_the_first_block_declines(self):
        """Fable's finding: no test probed a real position before the first block.

        `test_position_before_the_first_block_returns_none` passes position 0, which is
        refused by the `position < 1` guard before any block is consulted, so changing
        the no-block-found branch to return the first block's reference start left the
        suite green while every position before a genome's first block got a confident
        wrong coordinate.
        """
        probed = 0
        for acc in Crosswalk.available(CROSSWALK):
            cw = Crosswalk.load(acc, CROSSWALK)
            if not cw.blocks:
                continue
            first = min(b.query_start for b in cw.blocks)
            if first < 1:
                continue        # nothing lies before it in this genome
            # `query_start` is 0-based, so the first mapped base is `first + 1` in the
            # 1-based deposited frame and `first` is the base before it.
            inside = cw._to_deposited(first + 1)
            before = cw._to_deposited(first)
            probed += 1
            with self.subTest(accession=acc):
                self.assertIsNotNone(
                    cw.to_reference(inside),
                    f"{acc}: the first mapped base does not lift, so this test is not "
                    f"probing the boundary it names")
                self.assertIsNone(
                    cw.to_reference(before),
                    f"{acc}: a position before the first block returned a coordinate")
        self.assertGreater(probed, 20,
                           f"only {probed} genomes have anything before their first "
                           f"block, so this proves little")

    def test_every_block_file_states_the_length_of_its_own_genome(self):
        """`deposited_length` against the sequence, which nothing checked.

        Round 9: changing one forward genome's header from 150,670 to 1,000 left the
        suite green while that genome lost 148,809 of its 149,490 liftable positions.
        It escaped because the round-trip tests cover the first 40 accessions, and the
        agreement gate dropped a genome that fell under 300 scored positions while
        asserting only that 180 of 183 were checked. For the three reverse-complemented
        genomes this number is the pivot of every coordinate in them.
        """
        wrong = []
        for acc in Crosswalk.available(CROSSWALK):
            seq = self.genomes.get(acc)
            if not seq:
                continue
            cw = Crosswalk.load(acc, CROSSWALK)
            # `is not None` made a file that states no length pass a test named for
            # stating it, and `measure` falls back to `len(seq)`, so the baseline missed
            # it too. A missing frame already raises; a missing length must fail here.
            if cw.deposited_length != len(seq):
                wrong.append((acc, cw.deposited_length, len(seq)))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} block files state a deposited_length that is not the length "
            f"of the deposited sequence (accession, stated, actual): {wrong[:8]}")

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
        # Round 7: the threshold was 50 per cent against an observed worst of 3.99, so
        # 182 genomes at 49 per cent N would have passed, and only the single worst
        # genome was checked at all. Both are fixed: every genome is checked, and the
        # threshold sits near what the corpus actually contains. The count printed above
        # is not asserted, because it is a property of the deposited records rather than
        # of this repository, and a new deposit may legitimately change it.
        over = [r for r in worst if r[2] >= 10.0]
        self.assertEqual(over, [],
                         f"{len(over)} genomes are at least 10 per cent ambiguous, "
                         f"which is far outside this corpus and suggests a fetch or "
                         f"parsing fault rather than a deposit: {over}")
        skipped = sum(n for _, n, _ in worst)
        self.assertLess(
            skipped, sum(len(s) for s in self.genomes.values()) * 0.005,
            f"{skipped:,} ambiguous bases across the corpus, more than 0.5 per cent; "
            f"the sequence gates skip these, so at this scale they stop being evidence")


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
