"""The M1 gate, as written in ROADMAP.md before the dictionary was built.

Four criteria. Each is a test, and a failure here means the milestone has not been met,
not that the test needs adjusting.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import (AmbiguousName, Dictionary, RepeatPairedName,  # noqa: E402
                        normalise, repeat_pair_of)

TABLE = ROOT / "data" / "derived" / "aliases_literature.tsv"
CASES = ROOT / "tests" / "testcases_synonymy.tsv"
GAPS = ROOT / "tests" / "known_gaps.tsv"


def known_gaps() -> dict[str, dict[str, str]]:
    """Published name strings the dictionary does not yet resolve, each with a reason.

    A gap must be listed here to be tolerated. A new gap fails the suite, and a listed
    gap that has since been filled also fails it, so the list cannot rot.
    """
    with GAPS.open() as fh:
        return {r["name_string"]: r for r in csv.DictReader(fh, delimiter="\t")}


class Criterion1RoundTrip(unittest.TestCase):
    """Name to canonical to name is lossless for every entry."""

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)

    def test_every_alias_resolves(self):
        unresolved = []
        for entry in self.d.entries:
            hits = self.d.resolve_all(entry.alias)
            if entry.canonical not in hits:
                unresolved.append((entry.alias, entry.canonical, hits))
        self.assertEqual(unresolved, [], f"{len(unresolved)} aliases do not resolve back")

    def test_every_canonical_name_is_its_own_alias(self):
        """A canonical name must be resolvable as a name in its own right."""
        missing = [c for c in sorted(self.d.canonical_names)
                   if c not in self.d.resolve_all(c)]
        self.assertEqual(missing, [], f"{len(missing)} canonical names do not self-resolve")

    def test_round_trip_is_lossless(self):
        for canonical in sorted(self.d.canonical_names):
            for alias in self.d.aliases_for(canonical):
                with self.subTest(canonical=canonical, alias=alias):
                    self.assertIn(canonical, self.d.resolve_all(alias))


class Criterion2Citations(unittest.TestCase):
    """No entry carries an uncited alias."""

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)

    def test_no_uncited_row(self):
        uncited = [(e.canonical, e.alias) for e in self.d.entries if not e.citation]
        self.assertEqual(uncited, [], f"{len(uncited)} rows have no citation")

    def test_confidence_is_declared(self):
        allowed = {"established", "probable", "ambiguous"}
        bad = {e.confidence for e in self.d.entries} - allowed
        self.assertEqual(bad, set(), f"undeclared confidence values: {bad}")

    def test_inventory_is_declared(self):
        allowed = {"classical", "whisnant2020", "both"}
        bad = {e.inventory for e in self.d.entries} - allowed
        self.assertEqual(bad, set(), f"undeclared inventory values: {bad}")


class Criterion3Ambiguity(unittest.TestCase):
    """Where a name has been used for more than one gene, it is recorded, not resolved."""

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)

    def test_ambiguous_names_refuse_to_resolve(self):
        """A name used for two different loci must not silently resolve to one."""
        for entry in self.d.entries:
            if self.d.is_ambiguous(entry.alias):
                with self.subTest(alias=entry.alias):
                    with self.assertRaises(AmbiguousName):
                        self.d.resolve(entry.alias)

    def test_ambiguous_names_are_documented_in_the_table(self):
        """A genuine two-locus collision is recorded in the data, not only at runtime.

        The requirement is that the collision is documented, NOT that every row for the
        name is downgraded. `confidence` describes how well sourced one gene/alias pair
        is, which is a different axis from whether the alias collides across genes.
        `IE63` on UL54 is itself well established; what needs recording is that the same
        string has also been applied to US1. Flattening every row to `ambiguous` would
        discard that distinction.
        """
        ambiguous_aliases = {e.alias for e in self.d.entries
                             if self.d.is_ambiguous(e.alias)}
        undocumented = []
        for alias in sorted(ambiguous_aliases):
            rows = [e for e in self.d.entries if e.alias == alias]
            if not any(e.confidence == "ambiguous" for e in rows):
                undocumented.append(alias)
        self.assertEqual(undocumented, [],
                         f"{len(undocumented)} ambiguous names have no row recording "
                         f"the collision: {undocumented}")

    def test_ambiguous_names_carry_an_explanation(self):
        """The row recording a collision says what collides, so a reader can judge it."""
        for e in self.d.entries:
            if e.confidence == "ambiguous" and self.d.is_ambiguous(e.alias):
                with self.subTest(alias=e.alias, canonical=e.canonical):
                    self.assertTrue(e.notes.strip(),
                                    f"{e.alias} on {e.canonical} is marked ambiguous "
                                    f"with no explanation")

    def test_repeat_paired_names_are_not_treated_as_ambiguous(self):
        """One locus in both inverted-repeat copies is not an ambiguity.

        This is the distinction the whole repository exists to make. `TRL2.5` and
        `IRL2.5` are the same locus, present twice because it lies inside the inverted
        repeat, so a name mapping to both is fully determined apart from which copy.
        """
        self.assertIsNotNone(repeat_pair_of({"TRL2.5", "IRL2.5"}))
        self.assertIsNone(repeat_pair_of({"UL54", "US1"}), "different loci, not a pair")
        self.assertIsNone(repeat_pair_of({"TRL2.5", "TRL1"}), "same copy, different loci")
        self.assertIsNone(repeat_pair_of({"TRL2.5", "IRS1"}), "different genome segment")

        paired = [e.alias for e in self.d.entries if self.d.is_repeat_paired(e.alias)]
        for alias in sorted(set(paired)):
            with self.subTest(alias=alias):
                with self.assertRaises(RepeatPairedName):
                    self.d.resolve(alias)
                self.assertFalse(self.d.is_ambiguous(alias))


class Criterion4PublishedSynonymy(unittest.TestCase):
    """Reproduces the synonymy used in at least three published papers."""

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)
        with CASES.open() as fh:
            cls.cases = list(csv.DictReader(fh, delimiter="\t"))
        cls.gaps = known_gaps()

    def test_at_least_three_papers(self):
        papers = {c["pmid"] for c in self.cases if c["pmid"]}
        self.assertGreaterEqual(len(papers), 3, "fewer than three papers in the test set")

    def test_every_published_name_resolves_to_the_right_gene(self):
        failures = []
        for c in self.cases:
            if not c["expected_canonical"]:
                continue  # grouping row, covered by the test below
            if c["name_string"] in self.gaps:
                continue  # listed in known_gaps.tsv with a reason
            hits = self.d.resolve_all(c["name_string"])
            if c["expected_canonical"] not in hits:
                failures.append((c["case"], c["pmid"], c["name_string"],
                                 c["expected_canonical"], hits))
        if failures:
            lines = "\n".join(
                f"      {case} (PMID {pmid}): {name!r} should resolve to "
                f"{expected!r}, got {hits}"
                for case, pmid, name, expected, hits in failures)
            self.fail(f"{len(failures)} of {len(self.cases)} published name strings "
                      f"do not resolve correctly:\n{lines}")

    def test_grouping_row_names_resolve_to_something(self):
        """Names listed under a grouping label must at least be known to the dictionary."""
        unknown = sorted({c["name_string"] for c in self.cases
                          if not c["expected_canonical"]
                          and c["name_string"] not in self.gaps
                          and not self.d.resolve_all(c["name_string"])})
        self.assertEqual(unknown, [], f"{len(unknown)} published names are unknown: {unknown}")

    def test_every_known_gap_has_a_reason_and_a_category(self):
        allowed = {"missing_content", "out_of_scope_m1", "observed_variant"}
        for name, row in sorted(self.gaps.items()):
            with self.subTest(name=name):
                self.assertIn(row["category"], allowed)
                self.assertTrue(row["reason"].strip(), f"{name} has no reason recorded")

    def test_no_known_gap_has_been_silently_filled(self):
        """A gap that now resolves must be removed from the list, not left to rot."""
        stale = []
        for name, row in self.gaps.items():
            hits = self.d.resolve_all(name)
            expected = row["expected_canonical"]
            if hits and (expected == "(any)" or expected in hits):
                stale.append((name, hits))
        self.assertEqual(stale, [],
                         f"{len(stale)} known gaps now resolve and should be "
                         f"removed from known_gaps.tsv: {stale}")

    def test_greek_renderings_fold_together(self):
        """MEDLINE spells Greek letters out; the printed paper uses the character."""
        for a, b in [("alpha TIF", "α-TIF"), ("alphaTIF", "alpha-TIF"),
                     ("ICP0", "icp0")]:
            with self.subTest(pair=(a, b)):
                self.assertEqual(normalise(a), normalise(b))

    def test_normalise_does_not_overreach(self):
        """Names that differ meaningfully must not fold together."""
        for a, b in [("UL20", "UL20.6"), ("UL5", "UL50"), ("US1", "US11"),
                     ("ICP4", "ICP47")]:
            with self.subTest(pair=(a, b)):
                self.assertNotEqual(normalise(a), normalise(b))


if __name__ == "__main__":
    unittest.main(verbosity=2)
