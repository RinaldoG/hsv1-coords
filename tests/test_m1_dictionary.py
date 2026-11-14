"""The M1 gate, as written in ROADMAP.md before the dictionary was built.

Four criteria. Each is a test, and a failure here means the milestone has not been met,
not that the test needs adjusting.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import re
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
        """A name used for two different loci must not silently resolve to one.

        Either refusal counts: `AmbiguousName` when the data holds two canonicals, and
        `UnresolvableName` when it holds one and a row saying it must not be handed out.
        Round 10 added the second case, for ICP10 and AP.
        """
        from hsv1coords.resolver import UnresolvableName
        for entry in self.d.entries:
            if self.d.is_ambiguous(entry.alias):
                with self.subTest(alias=entry.alias):
                    with self.assertRaises((AmbiguousName, UnresolvableName)):
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

class RepeatPairNamingIsNotGuessed(unittest.TestCase):
    """Mutation survivors from review round 7, each with the answer it changed.

    `repeat_pair_of` deciding that two different loci are one repeat pair is the single
    thing `AmbiguousName` exists to prevent, and three mutations to it changed live
    answers with the whole suite green.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load()

    def test_three_repeat_names_are_not_a_pair(self):
        """Mutation: `len(parsed) != 2` to `< 2`, which then raises ValueError.

        A caller gets an unpacking error from inside the function instead of None.
        """
        self.assertIsNone(repeat_pair_of({"TRL1", "IRL1", "TRS1"}))
        self.assertIsNone(repeat_pair_of({"TRL1"}))
        self.assertIsNone(repeat_pair_of(set()))

    def test_two_different_loci_are_never_merged_into_one_pair(self):
        """Mutation: dropping the suffix and segment test returns 'S1' for TRS1/IRS2.

        That is two distinct loci reported as two copies of one, which is the error
        this whole repository exists to prevent.
        """
        self.assertIsNone(repeat_pair_of({"TRS1", "IRS2"}))
        self.assertIsNone(repeat_pair_of({"TRL1", "IRS1"}))
        self.assertIsNone(repeat_pair_of({"TRL1", "TRL1"}))
        self.assertEqual(repeat_pair_of({"TRL1", "IRL1"}), "L1")
        self.assertEqual(repeat_pair_of({"IRS1", "TRS1"}), "S1")

    def test_an_unknown_name_raises_rather_than_returning_none(self):
        """Mutation: `raise KeyError` to `return None` in `resolve`.

        A None that a caller treats as a gene name is worse than an exception, and
        nothing exercised this path.
        """
        with self.assertRaises(KeyError):
            self.d.resolve("not a gene")
        with self.assertRaises(KeyError):
            self.d.resolve("")

class SingleLetterGreekRenderingsResolve(unittest.TestCase):
    """Round 9: the module docstring said `normalise` folds `a-TIF` with `alpha-TIF`.

    It did not. `a-TIF` and `aTIF` both raised `KeyError`, which is a decline rather
    than a wrong answer, so the cost was a dropped record in a text pipeline that
    believed the docstring. Older MEDLINE renderings print a Greek letter as a single
    Latin letter, so those forms do arrive.

    The short keys are added when the dictionary loads, not by `normalise`, because
    folding a leading `a` to `alpha` in general would merge names that are not the same
    name.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load()

    def test_the_single_letter_form_resolves_like_the_spelled_out_form(self):
        pairs = [("alpha-TIF", "a-TIF"), ("alphaTIF", "aTIF")]
        checked = 0
        for long_form, short_form in pairs:
            with self.subTest(alias=short_form):
                self.assertEqual(self.d.resolve(short_form),
                                 self.d.resolve(long_form))
            checked += 1
        self.assertGreater(checked, 1)

    def test_it_adds_nothing_that_was_not_already_implied(self):
        """Every short key must come from a spelled-out alias that is in the table."""
        from hsv1coords.resolver import normalise
        spelled = {normalise(e.alias) for e in self.d.entries}
        words = {"alpha": "a", "beta": "b", "gamma": "g"}
        implied = {letter + key[len(word):]
                   for key in spelled for word, letter in words.items()
                   if key.startswith(word)}
        extra = sorted(k for k in self.d._by_key if k not in spelled
                       and k not in implied)
        self.assertEqual(extra, [],
                         f"the index holds {len(extra)} keys that no cited alias "
                         f"implies: {extra[:8]}")

    def test_it_creates_no_ambiguity_and_overrides_no_cited_alias(self):
        from hsv1coords.resolver import normalise, repeat_pair_of
        for key, targets in sorted(self.d._by_key.items()):
            if len(targets) <= 1 or repeat_pair_of(targets) is not None:
                continue
            with self.subTest(key=key):
                self.assertIn(key, {normalise(e.alias) for e in self.d.entries},
                              f"{key} resolves to {sorted(targets)} and is not a cited "
                              f"alias, so a derived key created an ambiguity")

    def test_a_name_that_merely_starts_with_a_is_untouched(self):
        """The control: this must not make every `a` name a Greek rendering."""
        with self.assertRaises(KeyError):
            self.d.resolve("a-NOTAGENE")


class TheTableDecidesWhatResolves(unittest.TestCase):
    """Round 10: `confidence` was loaded into the dictionary and read by nothing.

    Four rows flagged `ambiguous`, whose notes say in terms that they are not resolved
    here, resolved silently and confidently. ICP10 is the HSV-2 name for the positional
    homologue of UL39, verified against Chung et al. 1989 (PMID 2545912), whose title
    says type 2; it returned UL39 with no signal. AP is attributable to either UL26 or
    the nested UL26.5; it returned UL26. The module docstring says refusing to guess is
    the point, and the table recorded the refusal while the code guessed.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load()

    def test_a_row_marked_decline_does_not_resolve(self):
        from hsv1coords.resolver import UnresolvableName, normalise
        declining = sorted({e.alias for e in self.d.entries
                            if e.resolution == "decline"})
        self.assertGreaterEqual(len(declining), 2,
                                f"only {len(declining)} aliases are marked decline, so "
                                f"this gate is checking almost nothing: {declining}")
        for alias in declining:
            with self.subTest(alias=alias):
                with self.assertRaises(UnresolvableName) as caught:
                    self.d.resolve(alias)
                self.assertTrue(caught.exception.why,
                                f"{alias} declines without saying why")
                self.assertTrue(self.d.is_ambiguous(alias),
                                f"{alias} declines but is_ambiguous says False")
                self.assertTrue(self.d.why_declined(alias))

    def test_every_ambiguous_row_states_what_to_do_with_it(self):
        """The flag has to mean something to the code, not just to a reader."""
        unstated = sorted({e.alias for e in self.d.entries
                           if e.confidence == "ambiguous"
                           and e.resolution not in ("resolve", "decline")})
        self.assertEqual(unstated, [],
                         f"these ambiguous aliases state no resolution: {unstated}")

    def test_icp10_and_ap_specifically(self):
        """The two the reviewer demonstrated, named so a regression is obvious."""
        from hsv1coords.resolver import UnresolvableName
        for alias in ("ICP10", "AP"):
            with self.subTest(alias=alias):
                with self.assertRaises(UnresolvableName):
                    self.d.resolve(alias)

    def test_an_established_alias_still_resolves(self):
        """The control: the flag must not have turned the dictionary off."""
        for alias, canonical in (("ICP0", "RL2"), ("VP16", "UL48"), ("gB", "UL27"),
                                 ("ICP34.5", "RL1")):
            with self.subTest(alias=alias):
                self.assertEqual(self.d.resolve(alias), canonical)
        self.assertEqual(self.d.confidence_of("ICP0"), "established")

    def test_the_confidence_column_is_visible_to_a_caller(self):
        """It was in the data, in the dataclass, and unreachable through the API."""
        self.assertEqual(self.d.confidence_of("ICP10"), "ambiguous")
        self.assertEqual(self.d.confidence_of("nothing-like-this"), "")


class EveryCitationAttestsItsAlias(unittest.TestCase):
    """Round 10: 11 aliases cited a deposited record that does not contain them.

    Three claimed to be the RefSeq `/product` string and were not: the record says
    "ribonucleotide reductase subunit 1", not "large subunit". Those are corrected. The
    other eight are conventional abbreviations such as gB and RR1, which the record does
    not print and whose notes disclose the record's own wording. That is honest, but
    README.md says every alias carries a citation, so the distinction is now a column
    and this checks it.
    """

    RECORDS = {"NC_001806.2": "ref_NC_001806.2.gb",
               "NC_001806.1": "ref_NC_001806.1.gb",
               "BK012101.1": "alt_BK012101.1.gb",
               "JN555585.1": "base_JN555585.1.gb"}

    @classmethod
    def setUpClass(cls):
        cls.rows = [r for r in csv.DictReader(
            (ROOT / "data" / "derived" / "aliases_literature.tsv").open(),
            delimiter="\t")]
        cls.texts = {}
        for ident, filename in cls.RECORDS.items():
            path = ROOT / "data" / "raw" / filename
            if path.exists():
                cls.texts[ident] = path.read_text().lower()
        if not cls.texts:
            raise unittest.SkipTest("the deposited records are fetched, not committed")

    def test_the_verbatim_column_is_the_truth_about_the_record(self):
        wrong = []
        for r in self.rows:
            text = self.texts.get(r["citation_identifier"])
            if text is None:
                if r["verbatim_in_citation"]:
                    wrong.append((r["alias"], "marked for a record not held here"))
                continue
            present = r["alias"].lower() in text
            if (r["verbatim_in_citation"] == "yes") != present:
                wrong.append((r["alias"], r["citation_identifier"],
                              r["verbatim_in_citation"], present))
        self.assertEqual(wrong, [],
                         f"{len(wrong)} rows disagree with the record they cite: "
                         f"{wrong[:8]}")

    def test_a_row_not_verbatim_discloses_the_record_s_own_wording(self):
        """An abbreviation is fine. Claiming to be the record's string is not."""
        bare = []
        for r in self.rows:
            if r["verbatim_in_citation"] != "no":
                continue
            note = r["notes"]
            quoted = re.findall(r"'([^']+)'", note)
            text = self.texts.get(r["citation_identifier"], "")
            if any(q.lower() in text for q in quoted):
                continue
            if "NOT the" in note or "conventional" in note.lower():
                continue
            bare.append((r["canonical_systematic_name"], r["alias"], note[:60]))
        self.assertEqual(
            bare, [],
            f"{len(bare)} aliases are absent from the record they cite and their notes "
            f"quote nothing from it either: {bare[:6]}")

    def test_the_three_product_strings_are_the_record_s(self):
        """Named, because a bare 'RefSeq /product' note was false on exactly these."""
        for gene, product in (
                ("UL39", "ribonucleotide reductase subunit 1"),
                ("UL40", "ribonucleotide reductase subunit 2"),
                ("UL42", "DNA polymerase processivity subunit")):
            with self.subTest(gene=gene):
                rows = [r for r in self.rows
                        if r["canonical_systematic_name"] == gene
                        and r["alias"] == product
                        and r["citation_identifier"] == "NC_001806.2"]
                self.assertTrue(rows, f"{product} is not cited to the record for {gene}")
                self.assertEqual(rows[0]["verbatim_in_citation"], "yes")


if __name__ == "__main__":
    unittest.main(verbosity=2)
