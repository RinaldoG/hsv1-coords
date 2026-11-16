"""The M1 gate, as written in ROADMAP.md before the dictionary was built.

Four criteria. Each is a test, and a failure here means the milestone has not been met,
not that the test needs adjusting.
"""
from __future__ import annotations

import collections
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




def _annotation(record: str) -> str:
    """A record's annotation, lowercased, without its sequence or its translations."""
    return re.sub(r'/translation="[^"]*"', "", record.split("\nORIGIN")[0],
                  flags=re.S).lower()


def _features(record: str) -> list[str]:
    """One lowercased string per feature: its /gene, /locus_tag, /product and /note.

    Round 13 added this. The citation gate asked whether the alias appeared anywhere in
    the record, which attests the STRING and not the mapping, so repointing an alias at
    a different gene left the whole suite green: 138 of 293 (gene, alias) pairs could be
    sent to UL30, among them ICP8, VP22, gG, ICP6 and dUTPase, and deleting `ICP8` from
    the dictionary altogether passed too. What a record attests is that these two names
    sit on the same feature.
    """
    body = record.split("\nORIGIN")[0].split("\nFEATURES", 1)[1]
    blocks, current = [], None
    for line in body.splitlines():
        if re.match(r"^ {5}\S", line):
            if current:
                blocks.append(current)
            current = [line.strip()]
        elif current is not None and line.startswith(" " * 21):
            current.append(line.strip())
    if current:
        blocks.append(current)
    out = []
    for block in blocks:
        text = re.sub(r'/translation="[^"]*"', "", " ".join(block), flags=re.S)
        out.append(" ; ".join(
            re.findall(r'/(?:gene|locus_tag|product|note)="([^"]*)"', text)).lower())
    return out


def _says(name: str, text: str) -> bool:
    """Whether `text` states this name as a word.

    `RL2` is also stated by `TRL2` and `IRL2`: the two repeat copies of one gene, which
    is how BK012101.1 names them (`/gene="IRL2-IRL1"`, `/product="TRL2"`). That is the
    repository's own convention, not a fuzzy match: no other prefix is accepted.
    """
    patterns = [re.escape(name.lower())]
    repeat = re.match(r"^R([LS])(\d+(?:\.\d+)?)$", name, re.I)
    if repeat:
        patterns.append(f"[ti]r{repeat.group(1).lower()}{re.escape(repeat.group(2))}")
    return any(re.search(r"(?<![0-9a-z])" + p + r"(?![0-9a-z])", text)
               for p in patterns)


class WhatTheAnswerSaysBesidesTheName(unittest.TestCase):
    """Round 13: four signals beside `resolve`, each wrong in its own direction.

    A caller who gets a canonical name back needs to know three further things: is
    this gene present twice in the reference, is this string easy to confuse with
    another gene's, and how well attested the mapping is. Each had an accessor and
    each answered something other than the question.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)

    def test_a_two_copy_locus_says_so_whichever_naming_scheme_it_uses(self):
        """ICP0, ICP34.5 and ICP4 are the three loci this repository exists for.

        `resolve("pri-miRNA-H1")` raises `RepeatPairedName`, because the table holds
        both copies under their own names, while these three resolved to RL2, RL1 and
        RS1 with `is_repeat_paired` False and nothing saying the gene is diploid. The
        fact was already in `data/derived/diploid_pairs.tsv`.
        """
        for alias in ("ICP0", "ICP34.5", "ICP4", "IE110", "IE175"):
            with self.subTest(alias=alias):
                self.assertEqual(self.d.copies_in_the_reference(alias), 2)
        for alias in ("VP16", "gB", "UL48"):
            with self.subTest(alias=alias):
                self.assertEqual(self.d.copies_in_the_reference(alias), 1,
                                 "a unique gene is reported as diploid")
        self.assertIsNone(self.d.copies_in_the_reference("nothing-like-this"),
                          "an unknown name must be unknown, not single-copy")
        # And it is not folded into the refusal: these names still resolve.
        self.assertEqual(self.d.resolve("ICP0"), "RL2")

    def test_a_name_the_table_calls_ambiguous_carries_its_caution(self):
        """`is_ambiguous` means `resolve` refuses; the column means easy to confuse.

        Three of the six aliases flagged `ambiguous` resolve perfectly well, among
        them "capsid scaffolding protein" (UL26) and "capsid scaffold protein"
        (UL26.5), one letter apart, whose own note says a string match on either will
        silently pick the wrong gene. A caller branching on `is_ambiguous` was told
        nothing about them.
        """
        for alias in ("capsid scaffolding protein", "capsid scaffold protein"):
            with self.subTest(alias=alias):
                self.assertFalse(self.d.is_ambiguous(alias),
                                 "this name resolves, so the refusal flag is False")
                self.assertTrue(self.d.caution_for(alias),
                                "and nothing warns the caller about it")
        self.assertEqual(self.d.caution_for("VP16"), "",
                         "a name with nothing recorded against it carries no caution")
        flagged = {e.alias for e in self.d.entries if e.confidence == "ambiguous"}
        silent = sorted(a for a in flagged
                        if not self.d.caution_for(a) and not self.d.is_ambiguous(a))
        self.assertEqual(silent, [],
                         f"{len(silent)} aliases are flagged ambiguous in the table "
                         f"and say nothing to a caller: {silent}")

    def test_a_derived_greek_key_carries_what_its_source_key_carries(self):
        """`a-TIF` resolved to UL48 while `confidence_of` said '', which means unknown.

        And `_declined` was built before the short keys existed, so a declined alias
        spelled `alpha-...` was refused under its own spelling and answered under the
        short one.
        """
        self.assertEqual(self.d.resolve("a-TIF"), "UL48")
        self.assertEqual(self.d.confidence_of("a-TIF"),
                         self.d.confidence_of("alpha-TIF"))
        self.assertEqual(self.d.confidences_of("aTIF"),
                         self.d.confidences_of("alpha-TIF"))
        # The refusal travels too. Built here rather than asserted against the
        # committed data, because no alias currently declined begins with a Greek
        # word, which is what made this latent rather than live.
        from hsv1coords.resolver import (Dictionary as _D, Entry, UnresolvableName,
                                         normalise)
        rows = [Entry(canonical="UL48", alias="alpha-SYNTHETIC", alias_type="protein",
                      citation="synthetic", citation_identifier="", inventory="both",
                      confidence="ambiguous", resolution="decline",
                      notes="round 13 fixture: the table refuses this name")]
        built = _D(entries=rows)
        for e in rows:
            built._by_key.setdefault(normalise(e.alias), set()).add(e.canonical)
            if e.resolution == "decline":
                built._declined.setdefault(normalise(e.alias), []).append(e)
        built._add_single_letter_greek_keys()
        for spelling in ("alpha-SYNTHETIC", "a-SYNTHETIC", "aSYNTHETIC"):
            with self.subTest(spelling=spelling):
                with self.assertRaises(UnresolvableName):
                    built.resolve(spelling)

    def test_the_weakest_and_the_strongest_citation_are_both_reachable(self):
        """`confidence_of` returns the weakest, so 77 of 375 names read `probable`.

        A caller filtering on `established` discarded ICP4 and all 68 systematic UL
        and US names, because the 1988 and 1985 naming papers are paywalled and the
        rows say so. The pair is reachable now, and the weakest stays the default.
        """
        self.assertEqual(self.d.confidence_of("ICP4"), "probable")
        self.assertEqual(self.d.confidences_of("ICP4"), ("established", "probable"))
        self.assertIsNone(self.d.confidences_of("nothing-like-this"))
        self.assertEqual(self.d.confidence_of("nothing-like-this"), "")

    def test_a_string_this_dictionary_does_not_hold_names_no_repeat_copy(self):
        """`copy_named("trypsin")` was "terminal", and `copy_kept` said None.

        Any string starting `tr` or `ir` was read as a repeat copy, so the two
        accessors contradicted each other on caller input.
        """
        for word in ("trypsin", "transactivator", "iron", "trimmed"):
            with self.subTest(word=word):
                self.assertIsNone(self.d.copy_named(word))
                self.assertIsNone(self.d.copy_kept(word))
        # The control: the genuine repeat-copy names still answer.
        self.assertEqual(self.d.copy_named("IRL2"), "internal")
        self.assertEqual(self.d.copy_named("TRL2"), "terminal")
        self.assertFalse(self.d.copy_kept("IRL2"),
                         "IRL2 resolves to RL2, so the copy is dropped")



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
        cls.texts, cls.features = {}, {}
        for ident, filename in cls.RECORDS.items():
            path = ROOT / "data" / "raw" / filename
            if path.exists():
                record = path.read_text()
                cls.texts[ident] = _annotation(record)
                cls.features[ident] = _features(record)
        if not cls.texts:
            # A skip here is a hole: these three tests are the only gate tying an
            # alias to the record it cites, and `data/raw` is gitignored, so on a
            # clone they vanished while the suite reported OK and the discovery count
            # stayed at 38. Round 13 added
            # `test_a_checkout_holds_every_artefact_the_gates_skip_without` in
            # tests/test_m5_conventions.py, which fails in a checkout before this can
            # skip quietly.
            raise unittest.SkipTest("the deposited records are fetched, not committed")

    def test_the_verbatim_column_is_the_truth_about_the_record(self):
        wrong = []
        for r in self.rows:
            text = self.texts.get(r["citation_identifier"])
            if text is None:
                if r["verbatim_in_citation"]:
                    wrong.append((r["alias"], "marked for a record not held here"))
                continue
            # A word-boundary match, written here independently of the builder's.
            # Round 11: both used a bare substring test, so 15 rows were marked attested
            # because the record says "glycoprotein L" and "polymerase", and the gate
            # shared the defect it was meant to catch.
            #
            # Round 13: both then searched the whole FILE. The genome is 152,222 bp,
            # which is 60 x 2537 + 2, so `ORIGIN`'s last line is "   152221 gc", the
            # only occurrence of `gC` in either record, and `gI` occurs only as the
            # last two residues of a UL50 `/translation`. Three rows were certified
            # against the sequence. The sequence and the translations are not what the
            # depositors wrote, so neither is searched here.
            present = _says(r["alias"], text)
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
            # Round 13: there was a blanket exemption here for a note containing
            # "NOT the" or "conventional". Measured over the 26 `verbatim=no` rows, all
            # 26 passed on the quoted-string branch above and none used it, so it
            # admitted nothing and stood ready to wave through any row whose note
            # happened to use the word. Standing rule 46: an exemption inside a gate is
            # named and counted, and an uncounted one is deleted.
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


class EveryMappingIsAttestedAndTheUngatedOnesAreCounted(unittest.TestCase):
    """Round 13: nothing gated the (gene, alias) PAIR, only the alias string.

    `Criterion1RoundTrip` cannot fail by construction: `Dictionary.load` puts
    `canonical` into `_by_key[normalise(alias)]` for every row, so
    `entry.canonical in resolve_all(entry.alias)` is true whatever the row says.
    `EveryCitationAttestsItsAlias` asks whether the record prints the alias and never
    looks at the gene beside it. `testcases_synonymy.tsv` pins 74 of 366 pairs. So a
    reviewer repointed 138 of 293 pairs at UL30 with all 285 tests green, among them
    ICP8, VP22, gG, ICP6, dUTPase, OBP, UDG and every ICTV functional name, and
    deleting `ICP8` from the dictionary altogether also passed.

    What can be checked offline is checked, per source class, and what cannot is
    counted rather than left unsaid (standing rule 16):

    - a GenBank-cited row claiming the record prints its alias: the alias and the gene
      have to sit on the SAME feature. 205 of 205 repoints to another gene fail this.
    - a row cited to the Zenodo GTF: the alias has to be an identifier in the pinned
      file, and the gene has to be that same identifier or an identifier itself.
    - UniProt, PubMed, DOI and ISBN rows: not checkable without the network, which the
      suite does not use. Counted and printed, not assumed.
    """

    @classmethod
    def setUpClass(cls):
        cls.rows = list(csv.DictReader(
            (ROOT / "data" / "derived" / "aliases_literature.tsv").open(),
            delimiter="\t"))
        cls.features = {}
        for ident, filename in (("NC_001806.2", "ref_NC_001806.2.gb"),
                                ("BK012101.1", "alt_BK012101.1.gb")):
            path = ROOT / "data" / "raw" / filename
            if path.exists():
                cls.features[ident] = _features(path.read_text())
        gtf = ROOT / "data" / "raw" / "HSV1_Annotation.gtf"
        cls.gtf_ids = set()
        if gtf.exists():
            for line in gtf.read_text().splitlines():
                for field in ("gene_id", "transcript_id"):
                    found = re.search(rf'{field} "([^"]+)"', line)
                    if found:
                        cls.gtf_ids.add(found.group(1))
        if not cls.features or not cls.gtf_ids:
            raise unittest.SkipTest("the cited records are fetched, not committed")

    def _classes(self):
        """Each row with the class of source it cites. Every row lands in one."""
        for row in self.rows:
            ident = row["citation_identifier"]
            if ident in self.features:
                yield "GenBank record", row
            elif ident.startswith("10.5281/zenodo"):
                yield "Zenodo GTF", row
            elif ident.startswith("PMID:"):
                yield "PubMed", row
            elif ident.startswith("ISBN:"):
                yield "book", row
            elif ident.startswith("10."):
                yield "DOI", row
            else:
                yield "UniProt", row

    def _genbank_failures(self):
        out = []
        for kind, row in self._classes():
            if kind != "GenBank record" or row["verbatim_in_citation"] != "yes":
                continue
            gene = row["canonical_systematic_name"]
            if not any(_says(gene, f) and _says(row["alias"], f)
                       for f in self.features[row["citation_identifier"]]):
                out.append((gene, row["alias"], row["citation_identifier"]))
        return out

    def test_a_record_cited_for_an_alias_carries_it_on_that_gene_s_feature(self):
        failures = self._genbank_failures()
        self.assertEqual(
            failures, [],
            f"{len(failures)} rows cite a record that does not state the alias on any "
            f"feature naming that gene, so the record attests the string and not the "
            f"mapping: {failures}")

    def test_repointing_any_gated_row_at_another_gene_fails_that_gate(self):
        """The control, measured rather than assumed (standing rule 18).

        A gate that passes on the committed data says nothing until something it
        should reject is put through it.
        """
        caught = tried = 0
        for kind, row in self._classes():
            if kind != "GenBank record" or row["verbatim_in_citation"] != "yes":
                continue
            tried += 1
            elsewhere = "UL30" if row["canonical_systematic_name"] != "UL30" else "UL29"
            if not any(_says(elsewhere, f) and _says(row["alias"], f)
                       for f in self.features[row["citation_identifier"]]):
                caught += 1
        gtf_caught = gtf_tried = 0
        for kind, row in self._classes():
            if kind != "Zenodo GTF":
                continue
            gtf_tried += 1
            alias = row["alias"]
            elsewhere = "UL30" if row["canonical_systematic_name"] != "UL30" else "UL29"
            stem = alias
            for suffix in ("_gene", "_RNA"):
                stem = stem[:-len(suffix)] if stem.endswith(suffix) else stem
            if elsewhere not in (alias, stem):
                gtf_caught += 1
        self.assertGreater(tried, 0, "no GenBank-cited row was put through the gate")
        self.assertEqual(
            caught, tried,
            f"{tried - caught} of {tried} rows would still pass with their gene "
            f"replaced by another, so the gate is not reading the pair")
        self.assertEqual(gtf_caught, gtf_tried,
                         f"{gtf_tried - gtf_caught} GTF-cited rows would still pass "
                         f"repointed at another gene")
        # Measured over the whole table by replaying the reviewer's sweep, 16 November
        # 2026: repointing each of the 378 distinct (gene, alias) pairs in turn and
        # running this file, 320 fail and 58 pass. Every survivor is cited only to a
        # paper, a book or UniProt, which this suite cannot fetch. Before round 13 the
        # same sweep caught about half.
        print(f"\n    mapping gate: {caught} of {tried} GenBank-cited rows and "
              f"{gtf_caught} of {gtf_tried} GTF-cited rows fail when the gene is "
              f"replaced")

    def test_every_gtf_cited_alias_is_an_identifier_in_the_pinned_file(self):
        failures = []
        for kind, row in self._classes():
            if kind != "Zenodo GTF":
                continue
            gene, alias = row["canonical_systematic_name"], row["alias"]
            # The GTF carries `gene_id` and `transcript_id` and no synonym field, so
            # what it can attest is that the alias IS the gene, written as that file
            # writes it: `LAT_gene` for `LAT` is the one row where the two differ.
            # Round 13: requiring only that each be an identifier somewhere in the
            # file left a repointed row passing, because the gene it was repointed to
            # is an identifier too.
            stem = alias
            for suffix in ("_gene", "_RNA"):
                stem = stem[:-len(suffix)] if stem.endswith(suffix) else stem
            if alias not in self.gtf_ids:
                failures.append((gene, alias, "the alias is not an identifier"))
            elif gene not in (alias, stem):
                failures.append((gene, alias,
                                 f"the file states {alias!r}, which is not {gene!r} "
                                 f"and carries no synonym field that could say it is"))
        self.assertEqual(
            failures, [],
            f"{len(failures)} rows cite the Zenodo GTF for a name it does not carry: "
            f"{failures}")

    def test_the_share_of_mappings_a_gate_can_reach_is_measured_and_printed(self):
        """Printed per class, and asserted to cover every row.

        The number that matters is not how many rows pass; it is how many rows any
        gate could fail. Round 12 recorded that as the lesson of the documentation
        gate and the alias table was never measured the same way.
        """
        counted = collections.Counter(kind for kind, _ in self._classes())
        self.assertEqual(sum(counted.values()), len(self.rows),
                         "a row fell into no source class, so it is counted nowhere")
        checkable = {"GenBank record", "Zenodo GTF"}
        reachable = sum(n for kind, n in counted.items() if kind in checkable)
        # A GenBank row that does not claim the record prints its alias is out of
        # reach too: there is nothing in the record to compare against.
        unclaimed = sum(1 for kind, row in self._classes()
                        if kind == "GenBank record"
                        and row["verbatim_in_citation"] != "yes")
        print(f"\n    alias mappings a gate reads: {reachable - unclaimed} of "
              f"{len(self.rows)} rows, by class {dict(counted)}, with {unclaimed} "
              f"GenBank rows whose alias the record does not print and "
              f"{len(self.rows) - reachable} rows cited to a paper, a book or UniProt, "
              f"which this suite cannot fetch.")
        self.assertGreater(reachable - unclaimed, 0)



class TheDeclaredGapsAreRealGaps(unittest.TestCase):
    """A gap declared in `known_gaps.tsv` must still be a gap, or the file is stale.

    Round 10: `oriL` did not resolve and was not declared, unlike the miRNAs, while
    `oriS` resolved to the name of a reading frame rather than the origin. Both are
    declared now, and this checks the declaration against the dictionary and the repeat
    map so that it cannot quietly stop being true.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)
        cls.gaps = known_gaps()

    def test_the_origins_are_declared_and_their_coordinates_are_findable(self):
        from hsv1coords import RepeatMap
        regions = {r.name for r in RepeatMap.load().regions}
        # oriL does not resolve, so it is a declared gap. oriS does resolve, to the
        # name of a reading frame, so it is a documented collision instead and is
        # covered by the test below.
        self.assertIn("oriL", self.gaps, "oriL is not declared in known_gaps.tsv")
        self.assertIn("repeat_map.tsv", self.gaps["oriL"]["reason"])
        for feature in ("oriL palindrome", "oriS palindrome"):
            with self.subTest(feature=feature):
                self.assertIn(feature, regions,
                              f"{feature} is not in the repeat map, so the reason given "
                              f"in known_gaps.tsv points at nothing")

    def test_the_oris_collision_behaves_as_the_documents_say(self):
        """`oriS` resolves, so it is not a declared gap; it is a documented collision.

        A hyphen is a separator in the normaliser, so `oriS` and `Ori-S` fold to one
        key, and the answer is the 2020 re-annotation's reading frame spanning the TRS
        copy and TRS1, cited to Hubenthal-Voss 1987, not the oriS palindrome. Round 10
        found that behaviour undocumented. Writing it into `known_gaps.tsv` would have
        been wrong, since that file lists names which do NOT resolve, so it is stated in
        the module docstring and asserted here.
        """
        from hsv1coords import RepeatMap
        self.assertEqual(self.d.resolve("oriS"), "Ori-S")
        self.assertEqual(self.d.resolve("Ori-S"), "Ori-S")
        regions = {r.name: (r.start, r.end) for r in RepeatMap.load().regions}
        self.assertIn("oriS palindrome", regions,
                      "the palindrome is not in the repeat map, so the distinction "
                      "this test draws has nowhere to point")
        docstring = __import__("hsv1coords.resolver",
                               fromlist=["x"]).__doc__ or ""
        self.assertIn("oriS", docstring,
                      "the collision is not documented where a reader would meet it")
        self.assertIn("not the oriS palindrome", docstring)


class WhetherTheRepeatCopySurvives(unittest.TestCase):
    """Round 10: the canonical set mixes two schemes and the answer did not say which.

    `resolve("IRL2")` gives `RL2`, dropping the copy the caller named; `resolve("IRL2.5")`
    gives `IRL2.5`, keeping it. Both are correct, because the records differ: the
    reference annotates RL2 once and the 2020 re-annotation annotates TRL2 and IRL2
    separately, with copy-specific decimal names that have no copy-agnostic form. What
    was wrong is that a caller could not tell the two cases apart, on the distinction
    this repository exists to keep.

    This class was written in round 10 and lost before it was committed: an edit that
    rewrote this file as a slice of itself silently deleted everything after the cut,
    including the `__main__` block, and round 11 found the commit touching only the
    source file. See PITFALLS root cause 40.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Dictionary.load(TABLE)

    def test_a_copy_specific_input_reports_whether_the_copy_survived(self):
        cases = {"IRL2": False, "TRL2": False, "IRS1": False, "TRS1": False,
                 "IRL2.5": True, "TRL1.5": True, "IRL3": True, "TRL3": True}
        for name, kept in sorted(cases.items()):
            with self.subTest(name=name):
                self.assertIsNotNone(
                    self.d.copy_named(name),
                    f"{name} names a copy and copy_named says it does not")
                self.assertEqual(
                    self.d.copy_kept(name), kept,
                    f"{name} resolves to {self.d.resolve(name)}; copy_kept should be "
                    f"{kept}")

    def test_a_name_with_no_copy_reports_none(self):
        for name in ("ICP0", "VP16", "UL48", "LAT"):
            with self.subTest(name=name):
                self.assertIsNone(self.d.copy_named(name))
                self.assertIsNone(self.d.copy_kept(name))

    def test_both_outcomes_occur_so_the_report_is_not_a_constant(self):
        outcomes = set()
        for e in self.d.entries:
            if self.d.copy_named(e.alias) is not None:
                outcomes.add(self.d.copy_kept(e.alias))
        self.assertEqual(outcomes, {True, False},
                         f"copy_kept takes only {outcomes}, so it carries no "
                         f"information")

    def test_the_mixture_is_documented_where_a_reader_meets_it(self):
        source = (ROOT / "src" / "hsv1coords" / "resolver.py").read_text()
        self.assertIn("mixes two schemes", source,
                      "the inconsistency is not explained in the module that has it")


if __name__ == "__main__":
    unittest.main(verbosity=2)
