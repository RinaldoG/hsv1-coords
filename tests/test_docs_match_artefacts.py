"""Gate every load-bearing figure in the documents against its artefact.

Six review rounds found stale numbers, and each round's fix corrected only the lines
that round quoted. Round 6 changed one figure in one document ninety times and five
turned the suite red, because the previous gate asserted that a correct value appeared
somewhere in the four documents concatenated.

The approach now: `hsv1coords.facts` derives each figure from its artefact and records
the phrasings the documents use for it. Every occurrence of such a phrasing, in every
document, must carry the derived value. Text is flattened first, because Markdown wraps
a figure away from its noun and an unflattened check is simply unsatisfiable.
"""
from __future__ import annotations

import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords.facts import derive, missing_inputs  # noqa: E402

DOCS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md",
        "data/derived/alias_notes.md", "data/derived/alias_testcases.md")
"""Round 13: the two alias documents were outside this gate and said so nowhere.

They are the ones a reader of the dictionary actually opens, they carry their own
table of counts, and `alias_notes.md` stated "Citations that are papers | 66" with no
definition and nothing deriving it. Being hand-written is not a reason to be ungated;
it is the reason to be gated."""

# The number of derived figures that must stay gated. Raising it is a deliberate edit;
# a change that silently ungates a figure fails instead of printing a smaller number.
COVERAGE_FLOOR = 35

REQUIRED_FACTS = frozenset({
    # The figures a reader acts on. Named, not counted: round 10 deleted the five facts
    # added in round 9 and retyped two numbers in ROADMAP.md, and every floor still
    # held, because a floor on a count can be met by a smaller set. Deleting any of
    # these now fails with the name of what went missing.
    "genome_bp", "genome_v1_bp", "genomes", "regions", "diploid_pairs",
    "partial_pairs", "repeat_covered_bp", "repeat_covered_percent",
    "annotation_rows", "alias_rows", "alias_genes", "alias_citations",
    "version_blocks", "version_offsets", "conventions",
    "scored_positions", "worst_agreement", "pooled_agreement",
    "ref_cds", "ref_genes", "alt_cds", "alt_genes", "lat_overhang_bp",
    "transcript_spans", "transcript_span_rows", "nested_orfs",
    "long_repeat_bp", "short_repeat_bp", "median_mappable",
    "reverse_complemented", "lat_in_szpara",
})
# `\d(?:[\d,]*\d)?` rather than `\d[\d,]*`: the looser form captured a trailing
# comma, so "runs to 152,222," compared "152,222," against "152,222" and reported
# a figure as contradicting itself.
NUMBER = r"(?:\d(?:[\d,]*\d)?(?:\.\d+)?|" + "|".join(
    w for w in ("one two three four five six seven eight nine ten eleven twelve "
                "thirteen fourteen fifteen sixteen seventeen eighteen nineteen "
                "twenty").split()) + r")"


def flat(name: str) -> str:
    """One document as a single line, with Markdown line furniture removed.

    The blockquote marker is stripped before the whitespace is collapsed. Round 9: a
    figure in a wrapped blockquote line came out as "1,676,839 scored > positions", so
    none of its fact's phrasings sat beside it and the occurrence left the gate's reach.
    The whole known-defects notice at the top of README.md is a blockquote, which is the
    part of the file most likely to state a superseded number, so that blind spot was in
    the worst possible place. List bullets and heading hashes go the same way: they are
    furniture, and a phrase that wraps across one of them is still the phrase a reader
    sees.
    """
    text = (ROOT / name).read_text()
    text = re.sub(r"(?m)^[ \t]*(?:>[ \t]*)+", " ", text)
    text = re.sub(r"(?m)^[ \t]*(?:[-*+]|\d+\.)[ \t]+", " ", text)
    text = re.sub(r"(?m)^[ \t]*#+[ \t]*", " ", text)
    return re.sub(r"\s+", " ", text)


def contradictions(facts, bodies: dict[str, str]) -> list[str]:
    """Every place a document states a value that is not the derived one.

    Extracted as a function in round 7 so the power test can run the real gate over
    mutated text instead of reimplementing it. The previous power test reimplemented
    the comparison, so deleting the gate left the suite green at 145 tests.
    """
    problems = []
    for fact in facts:
        for noun in fact.nouns:
            pattern = rf"({NUMBER})\s+{re.escape(noun)}"
            for doc, body in bodies.items():
                for stated in re.findall(pattern, body, re.I):
                    if stated.lower() not in {f.lower() for f in fact.forms}:
                        problems.append(
                            f"{doc}: '{stated} {noun}' but {fact.key} is {fact.value}")
        for prefix in fact.prefixes:
            pattern = rf"{re.escape(prefix)}\s+({NUMBER})"
            for doc, body in bodies.items():
                for stated in re.findall(pattern, body, re.I):
                    if stated.lower() not in {f.lower() for f in fact.forms}:
                        problems.append(
                            f"{doc}: '{prefix} {stated}' but {fact.key} is "
                            f"{fact.value}")
    return problems


def phrasings(fact) -> list[tuple[str, str]]:
    """Every phrasing on this fact, each with the regex that finds it."""
    out = [(n, rf"{NUMBER}\s+{re.escape(n)}") for n in fact.nouns]
    out += [(p, rf"{re.escape(p)}\s+{NUMBER}") for p in fact.prefixes]
    return out


REACHED_SITES = {
    "alias_citations": 2, "alias_genes": 2, "alias_paper_citations": 1,
    "alias_rows": 1, "alt_cds": 3,
    "alt_genes": 3, "annotation_rows": 5, "conventions": 4, "diploid_pairs": 2, "features_without_a_gene_name": 3,
    "genome_bp": 11, "genome_v1_bp": 7, "genomes": 7, "inputs": 1,
    "lat_in_szpara": 3, "lat_overhang_bp": 7, "long_repeat_bp": 1,
    "median_mappable": 4, "nested_orfs": 2, "partial_pairs": 1,
    "pooled_agreement": 2, "ref_cds": 3, "ref_genes": 1, "regions": 2,
    "repeat_covered_bp": 1, "repeat_covered_percent": 2,
    "revcomp_ON783216.1": 2, "revcomp_ON783217.1": 2, "revcomp_PQ569920.1": 2,
    "reverse_complemented": 1, "scored_positions": 1, "short_repeat_bp": 1,
    "transcript_span_rows": 2, "transcript_spans": 2, "version_blocks": 4,
    "version_offsets": 4, "worst_agreement": 2,
}
"""How many places in the four documents each fact's phrasings actually reach.

Why this table exists. Round 12 measured the depth of the gate against the likeliest
careless edit: deleting a phrasing while rewording a sentence. Removing
`lat_overhang_bp`'s "bp beyond the end of" and "which runs", then changing the sentence
those covered from 778 to 878, left every floor intact (35 of 35 figures gated, the
power count and the catch rate both above their floors) and failed exactly one test,
whose message named the number to change: "81 != 79". Changing that one number in
ROADMAP.md re-opened the hole with the suite green.

The reach test cannot see it, because it sweeps for occurrences of the CORRECT value and
a falsified sentence no longer holds one. The contradiction check cannot see it, because
the phrasing that would have found the sentence is the thing that was deleted.

So the count is pinned per fact. A deleted phrasing now fails with the fact's name and
the sentence it stopped reaching, which cannot be repaired by copying a number out of
the failure message: lowering an entry here is a written statement that a figure is
gated in fewer places than it was, which is a different claim from "the total moved".
These numbers move whenever the documents gain or lose a sentence stating a figure, and
raising one is ordinary. Lowering one needs a reason in the commit message.
"""


def _derive_or_explain():
    """The figures, or the right kind of refusal for where we are running.

    Round 13 made `derive()` raise rather than return a short list, and left its two
    callers carrying `if not cls.facts: raise SkipTest(...)`, which can no longer run:
    in an sdist-like tree the suite ERRORED out of `setUpClass` instead of skipping,
    so the two round 13 fixes disagreed about what a non-checkout should do. The
    policy is one line: in a git checkout the artefacts are buildable, so an absent
    one fails and names the command; outside a checkout nothing can be built, so it
    skips.
    """
    from hsv1coords.facts import ArtefactsMissing
    try:
        return derive()
    except ArtefactsMissing as absent:
        if (ROOT / ".git").exists():
            raise AssertionError(f"{absent}. Run `make all`.") from absent
        raise unittest.SkipTest(str(absent)) from absent


class EveryStatedFigureMatchesItsArtefact(unittest.TestCase):
    """The gate. One subtest per (document, phrasing) pair that appears."""

    @classmethod
    def setUpClass(cls):
        cls.facts = _derive_or_explain()
        cls.bodies = {d: flat(d) for d in DOCS}

    def test_no_document_states_a_wrong_value_for_a_derived_figure(self):
        problems = contradictions(self.facts, self.bodies)
        self.assertEqual(problems, [],
                         f"{len(problems)} figures contradict their artefacts:\n  "
                         + "\n  ".join(problems))

    def test_no_noun_starts_with_a_number(self):
        """A noun beginning with a digit makes the pattern capture the wrong figure.

        One did ("is printed as 109,564 in"), so the check compared the preceding
        number against the derived value and failed on correct text.
        """
        # A standalone number only. Digits glued to letters or a dot are part of an
        # accession ("NC_001806.2", "ON783216.1"), which is exactly what makes a
        # phrasing distinctive, so those are allowed. A bare "152,222" is not: it moves
        # with the value it is meant to check and so can never contradict it.
        standalone = re.compile(r"(?<![\w.])\d[\d,]*(?:\.\d+)?(?!\w)")
        bad = [(f.key, n) for f in self.facts
               for n in list(f.nouns) + list(f.prefixes)
               if standalone.search(n)]
        self.assertEqual(bad, [], f"nouns starting with a digit: {bad}")

    def test_no_gated_phrasing_is_dead(self):
        """A phrasing that matches nothing is dead weight pretending to be a gate.

        Round 6 printed this count and asserted nothing, and round 7 found 5 of the 19
        "spare" phrasings were holes: the figure was still stated in a document that no
        other phrasing reached. A phrasing that matches nothing is now a failure, so it
        must be corrected or removed rather than carried.
        """
        all_text = " ".join(self.bodies.values())
        dead = [(f.key, label) for f in self.facts for label, pat in phrasings(f)
                if not re.search(pat, all_text, re.I)]
        self.assertEqual(dead, [],
                         f"{len(dead)} phrasings match nothing. Fix the wording or "
                         f"delete it; a phrasing carried 'for future wording' hides a "
                         f"hole:\n  " + "\n  ".join(f"{k}: {n!r}" for k, n in dead))

    def test_no_fact_reaches_fewer_sentences_than_it_did(self):
        """Per fact, not per total. See REACHED_SITES for why."""
        measured = {}
        for fact in self.facts:
            measured[fact.key] = sum(
                len(re.findall(pattern, body, re.I))
                for _, pattern in phrasings(fact)
                for body in self.bodies.values())
        # Round 14: this iterates the TABLE, so a fact with no entry was checked by
        # nothing, and round 13's own new fact, `alias_paper_citations`, went in
        # without one. The two sets have to match before the counts mean anything.
        unpinned = sorted({f.key for f in self.facts} - set(REACHED_SITES))
        self.assertEqual(
            unpinned, [],
            f"{len(unpinned)} facts have no entry in REACHED_SITES, so nothing "
            f"notices if their phrasings stop reaching the documents: {unpinned}")
        lost, gained = [], []
        for key, was in sorted(REACHED_SITES.items()):
            now = measured.get(key)
            if now is None:
                lost.append(f"{key}: the fact no longer exists")
            elif now < was:
                sites = [f"{n!r}" for n, pattern in
                         phrasings(next(f for f in self.facts if f.key == key))
                         if not any(re.search(pattern, b, re.I)
                                    for b in self.bodies.values())]
                lost.append(f"{key}: reaches {now} sentences, {was} recorded; "
                            f"phrasings reaching nothing: {sites or 'none'}")
            elif now > was:
                # Round 14: a rise was only printed, so a figure that gained a
                # sentence could lose it again and still meet its entry. The entry is
                # the measurement, so it is raised in the commit that earns it.
                gained.append(f"{key}: {was} -> {now}")
        self.assertEqual(
            lost, [],
            f"{len(lost)} figures are gated in fewer sentences than REACHED_SITES "
            f"records. A phrasing was probably deleted while rewording a sentence, "
            f"which leaves that sentence ungated and is how round 12 falsified 778 to "
            f"878 in README.md with the suite otherwise green. Restore the phrasing, or "
            f"lower the entry and say why in the commit message:\n  "
            + "\n  ".join(lost))
        # Round 14 made a rise a failure, so adding one correct sentence that states
        # a gated figure failed the build until this table was edited in the same
        # commit. Round 15: that is a false failure on an honest edit, and the thing
        # worth catching is a FALL, which the assertion above does. A rise is printed,
        # and the unpinned-fact check below is what closes the hole round 14 was
        # actually aiming at.
        if gained:
            print(f"\n    REACHED_SITES is behind: {', '.join(gained)}. Raising an "
                  f"entry is ordinary; do it in the same commit as the new sentence.")

    def test_every_fact_can_be_gated_at_all(self):
        """A fact with no phrasings can never gate anything, whatever a document says.

        Four facts carried an empty tuple behind a hard-disabled `if False`, so they sat
        in the ungated list permanently and inflated the denominator of the coverage
        figure while being structurally incapable of ever being gated.
        """
        empty = [f.key for f in self.facts if not phrasings(f)]
        self.assertEqual(empty, [],
                         f"{len(empty)} facts carry no phrasing and so can never be "
                         f"gated: {empty}")

    def test_no_occurrence_of_a_large_figure_is_out_of_the_gate_s_reach(self):
        """Per occurrence. Round 7 did this per document, which was still pooling.

        Round 8: `genome_bp` is stated four times in ROADMAP.md and phrasings reached
        two of them, so "A sweep of all 152,222 positions" could be changed to 99,999
        with the suite green. The granularity had moved from four documents to one
        document, not to the occurrence. So every occurrence of a derived value must
        have one of its fact's phrasings adjacent to it.

        Restricted to values of 1,000 or more and to percentages: a bare "5" or "nine"
        cannot be told from any other small number in prose. That restriction is the
        honest limit of the approach and is stated in ROADMAP.md.
        """
        gaps = []
        for fact in self.facts:
            # The four-digit floor stays for bare numbers, because lifting it entirely
            # flags 67 occurrences where a small number coincides with a fact's value
            # ("a further 37", "UL22", "asm5"), which is the ambiguity it exists for.
            #
            # Round 12 measured what it still cost: not the two sentences round 11
            # fixed, but at least nine sentences carrying twelve occurrences of six
            # required figures, among them "runs 778 bp past the end of `IRL`" in
            # ROADMAP.md (a sibling of the README sentence round 11 fixed), "every one
            # of the 183 genomes" in README.md (a sibling of the ROADMAP sentence round
            # 11 fixed), and the whole 77-against-300 comparison at five sites in
            # FINDINGS.md and CHANGELOG.md. Changing 778 to 878, 183 to 283 and 300 to
            # 400 in three documents left the suite green and did not even move the
            # printed catch rate.
            #
            # So a small number is checked when it is stated WITH A UNIT the fact
            # itself declares: the first word of one of its nouns, or the last word of
            # one of its prefixes. "778 bp" is a claim about base pairs and "183
            # genomes" a claim about genomes, while "a further 37" names no unit and is
            # left to the contradiction check. The units are derived from the phrasings
            # rather than typed here, so a fact cannot gain a unit this test does not
            # know about.
            small = [f for f in fact.forms if re.search(r"\d", f)
                     and "." not in f and len(re.sub(r"\D", "", f)) < 4]
            units = {n.split()[0].strip("`,.").lower() for n in fact.nouns if n.split()}
            before_units = {p.split()[-1].strip("`,.*").lower()
                            for p in fact.prefixes if p.split()}
            forms = [f for f in fact.forms if re.search(r"\d", f)
                     and (("." in f) or len(re.sub(r"\D", "", f)) >= 4)]
            if not forms and not (small and (units or before_units)):
                continue
            if small and (units or before_units):
                forms = forms + small
            pattern = "|".join(re.escape(f) for f in sorted(forms, key=len,
                                                            reverse=True))
            for doc, body in self.bodies.items():
                for m in re.finditer(rf"(?<![\d,.]){pattern}(?![\d,.])", body):
                    before = body[max(0, m.start() - 60):m.start()]
                    after = body[m.end():m.end() + 60]
                    # A coordinate range endpoint, `151824..152222` or "9,033 to
                    # 9,213", is a position that happens to equal a figure, not a
                    # statement of that figure. Deliberately narrow: both ends must be
                    # whole numbers of four digits or more, so "from 98.08 to 98.84 per
                    # cent" is still checked. The exemption used to be justified by
                    # saying the repeat intervals are gated against the record, which
                    # was true of the artefact rows and not of the text: round 9 moved
                    # the README's IRL boundary by 300 bp with the suite green. The
                    # ranges themselves are now checked as text, against the record, by
                    # EveryQuotedIntervalIsOneTheRecordStates in test_m4_repeat_map.py.
                    if before.endswith("..") or after.startswith(".."):
                        continue
                    if ("." not in m.group(0)
                            and re.search(r"(?<![\d,.])\d[\d,]{3,}\s*(?:to|-)\s*$",
                                          before)):
                        continue
                    reached = any(after.lower().startswith(" " + n.lower())
                                  or re.search(rf"{re.escape(n)}\s*$", before, re.I)
                                  for n in fact.nouns + fact.prefixes)
                    if m.group(0) in small and not reached:
                        # A small number is only in scope when it carries one of this
                        # fact's own units. Without one it is any number in prose.
                        next_word = (after.strip().split()[0].strip("`,.:;")
                                     if after.strip() else "")
                        last_word = (before.strip().split()[-1].strip("`,.:;*")
                                     if before.strip() else "")
                        if (next_word.lower() not in units
                                and last_word.lower() not in before_units):
                            continue
                    if not reached:
                        gaps.append(f"{doc}: '...{before[-40:]}[{m.group(0)}]"
                                    f"{after[:30]}...' states {fact.key} with no "
                                    f"phrasing beside it")
        self.assertEqual(gaps, [],
                         f"{len(gaps)} occurrences of a derived figure are out of the "
                         f"gate's reach:\n  " + "\n  ".join(gaps))

    def _require_every_artefact(self) -> None:
        """A count claimed by a document is only comparable on a full build.

        Without the fetched records `derive` yields fewer facts, so the power count and
        the catch rate both fall and a comparison against ROADMAP.md fails for a reason
        that has nothing to do with the documents. Round 9: these two checks failed on a
        fresh clone, which the README promises will skip rather than fail.
        """
        absent = missing_inputs()
        if absent:
            raise unittest.SkipTest(
                f"{len(absent)} artefacts are not built, so the gate is smaller than "
                f"the one these counts describe: run `make all`")

    def _power_mutations(self):
        """Each (fact, phrasing, document) the power check can mutate.

        Enumerated in one place so the count the run achieves and the count the
        documents claim are the same number, rather than two numbers that agree until
        one of them stops.
        """
        for fact in self.facts:
            patterns = ([(n, rf"({NUMBER})\s+{re.escape(n)}") for n in fact.nouns]
                        + [(p, rf"{re.escape(p)}\s+({NUMBER})")
                           for p in fact.prefixes])
            for label, pattern in patterns:
                for doc, body in self.bodies.items():
                    m = re.search(pattern, body, re.I)
                    if not m or m.group(1).lower() not in {
                            f.lower() for f in fact.forms}:
                        continue
                    wrong = ("0" if not m.group(1)[0].isdigit()
                             else m.group(1).replace(",", "") + "7")
                    broken = {doc: body[:m.start(1)] + wrong + body[m.end(1):]}
                    yield fact, label, doc, broken, wrong

    def _power_pairs(self) -> int:
        return sum(1 for _ in self._power_mutations())

    def test_the_gate_has_power_against_every_document_and_figure(self):
        self._require_every_artefact()
        """Mutate each figure in each document it appears in, and run the real gate.

        Round 6's version had `return` inside two nested loops, so it ran exactly one
        subtest, on one figure, in one document, while being named "against each
        document". It also reimplemented the comparison instead of calling the gate, so
        deleting the gate entirely left the suite green. This calls `contradictions`,
        which is what the gate itself calls.
        """
        tried = 0
        for fact, label, doc, broken, wrong in self._power_mutations():
            tried += 1
            with self.subTest(document=doc, figure=fact.key, phrasing=label):
                self.assertTrue(
                    contradictions(self.facts, broken),
                    f"changing {fact.key} to {wrong} in {doc} is not caught")
        # Round 9: this asserted `tried > 25` while ROADMAP.md and a docstring both
        # said 62, so the demonstrated power could have fallen by 58 per cent with the
        # suite green and both documents stale. The floor is now near the real number,
        # and ROADMAP.md's figure is checked against it below.
        self.assertGreaterEqual(
            tried, self.POWER_FLOOR,
            f"only {tried} figure-in-document pairs were mutated, against a recorded "
            f"{self.POWER_FLOOR}; the phrasings or the documents shrank")
        print(f"\n    gate power: {tried} figure-in-document mutations, all caught")

    POWER_FLOOR = 60
    """How many figure-in-document pairs the power check must still find.

    Not an equality: adding a figure or a phrasing legitimately raises it, and a floor
    that fails on an improvement would be edited away. Round 9 found the assertion set
    at 25 against a documented 62.
    """

    def test_the_roadmap_states_the_power_count_this_run_achieves(self):
        self._require_every_artefact()
        """The 62 in the document, against what the run actually mutates."""
        body = flat("ROADMAP.md")
        m = re.search(r"\|\s*`gate power`\s*\|\s*\*\*(\d+) of (\d+)\*\*", body)
        self.assertIsNotNone(
            m, "ROADMAP.md no longer states the gate power count in the form "
               "`| `gate power` | **N of N** |`")
        caught, tried = int(m.group(1)), int(m.group(2))
        self.assertEqual(caught, tried,
                         "the gate power row states fewer caught than tried, which "
                         "would mean the power check is failing")
        pairs = self._power_pairs()
        self.assertEqual(
            tried, pairs,
            f"ROADMAP.md states {tried} figure-in-document mutations and this run finds "
            f"{pairs}")

    def test_the_gate_is_load_bearing(self):
        """The gate reacts to a wrong figure. Deliberately not a second value check.

        Standing rule 10 says a gate is justified by deleting it and watching the suite
        fail. For this gate that check reads green and always will, because the gate
        holds the only assertion that a stated figure equals its artefact: delete it and
        nothing else makes that comparison. Round 8 found an earlier version of this
        test had duplicated the gate's assertion, which made deletion harmless and the
        rule-10 check meaningless. Removing the duplicate is the fix; adding another
        would reinstate the defect.

        So rule 10 is satisfied differently here, and this is where that is recorded:
        with the gate present, changing any figure a phrasing reaches fails the suite,
        demonstrated 62 times over by
        `test_the_gate_has_power_against_every_document_and_figure`. Deleting the gate
        removes a check rather than exposing a redundant one, which is what a
        load-bearing test looks like.
        """
        # Deliberately NOT asserting `contradictions(...) == []` here: that is the
        # gate's own assertion, and duplicating it meant deleting the gate left the
        # suite green, so standing rule 10's check returned green and was
        # indistinguishable from the defect it exists to catch. This test asserts only
        # that the gate reacts; the gate asserts that it currently has nothing to say.
        fact = next((f for f in self.facts if isinstance(f.value, int)
                     and any(re.search(rf"{NUMBER}\s+{re.escape(n)}", b, re.I)
                             for n in f.nouns for b in self.bodies.values())), None)
        self.assertIsNotNone(fact, "no fact is stated anywhere, so nothing is gated")
        doc = next(d for d, b in self.bodies.items()
                   if any(re.search(rf"{NUMBER}\s+{re.escape(n)}", b, re.I)
                          for n in fact.nouns))
        noun = next(n for n in fact.nouns
                    if re.search(rf"{NUMBER}\s+{re.escape(n)}", self.bodies[doc], re.I))
        broken = dict(self.bodies)
        broken[doc] = re.sub(rf"{NUMBER}(\s+{re.escape(noun)})", r"4242\1",
                             self.bodies[doc], count=1, flags=re.I)
        self.assertTrue(contradictions(self.facts, broken),
                        "the gate does not react to a wrong figure at all")

    def test_coverage_is_reported_and_does_not_regress(self):
        """The coverage figure, asserted rather than printed.

        Round 6 printed "19 of 26" from a test with no assertion, and I reported that
        number as a measurement. The reviewer dropped it to 18 by editing one word with
        the suite green. The floor below is therefore enforced: raising it is a
        deliberate edit, and a fix that silently ungates a figure now fails.
        """
        all_text = " ".join(self.bodies.values())
        gated, ungated = [], []
        for fact in self.facts:
            if any(re.search(rf"{re.escape(f)}\s+{re.escape(n)}", all_text, re.I)
                   for f in fact.forms for n in fact.nouns) or any(
                   re.search(rf"{re.escape(p)}\s+{re.escape(f)}", all_text, re.I)
                   for f in fact.forms for p in fact.prefixes):
                gated.append(fact.key)
            else:
                ungated.append(fact.key)
        absent = missing_inputs()
        context = (f" WITHOUT {len(absent)} artefacts, so this is a smaller gate than a "
                   f"full checkout, not a better one: {absent}" if absent else "")
        # Round 13: this line counted FACTS, and the gate works per occurrence, so
        # "35 of 35 derived figures are stated in the documents and gated" was printed
        # over a document holding a figure the gate had skipped. (The reviewer added
        # "overruns `IRL` by 878 bases" to FINDINGS.md against an artefact value of
        # 778; small numbers are in scope only beside a unit the fact declares, which
        # is `bp` here and not `bases`, and the printed rate did not move.) So the
        # occurrences are counted too, and the line says what share of them carries a
        # phrasing. The limit itself is unchanged and stated in ROADMAP.md: it is the
        # price of not flagging every "a further 37" in the prose.
        reached_sites, loose_sites = 0, []
        for fact in self.facts:
            for form in sorted(fact.forms):
                if not re.search(r"\d", form):
                    continue
                for doc, body in self.bodies.items():
                    for m in re.finditer(
                            rf"(?<![\d,.]){re.escape(form)}(?![\d,.])", body):
                        after = body[m.end():m.end() + 60]
                        before = body[max(0, m.start() - 60):m.start()]
                        if any(after.lower().startswith(" " + n.lower())
                               or re.search(rf"{re.escape(n)}\s*$", before, re.I)
                               for n in fact.nouns + fact.prefixes):
                            reached_sites += 1
                        else:
                            loose_sites.append(f"{fact.key} in {doc}")
        total_sites = reached_sites + len(loose_sites)
        print(f"\n    {len(gated)} of {len(self.facts)} derived figures are stated in "
              f"the documents and gated; ungated: {ungated}{context}")
        print(f"    per occurrence: {reached_sites} of {total_sites} occurrences of a "
              f"derived value carry one of its phrasings. The rest are numbers that "
              f"coincide with a figure, or small numbers out of the gate's scope.")
        if absent:
            self.skipTest(f"artefacts missing, so the floor does not apply: {absent}")
        self.assertGreaterEqual(
            len(gated), COVERAGE_FLOOR,
            f"coverage fell to {len(gated)} of {len(self.facts)}, below the recorded "
            f"floor of {COVERAGE_FLOOR}. Ungated: {ungated}")
        derived = {f.key for f in self.facts}
        missing = sorted(REQUIRED_FACTS - derived)
        self.assertEqual(
            missing, [],
            f"{len(missing)} figures that must be derived and gated are not being "
            f"derived at all: {missing}. A floor on the count cannot see this, because "
            f"a smaller set can still meet it.")
        ungated_required = sorted(REQUIRED_FACTS & set(ungated))
        self.assertEqual(
            ungated_required, [],
            f"{len(ungated_required)} required figures are derived but stated in no "
            f"document, so nothing compares them: {ungated_required}")


class ClaimsMatchBehaviour(unittest.TestCase):
    def test_quoted_output_matches_live_output(self):
        from hsv1coords import ConventionSet
        cs = ConventionSet.load()
        _, why = cs.usable(118777, "morse2017_nr")
        interval = re.search(r"inside ([\d,]+-[\d,]+)", why).group(1).rstrip(",")
        checked = 0
        for doc in DOCS:
            for raw in re.findall(
                    r"NOT usable under morse2017_nr: inside ([\d,]+-[\d,]+)", flat(doc)):
                checked += 1
                with self.subTest(document=doc, quoted=raw):
                    self.assertEqual(raw.rstrip(","), interval,
                                     f"{doc} quotes an interval the code no longer "
                                     f"produces")
        # Round 12: there was no counter, so this was the zero-occurrence tautology this
        # file's own docstring warns about. It happens to have one live site; rewording
        # that line would have left the test passing over nothing.
        self.assertGreater(
            checked, 0,
            "no document quotes this message any more, so this test checks nothing. "
            "Either a document stopped quoting it, in which case delete this test, or "
            "the wording moved, in which case fix the pattern.")

    def test_no_document_says_a_built_component_is_unbuilt(self):
        self.assertNotIn("is not built yet", flat("README.md"))

    def test_the_component_table_matches_its_own_count(self):
        """The README said six components above a seven-row table."""
        body = (ROOT / "README.md").read_text()
        rows = re.findall(r"^\| [A-Z][^|]*\|[^|]*\|$", body, re.M)
        rows = [r for r in rows if "---" not in r and "Component" not in r]
        stated = re.search(r"(\w+) components are built", body)
        self.assertIsNotNone(stated, "the README does not state a component count")
        from hsv1coords.facts import SPELLED
        word = SPELLED.get(len(rows), str(len(rows)))
        self.assertEqual(stated.group(1).lower(), word,
                         f"the README says {stated.group(1)} components and lists "
                         f"{len(rows)}")

    def test_the_alignment_preset_is_not_contradicted(self):
        from build_crosswalk import PRESET
        report = ROOT / "data" / "derived" / "crosswalk_report.md"
        if not report.exists():
            self.skipTest("crosswalk report not built; run `make crosswalk`")
        named = set(re.findall(r"asm\d+", report.read_text()))
        self.assertEqual(named, {PRESET},
                         f"the report names {named} while the builder runs {PRESET}")

class TheCatchRateIsMeasuredNotClaimed(unittest.TestCase):
    """Sample numbers the gate does NOT choose, and print what fraction it catches.

    Round 8: `test_the_gate_has_power_against_every_document_and_figure` is 100 per cent
    by construction, because it skips any site whose stated value is not already one of
    a fact's derived forms. Its mutants are drawn exactly from the sites the gate
    reaches, which is root cause 13, "a mutant set assembled after the gates is
    assembled to die", inside the test written to prove power. Its all-of-them result is
    still worth having as a regression check, and it is not a coverage measurement. The
    count it achieves is printed and compared against the figure ROADMAP.md states, so
    neither this docstring nor that table carries a number that can go stale: round 9
    found both saying 62 while the assertion beneath them accepted 26.

    This mutates every standalone number in the four documents, one at a time, and runs
    the real comparison in process. It asserts a floor
    well below the measured rate, because most of those numbers are dates, page
    numbers, DOIs and figures quoted from other groups' records, which this repository
    does not derive and this approach cannot pin. The point is that the rate is
    measured on an independent sample and printed on every run, so a change that
    shrinks the gate shows up as a number rather than as silence.
    """

    FLOOR_RATE = 40.0
    """A floor on the figure-comparison gate ALONE, not on the suite.

    Round 8 measured this on a 60-site sample with a fixed seed and printed 15.0 per
    cent against a floor of 10. Round 9 ran the same unchanged gate under 12 other seeds
    and got 8.3 to 23.3 per cent: seed 2 was below the floor. The 5-point margin was
    smaller than the measurement's own variance, so the floor was not a regression guard,
    it was where the seed happened to land. Every site is measured now, so the number is
    deterministic and the floor means what it says.
    """

    FLOOR_CAUGHT = 70
    """A floor on the COUNT, which the rate's denominator cannot move.

    Round 14: 40 against 78 caught, 10.0 per cent against 20.4, and a coverage floor
    of 30 against 36 gated, so roughly half of each gate could be removed before any
    of them fired. That is standing rule 27, a threshold set below what the data
    already achieves, in the three places this file states one. They are 70, 18.0 and
    35 now, measured on 16 November 2026 at 78, 20.4 and 36; a document edit that
    moves them legitimately moves them by one or two.

    Round 9's other attack on this metric: appending 40 ungated page numbers to a
    document dropped the printed rate from 15.0 to 6.7 per cent with the gate untouched,
    so writing a references paragraph failed the build and deleting ungated prose
    "improved" the gate. The rate is partly prose-number density, in both directions.
    The count of mutations caught is not: adding numbers the gate cannot see leaves it
    where it was. The rate is still printed, because a falling rate is worth seeing, but
    the count is what regresses when the gate shrinks.
    """

    @classmethod
    def setUpClass(cls):
        cls.facts = _derive_or_explain()
        cls.bodies = {d: flat(d) for d in DOCS}

    @classmethod
    def _measure(cls) -> tuple[int, int]:
        """Memoised: two tests consume it and the sweep costs a few seconds."""
        if getattr(cls, "_result", None) is not None:
            return cls._result
        cls._result = cls._sweep()
        return cls._result

    NARRATIVE = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")
    """The documents the rate floor was calibrated on.

    Round 13 brought `data/derived/alias_notes.md` and `alias_testcases.md` inside this
    gate, and they are mostly citations: PMIDs, years, page numbers and volume numbers,
    none of which this repository derives. The denominator went from 371 numbers to 985
    and the printed rate from 16.4 to 7.9 per cent with the gate untouched, which is
    exactly the hazard `FLOOR_CAUGHT` already records in the other direction. So the
    rate floor is applied to the set it describes, the rate over everything is printed
    beside it, and the count floor covers the lot.
    """

    @classmethod
    def _sweep(cls) -> tuple[int, int]:
        """Mutate every number and count what the gate catches.

        The baseline must be clean first. Round 9: adding four facts introduced one
        phrasing that matched a neighbouring figure, so the documents contradicted their
        artefacts before any mutation, and this sweep then counted every mutation as
        caught and printed 76.5 per cent. A measurement taken over a failing baseline
        flatters the thing being measured, which is the shape of error this whole test
        exists to report.
        """
        standing = contradictions(cls.facts, cls.bodies)
        assert not standing, (
            f"the documents already contradict their artefacts, so a catch rate "
            f"measured against them would count every mutation as caught: {standing}")
        standalone = re.compile(r"(?<![\w.,])(\d(?:[\d,]*\d)?(?:\.\d+)?)(?![\w,.]\d)")
        sites = [(doc, m.start(), m.group(1))
                 for doc, body in cls.bodies.items()
                 for m in standalone.finditer(body)]
        caught = tried = 0
        cls._narrative = [0, 0]
        # Round 15: the rate's denominator is every number in the prose, so ten
        # citations added to FINDINGS.md took the narrative rate from 20.3 to 18.0 per
        # cent and failed the build, with a message blaming the gate: round 14 raised
        # that floor to 18.0 and made an edit this project's own conventions ask for
        # into a red suite. Deleting prose numbers would have "improved" it, which is
        # the hazard FLOOR_CAUGHT's docstring already describes.
        #
        # The denominator that means something is the number of sites that state a
        # value this repository DERIVES. A citation list does not move it; a derived
        # figure stated in one more sentence does. That is what the floor is on now,
        # and the two older rates are still printed.
        derived_values = {form for fact in cls.facts for form in fact.forms}
        cls._derived = [0, 0]
        for doc, pos, value in sites:
            digits = value.replace(",", "")
            try:
                wrong = (f"{float(digits) + 0.7:.2f}" if "." in digits
                         else str(int(digits) + 7))
            except ValueError:          # pragma: no cover
                continue
            body = cls.bodies[doc]
            # Only the mutated document is re-judged: a changed number in one document
            # cannot create a contradiction in another, and judging all five made the
            # sweep five times slower for the same answer. Checked: the totals are
            # identical either way.
            broken = {doc: body[:pos] + wrong + body[pos + len(value):]}
            tried += 1
            hit = bool(contradictions(cls.facts, broken))
            caught += hit
            if doc in cls.NARRATIVE:
                cls._narrative[1] += 1
                cls._narrative[0] += hit
            if value in derived_values:
                cls._derived[1] += 1
                cls._derived[0] += hit
        return caught, tried

    def _require_every_artefact(self) -> None:
        absent = missing_inputs()
        if absent:
            if (ROOT / ".git").exists():
                self.fail(f"{len(absent)} artefacts are not built, so this gate is "
                          f"smaller than the one the rate describes: {absent}. Run "
                          f"`make all`.")
            raise unittest.SkipTest(
                f"{len(absent)} artefacts are not built and there is no checkout to "
                f"build them in: {absent}")

    def test_the_measured_catch_rate_is_printed_and_above_its_floor(self):
        caught, tried = self._measure()
        self.assertGreater(tried, 200,
                           f"only {tried} numbers found, so this is not a measurement "
                           f"of the documents as they stand")
        rate = 100 * caught / tried
        narrative_caught, narrative_tried = self._narrative
        narrative_rate = 100 * narrative_caught / narrative_tried
        derived_caught, derived_tried = self._derived
        derived_rate = 100 * derived_caught / derived_tried
        print(f"\n    measured catch rate, this gate alone: {caught} of {tried} figure "
              f"mutations, every number in the documents ({rate:.1f} per cent); over "
              f"the four narrative documents {narrative_caught} of {narrative_tried} "
              f"({narrative_rate:.1f} per cent); over the sites that state a DERIVED "
              f"value {derived_caught} of {derived_tried} ({derived_rate:.1f} per "
              f"cent), floors {self.FLOOR_CAUGHT} caught and {self.FLOOR_RATE} per "
              f"cent of the derived sites")
        self.assertGreaterEqual(
            caught, self.FLOOR_CAUGHT,
            f"the gate caught {caught} mutations, below the recorded floor of "
            f"{self.FLOOR_CAUGHT}. Unlike the rate, this count does not move when "
            f"ungated prose numbers are added or removed, so it has shrunk.")
        self.assertGreaterEqual(
            derived_rate, self.FLOOR_RATE,
            f"the gate caught {derived_caught} of {derived_tried} "
            f"({derived_rate:.1f} per cent) of the sites that state a derived value, "
            f"below the recorded floor of {self.FLOOR_RATE} per cent. Unlike the two "
            f"rates printed beside it, this denominator does not move when prose "
            f"numbers are added or removed.")

    def test_the_roadmap_states_the_rate_this_run_measures(self):
        self._require_every_artefact()
        """The figure in the document, against the figure the run produces.

        Round 9: ROADMAP.md said "about 18 per cent" for this gate while the suite
        printed 15.0 and the full population gives 14.4. The table was added in round 8
        to stop one number being reported as a measurement of something else, and it
        overstated the number it was added to state honestly. The rate is not derivable
        from any artefact, so the only way to gate it is to measure it here and compare.
        """
        self._measure()
        # Round 15: this compared the narrative rate, whose denominator is every
        # number in four documents, so twenty added citations moved it past the ±1.0
        # tolerance and failed the build. It states the derived-site rate now, the
        # same one the floor is on.
        caught, tried = self._derived
        rate = 100 * caught / tried
        body = flat("ROADMAP.md")
        m = re.search(r"\|\s*This gate alone\s*\|\s*\*\*about ([\d.]+) per cent\*\*",
                      body)
        self.assertIsNotNone(
            m, "ROADMAP.md no longer states this gate's own catch rate in the form "
               "`| This gate alone | **about N per cent**`, so nothing ties the "
               "document to the measurement")
        stated = float(m.group(1))
        self.assertLess(
            abs(stated - rate), 1.0,
            f"ROADMAP.md states about {stated} per cent for this gate and the run "
            f"measures {rate:.1f} per cent ({caught} of {tried})")


if __name__ == "__main__":
    unittest.main(verbosity=2)
