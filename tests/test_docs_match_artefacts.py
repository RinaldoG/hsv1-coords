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

DOCS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")

# The number of derived figures that must stay gated. Raising it is a deliberate edit;
# a change that silently ungates a figure fails instead of printing a smaller number.
COVERAGE_FLOOR = 25
# `\d(?:[\d,]*\d)?` rather than `\d[\d,]*`: the looser form captured a trailing
# comma, so "runs to 152,222," compared "152,222," against "152,222" and reported
# a figure as contradicting itself.
NUMBER = r"(?:\d(?:[\d,]*\d)?(?:\.\d+)?|" + "|".join(
    w for w in ("one two three four five six seven eight nine ten eleven twelve "
                "thirteen fourteen fifteen sixteen seventeen eighteen nineteen "
                "twenty").split()) + r")"


def flat(name: str) -> str:
    return re.sub(r"\s+", " ", (ROOT / name).read_text())


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


class EveryStatedFigureMatchesItsArtefact(unittest.TestCase):
    """The gate. One subtest per (document, phrasing) pair that appears."""

    @classmethod
    def setUpClass(cls):
        cls.facts = derive()
        if not cls.facts:
            raise unittest.SkipTest("no artefacts built; run `make all` first")
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
            forms = [f for f in fact.forms if re.search(r"\d", f)
                     and (("." in f) or len(re.sub(r"\D", "", f)) >= 4)]
            if not forms:
                continue
            pattern = "|".join(re.escape(f) for f in forms)
            for doc, body in self.bodies.items():
                for m in re.finditer(rf"(?<![\d,.]){pattern}(?![\d,.])", body):
                    before = body[max(0, m.start() - 60):m.start()]
                    after = body[m.end():m.end() + 60]
                    # A coordinate range endpoint, `151824..152222` or "9,033 to
                    # 9,213", is a position that happens to equal a figure, not a
                    # statement of that figure. Deliberately narrow: both ends must be
                    # whole numbers of four digits or more, so "from 98.08 to 98.84 per
                    # cent" is still checked, and the repeat intervals themselves are
                    # gated against the reference record in test_m4_repeat_map.py.
                    if before.endswith("..") or after.startswith(".."):
                        continue
                    if ("." not in m.group(0)
                            and re.search(r"(?<![\d,.])\d[\d,]{3,}\s*(?:to|-)\s*$",
                                          before)):
                        continue
                    reached = any(after.lower().startswith(" " + n.lower())
                                  or re.search(rf"{re.escape(n)}\s*$", before, re.I)
                                  for n in fact.nouns + fact.prefixes)
                    if not reached:
                        gaps.append(f"{doc}: '...{before[-40:]}[{m.group(0)}]"
                                    f"{after[:30]}...' states {fact.key} with no "
                                    f"phrasing beside it")
        self.assertEqual(gaps, [],
                         f"{len(gaps)} occurrences of a derived figure are out of the "
                         f"gate's reach:\n  " + "\n  ".join(gaps))

    def test_the_gate_has_power_against_every_document_and_figure(self):
        """Mutate each figure in each document it appears in, and run the real gate.

        Round 6's version had `return` inside two nested loops, so it ran exactly one
        subtest, on one figure, in one document, while being named "against each
        document". It also reimplemented the comparison instead of calling the gate, so
        deleting the gate entirely left the suite green. This calls `contradictions`,
        which is what the gate itself calls.
        """
        tried = 0
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
                    broken = dict(self.bodies)
                    broken[doc] = (body[:m.start(1)] + wrong + body[m.end(1):])
                    tried += 1
                    with self.subTest(document=doc, figure=fact.key, phrasing=label):
                        self.assertTrue(
                            contradictions(self.facts, broken),
                            f"changing {fact.key} to {wrong} in {doc} is not caught")
        self.assertGreater(tried, 25,
                           f"only {tried} figure-in-document pairs were mutated, so "
                           f"this proves little; the phrasings or documents changed")
        print(f"\n    gate power: {tried} figure-in-document mutations, all caught")

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
        print(f"\n    {len(gated)} of {len(self.facts)} derived figures are stated in "
              f"the documents and gated; ungated: {ungated}{context}")
        if absent:
            self.skipTest(f"artefacts missing, so the floor does not apply: {absent}")
        self.assertGreaterEqual(
            len(gated), COVERAGE_FLOOR,
            f"coverage fell to {len(gated)} of {len(self.facts)}, below the recorded "
            f"floor of {COVERAGE_FLOOR}. Ungated: {ungated}")


class ClaimsMatchBehaviour(unittest.TestCase):
    def test_quoted_output_matches_live_output(self):
        from hsv1coords import ConventionSet
        cs = ConventionSet.load()
        _, why = cs.usable(118777, "morse2017_nr")
        interval = re.search(r"inside ([\d,]+-[\d,]+)", why).group(1).rstrip(",")
        for doc in DOCS:
            for raw in re.findall(
                    r"NOT usable under morse2017_nr: inside ([\d,]+-[\d,]+)", flat(doc)):
                with self.subTest(document=doc, quoted=raw):
                    self.assertEqual(raw.rstrip(","), interval,
                                     f"{doc} quotes an interval the code no longer "
                                     f"produces")

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
    assembled to die", inside the test written to prove power. Its 62 of 62 is still
    worth having as a regression check, and it is not a coverage measurement.

    This samples every standalone number in the four documents with a fixed seed,
    mutates one at a time, and runs the real comparison in process. It asserts a floor
    well below the measured rate, because most of those numbers are dates, page
    numbers, DOIs and figures quoted from other groups' records, which this repository
    does not derive and this approach cannot pin. The point is that the rate is
    measured on an independent sample and printed on every run, so a change that
    shrinks the gate shows up as a number rather than as silence.
    """

    SAMPLE = 60
    FLOOR = 10.0
    """A floor on the figure-comparison gate ALONE, not on the suite.

    Measured, this gate catches about 18 per cent of mutations to arbitrary numbers in
    the documents, because most of them are dates, journal page numbers, DOIs and
    counts quoted from other groups' records. A full-suite sweep catches about 50 per
    cent, since the occurrence check, the component-table check, the quoted-output
    checks and the README example runner each catch things this comparison cannot. Both
    figures are in ROADMAP.md, labelled. Quoting the 50 for this gate would flatter it,
    which is how round 7's headline came to overstate it.
    """

    @classmethod
    def setUpClass(cls):
        cls.facts = derive()
        if not cls.facts:
            raise unittest.SkipTest("no artefacts built; run `make all` first")
        cls.bodies = {d: flat(d) for d in DOCS}

    def test_the_measured_catch_rate_is_printed_and_above_its_floor(self):
        import random

        standalone = re.compile(r"(?<![\w.,])(\d(?:[\d,]*\d)?(?:\.\d+)?)(?![\w,.]\d)")
        sites = [(doc, m.start(), m.group(1))
                 for doc, body in self.bodies.items()
                 for m in standalone.finditer(body)]
        self.assertGreater(len(sites), 200,
                           f"only {len(sites)} numbers found, so a sample of "
                           f"{self.SAMPLE} is not independent of the gate")
        random.Random(8).shuffle(sites)
        caught = tried = 0
        for doc, pos, value in sites[:self.SAMPLE]:
            digits = value.replace(",", "")
            try:
                wrong = (f"{float(digits) + 0.7:.2f}" if "." in digits
                         else str(int(digits) + 7))
            except ValueError:          # pragma: no cover
                continue
            body = self.bodies[doc]
            broken = dict(self.bodies)
            broken[doc] = body[:pos] + wrong + body[pos + len(value):]
            tried += 1
            if contradictions(self.facts, broken):
                caught += 1
        rate = 100 * caught / tried
        print(f"\n    measured catch rate, this gate alone: {caught} of {tried} "
              f"independently sampled figure mutations ({rate:.1f} per cent), floor "
              f"{self.FLOOR}")
        self.assertGreaterEqual(
            rate, self.FLOOR,
            f"the gate caught {caught} of {tried} ({rate:.1f} per cent), below the "
            f"recorded floor of {self.FLOOR} per cent")

if __name__ == "__main__":
    unittest.main(verbosity=2)
