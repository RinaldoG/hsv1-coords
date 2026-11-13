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
NUMBER = r"(?:\d[\d,]*(?:\.\d+)?|" + "|".join(
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

    def test_no_large_figure_is_stated_in_a_document_the_gate_cannot_see(self):
        """Per document, not pooled. This is the hole round 7 measured.

        A figure gated in one document was counted as gated everywhere, because the
        coverage check searched the four documents concatenated. `genome_bp` was stated
        five times and reachable at one of them. So: wherever a document states a
        derived value, some phrasing for that fact must match in THAT document.

        Restricted to values of 1,000 or more and to percentages, because a bare "5" or
        "nine" cannot be told from any other small number in prose. That restriction is
        the honest limit of this approach and is stated in ROADMAP.md.
        """
        gaps = []
        for fact in self.facts:
            forms = [f for f in fact.forms if re.search(r"\d", f)
                     and (("." in f) or len(re.sub(r"\D", "", f)) >= 4)]
            if not forms:
                continue
            for doc, body in self.bodies.items():
                present = any(re.search(rf"(?<![\d,.]){re.escape(f)}(?![\d,.])", body)
                              for f in forms)
                if not present:
                    continue
                reached = any(re.search(pat, body, re.I)
                              for _, pat in phrasings(fact))
                if not reached:
                    gaps.append(f"{doc} states {fact.key} ({fact.value}) but no "
                                f"phrasing reaches it there")
        self.assertEqual(gaps, [],
                         f"{len(gaps)} figures are stated in a document the gate "
                         f"cannot see:\n  " + "\n  ".join(gaps))

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
        """The gate's own assertion must be what fails, not something incidental.

        Standing rule 10 in PITFALLS says a gate is justified by deleting it and
        watching the suite fail. Round 7 did that and the suite stayed green, so this
        checks the inverse directly: a document with a wrong figure produces
        contradictions, and the unmutated documents produce none.
        """
        self.assertEqual(contradictions(self.facts, self.bodies), [])
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
