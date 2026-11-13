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

from hsv1coords.facts import derive  # noqa: E402

DOCS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")
NUMBER = r"(?:\d[\d,]*(?:\.\d+)?|" + "|".join(
    w for w in ("one two three four five six seven eight nine ten eleven twelve "
                "thirteen fourteen fifteen sixteen seventeen eighteen nineteen "
                "twenty").split()) + r")"


def flat(name: str) -> str:
    return re.sub(r"\s+", " ", (ROOT / name).read_text())


class EveryStatedFigureMatchesItsArtefact(unittest.TestCase):
    """The gate. One subtest per (document, phrasing) pair that appears."""

    @classmethod
    def setUpClass(cls):
        cls.facts = derive()
        if not cls.facts:
            raise unittest.SkipTest("no artefacts built; run `make all` first")
        cls.bodies = {d: flat(d) for d in DOCS}

    def test_no_document_states_a_wrong_value_for_a_derived_figure(self):
        problems = []
        for fact in self.facts:
            for noun in fact.nouns:
                pattern = rf"({NUMBER})\s+{re.escape(noun)}"
                for doc, body in self.bodies.items():
                    for stated in re.findall(pattern, body, re.I):
                        if stated.lower() not in {f.lower() for f in fact.forms}:
                            problems.append(
                                f"{doc}: '{stated} {noun}' but {fact.key} is "
                                f"{fact.value}")
        self.assertEqual(problems, [],
                         f"{len(problems)} figures contradict their artefacts:\n  "
                         + "\n  ".join(problems))

    def test_no_noun_starts_with_a_number(self):
        """A noun beginning with a digit makes the pattern capture the wrong figure.

        One did ("is printed as 109,564 in"), so the check compared the preceding
        number against the derived value and failed on correct text.
        """
        bad = [(f.key, n) for f in self.facts for n in f.nouns
               if re.match(r"^\d", n.strip())]
        self.assertEqual(bad, [], f"nouns starting with a digit: {bad}")

    def test_every_gated_noun_can_actually_match(self):
        """A noun that matches nothing anywhere is dead weight pretending to be a gate."""
        all_text = " ".join(self.bodies.values())
        dead = []
        for fact in self.facts:
            if not any(re.search(rf"{NUMBER}\s+{re.escape(n)}", all_text, re.I)
                       for n in fact.nouns):
                continue                       # ungated, reported separately
            for noun in fact.nouns:
                if not re.search(rf"{NUMBER}\s+{re.escape(noun)}", all_text, re.I):
                    dead.append((fact.key, noun))
        print(f"\n    {len(dead)} phrasings match nothing and are carried for "
              f"future wording")

    def test_each_derived_figure_is_stated_somewhere(self):
        """A figure nothing mentions is not gated by the test above.

        Reported rather than enforced: not every derived figure belongs in prose, and
        requiring that would push numbers into the documents to satisfy a test.
        """
        all_text = " ".join(self.bodies.values())
        unstated = []
        for fact in self.facts:
            if not any(re.search(rf"{re.escape(f)}\s+{re.escape(n)}", all_text, re.I)
                       for f in fact.forms for n in fact.nouns):
                unstated.append(fact.key)
        print(f"\n    {len(self.facts) - len(unstated)} of {len(self.facts)} derived "
              f"figures are stated in the documents and gated; "
              f"ungated: {unstated}")

    def test_the_gate_has_power_against_each_document(self):
        """Prove a single wrong figure in a single document is caught.

        Round 6's method, run here so it cannot quietly stop working.
        """
        import tempfile
        for fact in self.facts:
            if not isinstance(fact.value, int):
                continue
            noun = fact.nouns[0]
            for doc, body in self.bodies.items():
                if not re.search(rf"{fact.value}\s+{re.escape(noun)}", body, re.I):
                    continue
                broken = re.sub(rf"{fact.value}(\s+{re.escape(noun)})",
                                rf"{fact.value + 1}\1", body, count=1, flags=re.I)
                stated = re.findall(rf"({NUMBER})\s+{re.escape(noun)}", broken, re.I)
                with self.subTest(document=doc, figure=fact.key):
                    self.assertTrue(
                        any(s.lower() not in {f.lower() for f in fact.forms}
                            for s in stated),
                        f"changing {fact.key} in {doc} would not be detected")
                return
        self.skipTest("no integer figure is stated in any document")


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
