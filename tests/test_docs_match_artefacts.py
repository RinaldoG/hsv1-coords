"""Gate the documentation against the artefacts it describes.

Four review rounds found stale figures in the documents, and each round's fix went to
the lines the blocker quoted rather than to the whole file, so the next round found
more. Prose cannot be kept true by attention. These checks fail when a number in a
document stops matching the artefact it claims to describe.

Only load-bearing figures are checked. A test that pins every number would make
ordinary editing painful and would itself go stale.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords import ConventionSet, RepeatMap, VersionCrosswalk  # noqa: E402

DOCS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")


def text(name: str) -> str:
    return (ROOT / name).read_text()


def all_docs() -> str:
    """All four documents, with whitespace flattened.

    Markdown wraps prose, so a figure and its noun are routinely split across a line
    break and an indent. Matching on the raw text made several checks silently
    unsatisfiable: they looked strict and could never find a correct value either.
    """
    import re as _re
    return _re.sub(r"\s+", " ", "\n".join(text(d) for d in DOCS))


def docs_stating(pattern: str) -> list[str]:
    """Which documents mention a figure at all."""
    import re as _re
    return [d for d in DOCS if _re.search(pattern, text(d))]


def assert_no_doc_contradicts(case, label: str, correct: str, wrong: list[str]) -> None:
    """Fail if any single document still carries a superseded figure.

    The first version of this module asserted only that the CORRECT value appeared
    somewhere in the four documents concatenated. Each figure is stated in two to four
    of them, so changing one copy left the check green: it caught rebuild drift, where
    an artefact moves and the prose does not, and missed editing drift, which is the
    failure mode of every review round so far. Checking for the absence of known-wrong
    values is what catches that.
    """
    import re as _re
    for doc in DOCS:
        body = text(doc)
        for bad in wrong:
            # Word-bounded, so "2 diploid locus pairs" does not match inside "32 ...".
            if _re.search(r"(?<![\d.])" + _re.escape(bad), body):
                case.fail(f"{doc} still states {bad} for {label}; it is {correct}")


class FiguresMatch(unittest.TestCase):
    def test_repeat_map_counts(self):
        m = RepeatMap.load(ROOT / "data" / "derived")
        partial = sum(1 for p in m.pairs if not p.exact)
        docs = all_docs()
        self.assertIn(f"{len(m.pairs)} diploid locus pairs", docs,
                      f"the documents do not state the live pair count of {len(m.pairs)}")
        self.assertRegex(docs, rf"{partial} of them\s*marked partial",
                         f"the live partial-pair count is {partial}")
        assert_no_doc_contradicts(
            self, "diploid pairs", f"{len(m.pairs)}",
            [f"{n} diploid locus pairs" for n in range(1, 60) if n != len(m.pairs)])
        assert_no_doc_contradicts(
            self, "partial pairs", f"{partial}",
            [f"{n} of them marked partial" for n in range(1, 20) if n != partial])

    def test_no_document_carries_a_superseded_coverage_figure(self):
        """FINDINGS.md kept 85.8/94.7/99.4 through two rounds that named the line."""
        path = ROOT / "data" / "derived" / "crosswalk_coverage.tsv"
        if not path.exists():
            self.skipTest("crosswalk coverage not built; run `make crosswalk`")
        with path.open() as fh:
            live = {r["accession"]: r["mapped_percent"]
                    for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("orientation") == "reverse_complemented"}
        superseded = ["85.8", "94.7", "99.4"]
        for doc in DOCS:
            body = text(doc)
            for bad in superseded:
                if bad in body and bad not in live.values():
                    self.fail(f"{doc} still carries the superseded figure {bad}; the "
                              f"live values are {sorted(live.values())}")

    def test_the_mutant_count_matches_the_mutants(self):
        """A line rewritten to fix one figure introduced another."""
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "t_m5", ROOT / "tests" / "test_m5_conventions.py")
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except Exception:
            self.skipTest("test module needs fetched inputs to import")
        from hsv1coords import VersionCrosswalk
        v = VersionCrosswalk.load(ROOT / "data" / "derived" / "version_crosswalk.tsv")
        n = len(mod.VersionCrosswalkGate._mutants(v.blocks))
        docs = all_docs().lower()
        words = {10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen"}
        forms = {str(n), words.get(n, str(n))}
        self.assertTrue(
            any(f"{f} deliberately broken crosswalks" in docs for f in forms),
            f"the suite builds {n} mutants and no document says so")
        for other, word in words.items():
            if other == n:
                continue
            self.assertNotIn(f"{word} deliberately broken crosswalks", docs,
                             f"a document states {other} mutants; there are {n}")

    WORDS = {5: "five", 9: "nine", 22: "22"}

    def test_version_crosswalk_counts(self):
        v = VersionCrosswalk.load(ROOT / "data" / "derived" / "version_crosswalk.tsv")
        docs = all_docs().lower()
        for n, label in ((len(v.blocks), "blocks"),
                         (len(v.distinct_offsets), "distinct offsets")):
            forms = {str(n), self.WORDS.get(n, str(n))}
            with self.subTest(figure=label):
                self.assertTrue(any(f"{f} {label}" in docs for f in forms),
                                f"no document states {n} {label}")

    def test_convention_count_and_refusal_state(self):
        cs = ConventionSet.load()
        docs = all_docs().lower()
        n = len(cs.conventions)
        forms = {str(n), self.WORDS.get(n, str(n))}
        self.assertTrue(any(f"{f} handlings" in docs for f in forms),
                        f"no document states {n} handlings")
        refused = cs.report(60000).refused
        if not refused:
            self.assertNotIn("one of them refused", docs,
                             "no convention is refused any more")

    def test_reverse_complemented_percentages(self):
        path = ROOT / "data" / "derived" / "crosswalk_coverage.tsv"
        if not path.exists():
            self.skipTest("crosswalk coverage not built; run `make crosswalk`")
        with path.open() as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("orientation") == "reverse_complemented"]
        docs = all_docs()
        for r in rows:
            with self.subTest(accession=r["accession"]):
                self.assertIn(r["mapped_percent"], docs,
                              f"{r['accession']} maps at {r['mapped_percent']} per cent "
                              f"and no document says so")

    def test_the_alignment_preset_is_not_contradicted(self):
        """The generated report once named a preset the builder had stopped using."""
        from build_crosswalk import PRESET
        report = ROOT / "data" / "derived" / "crosswalk_report.md"
        if not report.exists():
            self.skipTest("crosswalk report not built; run `make crosswalk`")
        named = set(re.findall(r"asm\d+", report.read_text()))
        self.assertEqual(named, {PRESET},
                         f"the report names {named} while the builder runs {PRESET}")


class ClaimsMatchBehaviour(unittest.TestCase):
    def test_quoted_output_matches_live_output(self):
        """Every `NOT usable under morse2017_nr: inside X` in a document must be live."""
        cs = ConventionSet.load()
        _, why = cs.usable(118777, "morse2017_nr")
        interval = re.search(r"inside ([\d,]+-[\d,]+)", why).group(1).rstrip(",")
        for doc in DOCS:
            body = text(doc)
            for raw in re.findall(r"NOT usable under morse2017_nr: inside "
                                  r"([\d,]+-[\d,]+)", body):
                quoted = raw.rstrip(",")
                with self.subTest(document=doc, quoted=quoted):
                    self.assertEqual(quoted, interval,
                                     f"{doc} quotes an interval the code no longer "
                                     f"produces")

    def test_no_document_says_a_built_component_is_unbuilt(self):
        """The README said the mask was not built, and listed it as provided below."""
        readme = text("README.md")
        for phrase in ("mask is not built yet", "is not built yet"):
            self.assertNotIn(phrase, readme,
                             "every component listed in the table below is built")

if __name__ == "__main__":
    unittest.main(verbosity=2)
