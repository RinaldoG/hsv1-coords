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
    return "\n".join(text(d) for d in DOCS)


class FiguresMatch(unittest.TestCase):
    def test_repeat_map_counts(self):
        m = RepeatMap.load(ROOT / "data" / "derived")
        partial = sum(1 for p in m.pairs if not p.exact)
        docs = all_docs()
        self.assertIn(f"{len(m.pairs)} diploid locus pairs", docs,
                      f"the documents do not state the live pair count of {len(m.pairs)}")
        self.assertRegex(docs, rf"{partial} of them\s*\n?\s*marked partial",
                         f"the live partial-pair count is {partial}")

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

    def test_the_retraction_log_is_referenced_and_present(self):
        parent = ROOT.parent / "PITFALLS.md"
        if not parent.exists():
            self.skipTest("PITFALLS.md lives in the private parent folder")
        body = parent.read_text()
        roots = len(re.findall(r"^### ", body, re.M))
        rows = [l for l in body.splitlines()
                if l.startswith("| ") and not l.startswith("|---")
                and not re.match(r"\|\s*(Claim|Gate|Instance|What was wrong|"
                                 r"Earlier claim)\s*\|", l)]
        self.assertIn(f"**{roots} root causes**", body,
                      f"the stated root-cause count does not match the {roots} "
                      f"headings in the file")
        self.assertIn(f"**{len(rows)} retracted claims", body,
                      f"the stated tabulated count does not match the {len(rows)} "
                      f"rows in the file")


if __name__ == "__main__":
    unittest.main(verbosity=2)
