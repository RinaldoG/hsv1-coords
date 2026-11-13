"""The M2 gate: the annotation is versioned, classed, sourced and documented.

The roadmap's M2 criteria were written before the build. Two of the three could not be
met as originally worded, because the premise changed: see FINDINGS.md and the note on
test_discrepancy_is_documented_not_resolved below.
"""
from __future__ import annotations

import csv
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
VERSION = "0.1.0"
TABLE = ROOT / "data" / "derived" / f"annotation_v{VERSION}.tsv"
SUMMARY = ROOT / "data" / "derived" / "annotation_summary.md"
CHANGELOG = ROOT / "CHANGELOG.md"


def rows() -> list[dict[str, str]]:
    with TABLE.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Versioning(unittest.TestCase):
    def test_the_release_is_versioned(self):
        self.assertTrue(TABLE.exists(), f"no annotation at {TABLE}")
        self.assertRegex(TABLE.name, r"annotation_v\d+\.\d+\.\d+\.tsv")

    def test_the_changelog_covers_this_version(self):
        self.assertIn(f"## v{VERSION}", CHANGELOG.read_text())

    def test_the_changelog_records_divergences(self):
        text = CHANGELOG.read_text()
        section = text.split(f"## v{VERSION}")[1].split("## v0.0.0")[0]
        self.assertIn("Divergences", section,
                      "a release with no stated divergences is not auditable")


class Provenance(unittest.TestCase):
    def test_every_feature_carries_a_source_and_its_provenance(self):
        missing = [(r["name"], r["source"]) for r in rows()
                   if not r["source"].strip() or not r["provenance"].strip()]
        self.assertEqual(missing, [], f"{len(missing)} features lack provenance")

    def test_every_feature_carries_a_class(self):
        missing = [r["name"] for r in rows() if not r["feature_class"].strip()]
        self.assertEqual(missing, [], f"{len(missing)} features lack a class")

    def test_every_feature_carries_a_location(self):
        missing = [r["name"] for r in rows() if not r["location"].strip()]
        self.assertEqual(missing, [], f"{len(missing)} features lack a location")


class Classification(unittest.TestCase):
    def test_classes_are_from_the_declared_set(self):
        declared = {"sORF", "uORF", "uoORF", "iORF", "dORF", "transcript_span",
                    "novel_orf_2020", "primary_mirna", "named_gene", "other"}
        found = {r["feature_class"] for r in rows()}
        self.assertEqual(found - declared, set(), f"undeclared classes: {found - declared}")

    def test_transcript_spans_are_not_counted_as_genes(self):
        """A span across two genes is not a gene, and the class must say so."""
        for r in rows():
            if re.match(r"^[A-Za-z0-9.]+-[A-Za-z0-9.]+$", r["name"]):
                with self.subTest(name=r["name"]):
                    self.assertNotEqual(r["feature_class"], "named_gene")

    def test_every_class_in_the_summary_appears_in_the_table(self):
        summary = SUMMARY.read_text()
        for cls in sorted({r["feature_class"] for r in rows()}):
            self.assertIn(f"`{cls}`", summary, f"{cls} missing from the summary")


class DocumentedDiscrepancy(unittest.TestCase):
    """The M2 criterion: the 77-against-300 difference is documented, not resolved.

    The criterion as written anticipated arbitrating between two conflicting
    annotations. There is no conflict to arbitrate: one record is the other's successor,
    deposited on the same sequence with a wider scope. So the criterion is met by
    recording what the difference actually is, which is a difference of scope, and by
    making the class column prevent the two counts being compared directly again.
    """

    def test_the_counts_differ_by_class_not_by_contradiction(self):
        by_source_class: dict[tuple[str, str], int] = {}
        for r in rows():
            key = (r["source"], r["feature_class"])
            by_source_class[key] = by_source_class.get(key, 0) + 1
        spans = by_source_class.get(("BK012101.1", "transcript_span"), 0)
        self.assertGreater(spans, 100,
                           "the 2020 record's transcript spans are the bulk of the "
                           "count difference; if this is small, the explanation is wrong")
        self.assertEqual(by_source_class.get(("NC_001806.2", "transcript_span"), 0), 0,
                         "the reference annotation does not describe transcript spans")

    def test_the_discrepancy_is_explained_in_the_changelog(self):
        section = CHANGELOG.read_text().split(f"## v{VERSION}")[1]
        for required in ("77", "300", "transcript"):
            self.assertIn(required, section,
                          f"the changelog does not mention {required!r}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
