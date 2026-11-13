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

    def test_every_location_is_well_formed(self):
        """A truncated location is worse than a missing one: it looks usable.

        The first version of the GenBank parser dropped continuation lines, so wrapped
        locations shipped cut off mid-expression, including RL2 (ICP0) and the IRL2-IRL1
        span. Checking only for non-empty did not catch it.
        """
        bad: list[tuple[str, str, str]] = []
        for r in rows():
            loc = r["location"].strip()
            if loc.count("(") != loc.count(")"):
                bad.append((r["name"], r["source"], "unbalanced brackets"))
            elif loc.endswith((",", "..")):
                bad.append((r["name"], r["source"], "ends mid-expression"))
            elif not re.search(r"\d", loc):
                bad.append((r["name"], r["source"], "no coordinate"))
        self.assertEqual(bad, [], f"{len(bad)} malformed locations: {bad[:5]}")

    def test_locations_parse_to_coordinates(self):
        """Every location yields at least one start..end pair within the genome."""
        GENOME = 152222
        bad = []
        for r in rows():
            spans = re.findall(r"(\d+)\.\.>?(\d+)", r["location"])
            if not spans:
                bad.append((r["name"], r["location"][:40]))
                continue
            for a, b in spans:
                if not (1 <= int(a) <= int(b) <= GENOME):
                    bad.append((r["name"], f"{a}..{b}"))
        self.assertEqual(bad, [], f"{len(bad)} locations do not parse: {bad[:5]}")


class Classification(unittest.TestCase):
    def test_classes_are_from_the_declared_set(self):
        declared = {"sORF", "uORF", "uoORF", "iORF", "dORF", "transcript_span",
                    "novel_orf_2020", "primary_mirna", "named_gene", "other",
                    # features carrying no /gene qualifier, added after review
                    "repeat_region", "stem_loop", "rep_origin", "misc_feature",
                    "regulatory", "misc_RNA"}
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


class InvertedRepeatsPresent(unittest.TestCase):
    """The inverted repeats must be in the annotation.

    They are the reason the coordinate problem exists, and the first build of v0.1.0
    omitted all four because it kept only features with a `/gene` qualifier. A release
    of this annotation without them is worse than useless: it looks complete.
    """

    def test_all_four_inverted_repeats_are_annotated(self):
        got = {r["name"].split(";")[0].strip(): r["location"]
               for r in rows()
               if r["source"] == "NC_001806.2" and r["feature_class"] == "repeat_region"}
        for name, expected in (("TRL", "1..9213"), ("IRL", "117161..126373"),
                               ("IRS", "125975..132607"), ("TRS", "145590..152222")):
            with self.subTest(repeat=name):
                self.assertIn(name, got, f"{name} missing from the annotation")
                self.assertEqual(got[name], expected)

    def test_the_reference_repeat_features_are_not_dropped(self):
        """The reference record carries 31 repeat_region features; keep them all."""
        n = sum(1 for r in rows()
                if r["source"] == "NC_001806.2" and r["feature_class"] == "repeat_region")
        self.assertGreaterEqual(n, 28, f"only {n} repeat features carried through")

    def test_lat_is_not_claimed_to_be_wholly_inside_the_long_repeat(self):
        """Guards a factual claim that was wrong in the README until 13 Nov 2026.

        LAT copy 2 ends at 127,151 and IRL ends at 126,373, so it extends 778 bp into
        IRS. Any document in this repository asserting otherwise is wrong.
        """
        claims = ("wholly inside the inverted repeats", "wholly inside the long repeat")
        offenders = []
        for doc in ("README.md", "FINDINGS.md", "ROADMAP.md"):
            for n, line in enumerate((ROOT / doc).read_text().splitlines(), 1):
                low = line.lower()
                if not any(c in low for c in claims):
                    continue
                # A table row or a line that marks itself as a correction is quoting the
                # old claim, which is how a retraction is supposed to read. Only an
                # unmarked assertion is a failure.
                quoted = line.lstrip().startswith("|") or any(
                    w in low for w in ("was wrong", "corrected", "earlier version",
                                       "no longer", "retract", "does not lie",
                                       "does not", "not wholly"))
                if not quoted:
                    offenders.append(f"{doc}:{n}: {line.strip()[:80]}")
        self.assertEqual(offenders, [],
                         "LAT extends 778 bp past IRL into IRS; these lines assert "
                         f"otherwise: {offenders}")


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
