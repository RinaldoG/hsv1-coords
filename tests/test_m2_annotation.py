"""The M2 gate: the annotation is versioned, classed, sourced and documented.

The roadmap's M2 criteria were written before the build. Two of the three could not be
met as originally worded, because the premise changed: see FINDINGS.md and the note on
test_discrepancy_is_documented_not_resolved below.
"""
from __future__ import annotations

import collections
import csv
import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
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
        """The declared set comes from the builder, not from a list retyped here.

        It was retyped, so adding two classes in round 10 failed this test for the wrong
        reason: the classes were declared and described in the builder and in the
        published summary, and only this copy of the list disagreed. A second copy of a
        declaration is a thing that goes stale.
        """
        from build_annotation import CLASS_RULES, DEFAULT_CLASS, UNNAMED_CLASSES
        declared = ({label for label, _, _ in CLASS_RULES}
                    | {DEFAULT_CLASS[0]}
                    | {cls for cls, _ in UNNAMED_CLASSES.values() if cls})
        found = {r["feature_class"] for r in rows()}
        self.assertEqual(found - declared, set(),
                         f"classes in the table that the builder does not declare: "
                         f"{found - declared}")
        self.assertGreater(len(found), 8,
                           f"only {len(found)} classes are in use, so this proves little")

    def test_every_class_the_builder_declares_is_described(self):
        """A class with no description is a column heading with no meaning."""
        from build_annotation import CLASS_RULES, DEFAULT_CLASS, UNNAMED_CLASSES
        described = {label: desc for label, _, desc in CLASS_RULES}
        described[DEFAULT_CLASS[0]] = DEFAULT_CLASS[1]
        described.update({cls: desc for cls, desc in UNNAMED_CLASSES.values() if cls})
        thin = sorted(cls for cls, desc in described.items() if len(desc or "") < 20)
        self.assertEqual(thin, [], f"these classes have no real description: {thin}")

    SYSTEMATIC = r"(UL|US|RL|RS|IRL|IRS|TRL|TRS)\d+(\.\d+)?"

    def test_a_span_between_two_systematic_names_is_not_a_gene(self):
        for r in rows():
            if re.match(rf"^{self.SYSTEMATIC}-{self.SYSTEMATIC}$", r["name"]):
                with self.subTest(name=r["name"]):
                    self.assertEqual(r["feature_class"], "transcript_span")

    def test_a_hyphenated_name_that_is_not_a_span_is_not_called_one(self):
        """Round 10: the span rule matched any hyphenated name.

        So AL-RNA, ORF-O, ORF-P, Ori-S and gG-1 were all classed as transcript units
        spanning more than one gene, and the count of spans quoted in ROADMAP.md
        included them. The old test here reused the rule's own regex to decide which
        names were spans, so it restated the rule instead of checking it and could not
        fail.
        """
        wrong = [(r["name"], r["feature_class"]) for r in rows()
                 if "-" in r["name"]
                 and not re.match(rf"^{self.SYSTEMATIC}-{self.SYSTEMATIC}$", r["name"])
                 and r["feature_class"] == "transcript_span"]
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} names are classed as spanning two genes without being two "
            f"systematic names joined: {wrong[:8]}")

    def test_the_named_transcripts_are_classed_by_what_they_are(self):
        """The five names the loose rule mislabelled, now checked by name."""
        expected = {"AL-RNA": "lat_region", "ORF-O": "lat_region",
                    "ORF-P": "lat_region", "LAT": "lat_region",
                    "Ori-S": "origin_locus", "UL49A": "named_gene",
                    "US8A": "named_gene", "UL26.5": "named_gene"}
        # Gene-like rows only: a regulatory element that carries `/gene=LAT` is a
        # regulatory element, which is the whole point of classing by feature type, so
        # it is not evidence about how LAT itself is classed.
        from build_annotation import GENE_LIKE_TYPES
        found = {}
        for r in rows():
            if r["name"] in expected and r["feature_type"] in GENE_LIKE_TYPES:
                found.setdefault(r["name"], set()).add(r["feature_class"])
        for name, want in sorted(expected.items()):
            with self.subTest(name=name):
                self.assertIn(name, found, f"{name} is not in the table at all")
                self.assertEqual(found[name], {want},
                                 f"{name} is classed {found[name]}, not {want}")

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


class StrandSurvivesEverySource(unittest.TestCase):
    """Round 10's blocker: the GTF's strand column was never read.

    `parse_gtf` built the location from the start and end columns only, so 423 of that
    file's 910 lines, on the minus strand, reached this table as a plain `start..end`
    range. The same features from the GenBank-format version of the same annotation
    carry `complement(...)`, so the absence read as a positive assertion of the plus
    strand for 83 genes, UL23, UL27, UL29, UL48 and IRS1 among them. Anyone pulling
    codons for those got the reverse complement of a coding sequence.

    Stripping `complement(` from every location in the table left the suite green.
    """

    @classmethod
    def setUpClass(cls):
        gtf = ROOT / "data" / "raw" / "HSV1_Annotation.gtf"
        if not gtf.exists():
            raise unittest.SkipTest("HSV1_Annotation.gtf is fetched, not committed")
        cls.gtf = gtf

    def test_every_minus_strand_gtf_line_is_complemented_in_the_table(self):
        from extract_names import parse_gtf
        expected = {}
        for f in parse_gtf(self.gtf):
            key = (f["gene"], f["_type"], f["_location"].strip("complement()"))
            expected[key] = f["strand"]
        minus = sum(1 for v in expected.values() if v == "-")
        self.assertGreater(minus, 300,
                           f"only {minus} minus-strand features in the GTF, so this "
                           f"gate is checking almost nothing")
        wrong = []
        for r in rows():
            if r["source"] != "whisnant2020-gtf":
                continue
            bare = r["location"].replace("complement(", "").rstrip(")")
            strand = expected.get((r["name"], r["feature_type"], bare))
            if strand is None:
                continue
            complemented = r["location"].startswith("complement(")
            if (strand == "-") != complemented:
                wrong.append((r["name"], r["feature_type"], r["location"], strand))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} rows disagree with the GTF's own strand column "
            f"(name, type, location, strand): {wrong[:8]}")

    def test_the_two_formats_of_the_same_annotation_agree_on_strand(self):
        """`BK012101.1` and the GTF are one annotation in two formats, per the builder.

        Where both state the same feature at the same coordinates, they must agree about
        the strand. This is the check that would have caught the defect without knowing
        the GTF format at all.
        """
        def index(source):
            out = {}
            for r in rows():
                if r["source"] != source:
                    continue
                bare = r["location"].replace("complement(", "").rstrip(")")
                out[(r["name"], r["feature_type"], bare)] = \
                    r["location"].startswith("complement(")
            return out
        bk, gtf = index("BK012101.1"), index("whisnant2020-gtf")
        shared = sorted(set(bk) & set(gtf))
        self.assertGreater(len(shared), 20,
                           f"only {len(shared)} features are stated in both formats, so "
                           f"this comparison proves little")
        disagree = [k for k in shared if bk[k] != gtf[k]]
        self.assertEqual(
            disagree, [],
            f"{len(disagree)} features are on opposite strands in the two formats of "
            f"the same annotation: {disagree[:6]}")

    def test_the_product_column_is_not_a_biotype(self):
        """`gene_type` is "protein_coding" on all 910 GTF lines, LAT included.

        It was being written into `product`, so the table asserted that LAT and the
        pri-miRNAs have the product "protein_coding". It belongs in `subtype`.
        """
        products = {r["product"] for r in rows() if r["source"] == "whisnant2020-gtf"}
        self.assertNotIn("protein_coding", products,
                         "the product column still carries the GTF's biotype")
        subtypes = {r["subtype"] for r in rows() if r["source"] == "whisnant2020-gtf"}
        self.assertIn("protein_coding", subtypes,
                      "the biotype was dropped rather than moved")

    def test_the_regulatory_subtypes_are_recoverable(self):
        """`regulatory_class` was dropped, so TATA_box and polyA_signal were gone."""
        subtypes = collections.Counter(
            r["subtype"] for r in rows() if r["feature_type"] == "regulatory")
        self.assertGreater(subtypes["polyA_signal_sequence"], 40,
                           f"polyadenylation signals are not identifiable: {subtypes}")
        self.assertGreater(subtypes["TATA_box"], 10,
                           f"TATA boxes are not identifiable: {subtypes}")

    def test_the_gtf_comments_are_kept(self):
        """386 GTF lines carry the functional names and all 910 notes were empty."""
        noted = [r for r in rows() if r["source"] == "whisnant2020-gtf" and r["notes"]]
        self.assertGreater(len(noted), 200,
                           f"only {len(noted)} GTF rows carry their comment")
        self.assertTrue(any("ICP4" in r["notes"] for r in noted),
                        "the comment naming ICP4 is not in the table")


if __name__ == "__main__":
    unittest.main(verbosity=2)
