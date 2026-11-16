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

    def test_an_unknown_feature_type_is_not_classed_from_its_name(self):
        """Round 11, latent: `classify` fell through for an untyped feature type.

        A record carrying `mat_peptide`, `3'UTR` or `polyA_site` with a `/gene`
        qualifier would have been classed `named_gene`, which is exactly what the
        type-first branch was added to prevent. No source in use carries such a type,
        so nothing reached it; this reaches it.
        """
        from build_annotation import classify
        for feature_type in ("mat_peptide", "3'UTR", "polyA_site", "tRNA"):
            with self.subTest(feature_type=feature_type):
                self.assertEqual(
                    classify("UL23", feature_type), "other",
                    f"a {feature_type} feature named UL23 is being classed from its "
                    f"name")
        self.assertEqual(classify("UL23", "CDS"), "named_gene",
                         "a gene-like type must still be classed by name")

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
        """Every one of them, counted from the record rather than from a floor.

        Round 15: this was `assertGreaterEqual(n, 28)` against a record that carries
        31, a budget of three, and dropping three repeat rows while padding the total
        to keep 1,779 passed (standing rule 47: a count is compared against the thing
        it counts).
        """
        record = ROOT / "data" / "raw" / "ref_NC_001806.2.gb"
        if not record.exists():
            raise unittest.SkipTest("ref_NC_001806.2.gb is fetched, not committed")
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        stated = sum(1 for f in parse_genbank(record)
                     if f["_type"] == "repeat_region")
        carried = sum(1 for r in rows()
                      if r["source"] == "NC_001806.2"
                      and r["feature_type"] == "repeat_region")
        self.assertEqual(
            carried, stated,
            f"the record states {stated} repeat_region features and the table carries "
            f"{carried}")

    def test_lat_is_not_claimed_to_be_wholly_inside_the_long_repeat(self):
        """Guards a factual claim that was wrong in the README until 13 Nov 2026.

        LAT copy 2 ends at 127,151 and IRL ends at 126,373, so it extends 778 bp into
        IRS. Any document in this repository asserting otherwise is wrong.
        """
        # Round 16: `claims` was two fixed phrases, and FINDINGS.md itself writes a
        # third wording, "wholly inside the long inverted repeat", which neither
        # matched, so that wording was outside the gate entirely. The rule is a
        # property of the sentence now: "wholly inside" within a few words of LAT and
        # of a repeat, which is what the false claim says however it is phrased
        # (cause 51).
        claim = re.compile(r"wholly inside[^.]{0,60}\brepeat", re.I)
        # Round 15: the exemption was any line containing the bare phrase "does not",
        # or starting with `|`, so "Because LAT sits wholly inside the inverted
        # repeats, the mask does not need IRS" was exempt, and both lines in the live
        # documents were exempt, which left the gate evaluating nothing (cause 60: an
        # exemption with no bound). The retraction markers are now phrases that can
        # only belong to a retraction, and the two lines that carry the old claim are
        # listed and counted, so a third one fails until somebody looks at it.
        markers = ("was wrong", "corrected", "earlier version", "no longer",
                   "retract", "does not lie", "not wholly")
        # Measured 16 November 2026: one line, the retraction table's own row.
        # Round 16: this was a SET of (document, marker), so a second line in the
        # same document with the same marker collapsed into the existing entry and
        # the equality still held: a new table row asserting the retracted claim as
        # current passed (cause 30, a check that can only see what is added). The
        # lines are counted per document now.
        # Measured 16 November 2026: one line, the retraction table's own row, which
        # carries the word "corrected".
        expected_quotations = {("FINDINGS.md", "corrected"): 1}
        import collections as _collections
        offenders, quotations = [], _collections.Counter()
        for doc in ("README.md", "FINDINGS.md", "ROADMAP.md"):
            for n, line in enumerate((ROOT / doc).read_text().splitlines(), 1):
                low = line.lower()
                if "lat" not in low or not claim.search(line):
                    continue
                if re.search(r"\bIRS\b", line) and "778" in line:
                    # The true sentence: LAT runs 778 bp past IRL into IRS.
                    continue
                marker = next((w for w in markers if w in low), None)
                if marker:
                    quotations[(doc, marker)] += 1
                elif line.lstrip().startswith("|"):
                    # A table row in the retraction table quotes the old claim, which
                    # is how a retraction reads.
                    quotations[(doc, "table row")] += 1
                else:
                    offenders.append(f"{doc}:{n}: {line.strip()[:80]}")
        self.assertEqual(offenders, [],
                         "LAT extends 778 bp past IRL into IRS; these lines assert "
                         f"otherwise: {offenders}")
        self.assertEqual(
            dict(quotations), expected_quotations,
            f"the lines quoting the old LAT claim have changed: {dict(quotations)} "
            f"against the recorded {expected_quotations}. Each one is a place a "
            f"reader meets the retracted claim, so they are counted rather than "
            f"collected into a set.")


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


class TheSummaryIsTheTableCountedUp(unittest.TestCase):
    """Round 15: `annotation_summary.md` was opened once, to check class names.

    None of its roughly 80 counts was compared with the table it summarises, which is
    root cause 62's "a summary is checked against what it summarises". Falsifying a
    load-bearing figure in it, `named_gene` 74/6/77 to 7/6/777 and "53 polyadenylation
    signals and 19 TATA boxes" to 9 and 63, left all 307 tests green.

    Both of its tables are recomputed here from `annotation_v0.1.0.tsv`: the first
    counts distinct names per (class, source), the second rows per (type, source).
    """

    CELL = re.compile(r"^\|\s*`([^`]+)`\s*\|([^|]*)\|([^|]*)\|([^|]*)\|")
    SOURCES = ("NC_001806.2", "BK012101.1", "whisnant2020-gtf")

    @classmethod
    def setUpClass(cls):
        summary = ROOT / "data" / "derived" / "annotation_summary.md"
        if not summary.exists():
            raise unittest.SkipTest("annotation_summary.md is built, not committed")
        cls.text = summary.read_text()
        cls.rows = rows()

    def _section(self, heading: str) -> dict:
        """Every labelled row of a section, and a failure for one it cannot read.

        Round 16: this skipped a row whose three cells were not all `isdigit()`, and
        the tests iterated what it returned, so deleting the ``| `named_gene` | 74 |
        6 | 77 |`` row and writing `| `CDS` | **7** | 300 | 309 |` both passed: the
        published summary then said the reference carries 7 CDS features and nothing
        at all about `named_gene`. That is root cause 63, which round 15 logged one
        commit earlier, in round 15's own new gate.
        """
        body = self.text.split(heading, 1)[1].split("\n## ", 1)[0]
        out, unreadable = {}, []
        for line in body.splitlines():
            found = self.CELL.match(line.strip())
            if not found:
                continue
            counts = []
            for cell in found.groups()[1:4]:
                cell = cell.strip()
                counts.append(int(cell) if cell.isdigit() else None)
            if all(c is not None for c in counts):
                out[found.group(1)] = counts
            else:
                unreadable.append((found.group(1), found.groups()[1:4]))
        self.assertEqual(
            unreadable, [],
            f"{len(unreadable)} rows of '{heading}' have a cell that is not a plain "
            f"count, so nothing can compare them with the table: {unreadable}")
        return out

    def test_the_distinct_name_table_is_the_table_counted_up(self):
        stated = self._section("## Distinct names by class")
        self.assertTrue(stated, "no class rows were read out of the summary")
        # Both directions, round 16: a deleted row left the gate rather than failing
        # it (standing rule 62).
        self.assertEqual(
            sorted(stated), sorted({r["feature_class"] for r in self.rows}),
            "the summary's class table and the annotation table name different "
            "classes")
        wrong = []
        for label, counts in stated.items():
            for source, number in zip(self.SOURCES, counts):
                measured = len({r["name"] for r in self.rows
                                if r["feature_class"] == label
                                and r["source"] == source})
                if measured != number:
                    wrong.append((label, source, number, measured))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} cells of the summary's class table are not what the "
            f"annotation table holds (class, source, stated, counted): {wrong[:6]}")

    def test_the_feature_type_table_is_the_table_counted_up(self):
        stated = self._section("## Rows by feature type")
        self.assertTrue(stated, "no feature-type rows were read out of the summary")
        self.assertEqual(
            sorted(stated), sorted({r["feature_type"] for r in self.rows}),
            "the summary's feature-type table and the annotation table name "
            "different types")
        wrong = []
        for label, counts in stated.items():
            for source, number in zip(self.SOURCES, counts):
                measured = sum(1 for r in self.rows
                               if r["feature_type"] == label
                               and r["source"] == source)
                if measured != number:
                    wrong.append((label, source, number, measured))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} cells of the summary's feature-type table are not what the "
            f"annotation table holds (type, source, stated, counted): {wrong[:6]}")

    def test_the_regulatory_counts_in_the_prose_are_the_record_s(self):
        """"53 polyadenylation signals and 19 TATA boxes" is stated twice, typed both
        times: in `build_annotation.py`'s class description and in this summary.
        """
        counts = {}
        for row in self.rows:
            if row["source"] == "NC_001806.2" and row["feature_type"] == "regulatory":
                counts[row["subtype"]] = counts.get(row["subtype"], 0) + 1
        polya = counts.get("polyA_signal_sequence", 0)
        tata = counts.get("TATA_box", 0)
        self.assertIn(f"{polya} are polyadenylation signals and {tata} are TATA boxes",
                      self.text,
                      f"the summary's prose does not state the {polya} polyA signals "
                      f"and {tata} TATA boxes the table holds")


class TheNameArtefactsAgreeWithEachOtherAndTheirReports(unittest.TestCase):
    """Round 16: three artefacts round 15 unstaled still had no gate.

    Nothing read `names_by_source.tsv`, `name_reconciliation.tsv` or
    `mechanical_report.md` beyond an existence check, so four mutations in one run
    passed with 328 tests green: every `whisnant2020` row given
    `products = protein_coding`, which is the exact cause-66 state round 15 fixed;
    `name_reconciliation.tsv` made to claim `UL23`, the aciclovir target, absent from
    NC_001806.2; the report's distinct-name count 232 -> 199; and
    `strain_inventory.md`'s four headline counts replaced by the figures round 15's
    own commit message names as the damage from a hyphen-truncating `strain_of`.

    Cause 66's lesson has two halves, and round 15 did the rebuild half only. These
    are cross-artefact and summary-versus-artefact checks, which need no parser:
    `make names` byte-identity is what CI adds on top.
    """

    @classmethod
    def setUpClass(cls):
        derived = ROOT / "data" / "derived"
        for name in ("names_by_source.tsv", "name_reconciliation.tsv",
                     "mechanical_report.md", "strain_inventory.tsv",
                     "strain_inventory.md"):
            if not (derived / name).exists():
                raise unittest.SkipTest(f"{name} is built, not committed")
        cls.by_source = list(csv.DictReader(
            (derived / "names_by_source.tsv").open(), delimiter="\t"))
        cls.reconciliation = list(csv.DictReader(
            (derived / "name_reconciliation.tsv").open(), delimiter="\t"))
        cls.report = (derived / "mechanical_report.md").read_text()
        cls.inventory = list(csv.DictReader(
            (derived / "strain_inventory.tsv").open(), delimiter="\t"))
        cls.inventory_md = (derived / "strain_inventory.md").read_text()

    def test_the_reconciliation_counts_the_rows_the_other_table_holds(self):
        per_gene: dict[str, set[str]] = {}
        for row in self.by_source:
            per_gene.setdefault(row["gene"], set()).add(row["source"])
        wrong = []
        for row in self.reconciliation:
            gene = row["gene"]
            sources = per_gene.get(gene, set())
            if int(row["n_sources"]) != len(sources):
                wrong.append((gene, row["n_sources"], sorted(sources)))
            for source in ("NC_001806.2", "BK012101.1", "JN555585.1",
                           "whisnant2020"):
                stated = row[source] == "yes"
                if stated != (source in sources):
                    wrong.append((gene, source, row[source], source in sources))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} reconciliation cells disagree with names_by_source.tsv: "
            f"{wrong[:6]}")
        self.assertEqual(
            sorted(per_gene), sorted(r["gene"] for r in self.reconciliation),
            "the two name tables hold different genes")

    def test_the_mechanical_report_is_those_tables_counted_up(self):
        per_gene: dict[str, set[str]] = {}
        for row in self.by_source:
            per_gene.setdefault(row["gene"], set()).add(row["source"])
        in_all = sum(1 for g, s in per_gene.items() if len(s) == 4)
        # The column states the disagreement in full, prefixed `DIFFERS:`.
        differing = sum(1 for r in self.reconciliation
                        if r["product_agreement"].startswith("DIFFERS"))
        for label, number in (
                ("Distinct gene names across all sources", len(per_gene)),
                ("Present in all 4 sources", in_all),
                ("Genes whose `product` text differs between sources", differing)):
            with self.subTest(figure=label):
                self.assertIn(f"{label}: **{number}**", self.report,
                              f"mechanical_report.md does not state {number} for "
                              f"{label!r}, which is what the tables hold")

    def test_the_strain_report_is_the_inventory_counted_up(self):
        # The report counts every distinct value, the unlabelled placeholder
        # included: `fetch_genomes.py` writes `len(by_label)` over all records.
        labels = {r["strain"] for r in self.inventory}
        groups = {r["strain_group"] for r in self.inventory}
        sizes: dict[str, set[str]] = {}
        lengths: dict[str, set[int]] = {}
        for row in self.inventory:
            if row["strain_group"] == "(UNLABELLED)":
                continue
            sizes.setdefault(row["strain_group"], set()).add(row["accession"])
            lengths.setdefault(row["strain_group"], set()).add(int(row["length_bp"]))
        repeated = sum(1 for g, a in sizes.items() if len(a) > 1)
        disputed = sum(1 for g, seen in lengths.items() if len(seen) > 1)
        for label, number in (
                ("Records", len(self.inventory)),
                ("Distinct strain labels, as deposited", len(labels)),
                ("Distinct laboratory strains, after dropping qualifiers such as "
                 "`isolate` and `clone`", len(groups)),
                ("Laboratory strains deposited more than once", repeated)):
            with self.subTest(figure=label):
                self.assertIn(f"{label}: **{number}**", self.inventory_md,
                              f"strain_inventory.md does not state {number} for "
                              f"{label!r}")
        self.assertIn(f"**disagree on length**: **{disputed}**", self.inventory_md)

    def test_no_product_column_carries_a_biotype(self):
        """Cause 66's state, asserted rather than left to a rebuild.

        `protein_coding` is the GTF's `gene_type`, a biotype on all 910 lines, and it
        spent nine rounds in a `products` column where a reader would take it for a
        product.
        """
        offenders = [(r["gene"], r["source"]) for r in self.by_source
                     if "protein_coding" in r["products"]]
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} rows carry the GTF's biotype in their products "
            f"column: {offenders[:6]}")


class TheGenbankParserHandlesTheFormatSOddities(unittest.TestCase):
    """Round 15: two parser cases nothing exercised, both found by reading.

    A record none of the four committed ones resembles reaches them immediately, and
    the project's rule is that a branch nothing exercises is where the next defect
    lives. These are checked against synthetic records written here, which is the only
    way to reach a case the real inputs do not contain.
    """

    @staticmethod
    def _parse(text: str):
        import tempfile
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        path = pathlib.Path(tempfile.mkdtemp()) / "synthetic.gb"
        path.write_text(text)
        return parse_genbank(path)

    def test_a_qualifier_with_no_value_does_not_join_the_previous_one(self):
        features = self._parse(
            "LOCUS       SYNTHETIC   10 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     gene            1..10\n"
            '                     /gene="UL1"\n'
            "                     /pseudo\n"
            '                     /note="after a value-less qualifier"\n'
            "ORIGIN\n        1 acgtacgtac\n//\n")
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "UL1",
                         "a value-less qualifier was appended to the gene name")
        self.assertEqual(features[0].get("pseudo"), "")
        self.assertEqual(features[0]["note"], "after a value-less qualifier")

    def test_a_repeated_qualifier_keeps_both_values(self):
        """Round 15: the second `/rpt_type` overwrote the first on 26 real features."""
        features = self._parse(
            "LOCUS       SYNTHETIC   10 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     repeat_region   1..10\n"
            '                     /rpt_type=direct\n'
            '                     /rpt_type=terminal\n'
            "ORIGIN\n        1 acgtacgtac\n//\n")
        self.assertEqual(features[0]["rpt_type"], "direct; terminal")


class EveryColumnIsTheSourceSOwnValue(unittest.TestCase):
    """Round 15: four columns were checked by nothing, or by a typed floor.

    Each of these left all 307 tests green, with the whole table rebuilt where the
    mutation was in the builder:

    - every one of the 31 reference repeat rows' `subtype` set to `inverted`, so the
      24 tandem reiteration sets and the 2 terminal repeats claim to be inverted
      repeats;
    - the `subtype` column reversed over the 72 `regulatory` rows, which preserves the
      53 polyA and 19 TATA counts exactly and mislabels 18 features, because the only
      gate was two floors (`polyA > 40`, `TATA > 10`);
    - all 599 non-empty `product` values set to `glycoprotein B`;
    - all 693 `transcript_id` values set to one constant;
    - all 72 reference `regulatory` rows reclassed `named_gene`, which is root cause
      36 in the artefact: 53 polyadenylation signals published under a class whose own
      description is "a classically named gene". Round 10 fixed that at `classify()`
      and nothing held the TABLE to it (cause 32).

    So the class is re-derived for every row by calling the classifier, and the three
    value columns are compared with the qualifier the source states for that feature.
    Measured when written: 1,779 of 1,779 rows agree, with no exemption.
    """

    @classmethod
    def setUpClass(cls):
        for name in ("ref_NC_001806.2.gb", "alt_BK012101.1.gb"):
            if not (ROOT / "data" / "raw" / name).exists():
                raise unittest.SkipTest(f"{name} is fetched, not committed")
        sys.path.insert(0, str(ROOT / "src"))
        from build_annotation import SOURCES, pre_2020_names
        from extract_names import parse_genbank, parse_gtf
        cls.known = pre_2020_names(SOURCES, ROOT / "data" / "raw")
        cls.features: dict[tuple, list[dict]] = {}
        for accession, filename, kind, _ in SOURCES:
            path = ROOT / "data" / "raw" / filename
            if not path.exists():
                raise unittest.SkipTest(f"{filename} is fetched, not committed")
            parsed = parse_genbank(path) if kind == "genbank" else parse_gtf(path)
            for feature in parsed:
                key = (accession, feature["_type"], feature["_location"])
                cls.features.setdefault(key, []).append(feature)
        cls.rows = rows()

    def test_every_row_s_class_is_what_the_classifier_gives(self):
        from build_annotation import UNNAMED_CLASSES, classify
        wrong = []
        for row in self.rows:
            name, kind = row["name"], row["feature_type"]
            if name.startswith(f"{kind}@"):
                # The builder's placeholder for a feature the record leaves unnamed.
                expected = UNNAMED_CLASSES.get(kind, ("", ""))[0]
            else:
                expected = classify(name, kind, self.known)
            if expected != row["feature_class"]:
                wrong.append((name, kind, row["feature_class"], expected))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} of {len(self.rows)} rows carry a class the classifier does "
            f"not give for that name and type (name, type, stated, classifier): "
            f"{wrong[:6]}")

    NON_GENE_TYPES = frozenset({"regulatory", "repeat_region", "stem_loop",
                                "rep_origin", "misc_feature"})
    GENE_LIKE_CLASSES = frozenset({"named_gene", "novel_orf_2020", "transcript_span",
                                   "lat_region", "origin_locus", "sORF", "uORF",
                                   "uoORF", "iORF", "dORF", "primary_mirna"})

    def test_a_feature_s_type_decides_whether_its_class_can_be_gene_like(self):
        """Round 16: the class check called the builder's own expression.

        `test_every_row_s_class_is_what_the_classifier_gives` re-derives through
        `classify`, which is character for character what the builder calls, so a
        one-word change to the table BOTH paths read passed: setting
        `UNNAMED_CLASSES["regulatory"]` to `named_gene` and rebuilding published all
        72 reference regulatory rows, 53 of them polyadenylation signals, under a
        class described as "a classically named gene", with 328 tests green. That is
        root cause 36 in the artefact and cause 53's lesson about a gate sharing the
        builder's rule.

        This is written from the artefact instead: a feature whose TYPE is not
        gene-like cannot carry a gene-like class, whatever any table says, and its
        class is its own type. The two tests are kept together deliberately, one
        agreeing with the classifier and one not depending on it.
        """
        wrong = [(r["name"], r["feature_type"], r["feature_class"]) for r in self.rows
                 if r["feature_type"] in self.NON_GENE_TYPES
                 and r["feature_class"] != r["feature_type"]]
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} rows whose feature type is not gene-like carry a class "
            f"that is not their own type: {wrong[:6]}")
        gene_like = {r["feature_class"] for r in self.rows
                     if r["feature_type"] not in self.NON_GENE_TYPES}
        self.assertEqual(
            sorted(gene_like - self.GENE_LIKE_CLASSES), [],
            f"a gene-like feature carries a class this test does not know: "
            f"{sorted(gene_like - self.GENE_LIKE_CLASSES)}. If it is a new class, "
            f"add it here and say in the commit what kind of thing it is.")
        self.assertEqual(
            sorted(self.NON_GENE_TYPES & gene_like), [],
            "a class name is being used both for a non-gene feature type and for a "
            "gene-like feature, so the column no longer separates them")

    def test_no_repeat_region_is_named_after_a_gene_it_sits_inside(self):
        """M4's rule, applied to M2's rows: round 16 found only M4 had it.

        `RL1 repeat_region 988..1040`, `RL2 repeat_region 2466..2624` and `US1
        repeat_region 132391..132516` were published under the names of the genes
        they sit inside, while `repeat_map.tsv` names the identical three features
        from their own notes. M4 has a test forbidding it (root cause 46) and M2 did
        not share the rule, in the one subject both artefacts describe.
        """
        genes = {r["name"] for r in self.rows
                 if r["feature_type"] not in self.NON_GENE_TYPES}
        offenders = []
        for row in self.rows:
            if row["feature_type"] != "repeat_region":
                continue
            if row["name"] in genes:
                offenders.append((row["name"], row["location"], row["notes"][:40]))
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} repeat regions are named after a gene that is also a "
            f"gene-like row, so a reader cannot tell the repeat from the gene inside "
            f"it: {offenders}")

    def test_each_class_description_matches_the_rows_filed_under_it(self):
        """A description that contradicts its own rows is how the reclass published.

        The summary prints one description per class, and after the mutation the
        `named_gene` row carried "A regulatory element. In this reference record 53
        are polyadenylation signals and 19 are TATA boxes".
        """
        from build_annotation import UNNAMED_CLASSES
        for kind, (label, description) in UNNAMED_CLASSES.items():
            if not label:
                # A type the builder deliberately drops, such as the whole-genome
                # `source` record, carries an empty class and no rows.
                self.assertEqual(
                    [r["name"] for r in self.rows if r["feature_type"] == kind], [],
                    f"{kind} is declared as dropped and has rows in the table")
                continue
            with self.subTest(feature_type=kind):
                self.assertEqual(
                    label, kind,
                    f"the class for feature type {kind!r} is {label!r}; an unnamed "
                    f"feature's class is its own type, so that a description cannot "
                    f"end up over rows of a different kind")
                self.assertTrue(description.strip())

    def test_subtype_product_and_transcript_are_the_source_s_qualifiers(self):
        from build_annotation import squeeze, subtype_of
        wrong, checked = [], 0
        for row in self.rows:
            stated = self.features.get(
                (row["source"], row["feature_type"], row["location"]))
            if not stated:
                wrong.append((row["name"], row["feature_type"], row["location"],
                              "the source states no such feature"))
                continue
            checked += 1
            for column, derive in (
                    ("subtype", subtype_of),
                    ("product", lambda f: squeeze(f.get("product") or "")),
                    ("transcript_id", lambda f: squeeze(f.get("transcript") or ""))):
                if not any(derive(f) == row[column] for f in stated):
                    wrong.append((row["name"], column, row[column],
                                  [derive(f) for f in stated][:2]))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} cells are not what the source states for that feature "
            f"(name, column, stated, source): {wrong[:6]}")
        self.assertEqual(
            checked, len(self.rows),
            f"{len(self.rows) - checked} rows were not compared with a source feature")


class EveryGtfRowIsALineOfTheGtf(unittest.TestCase):
    """Round 15: the 910 GTF-sourced rows, 51 per cent of the table, were ungated.

    The strand gate looked each row up by (name, type, bare location) and skipped any
    row it could not find, so a wrong name or a moved location left the gate rather
    than failing it. The cross-format check sees only the 50 rows whose triple also
    appears in `BK012101.1`, and 431 rows are `exon` or `transcript`, which that
    record never states, so they could never enter it. 848 rows spanning 2,005,286
    bases, 277 of them CDS, rested on an eight-name hand list.

    What that cost, each replayed with the whole suite and each green before this:
    swapping the names `UL23` and `UL27` over 12 rows, which puts gB's name on 1,306
    bp of thymidine kinase; moving four `UL23` CDS rows 5,000 bp; and two builder
    mutations, one shifting every exon and transcript by 137 bp (431 rows) and one
    taking each row's `gene_id` from the previous line (169 rows pointing at the
    neighbouring gene). The builder cases are the serious ones: CI rebuilds the
    artefact, so its rebuild-and-diff step reproduces the corrupted table and passes
    too.

    This reads the GTF's own columns rather than `extract_names.parse_gtf`, because a
    gate written from the parser it checks shares the parser's defects (root cause
    53). Measured when written: 910 lines, 910 rows, 0 mismatches, no exemption.
    """

    GENE_ID = re.compile(r'gene_id "([^"]+)"')
    SPAN = re.compile(r"^(\d+)\.\.(\d+)$")

    @classmethod
    def setUpClass(cls):
        gtf = ROOT / "data" / "raw" / "HSV1_Annotation.gtf"
        if not gtf.exists():
            raise unittest.SkipTest("HSV1_Annotation.gtf is fetched, not committed")
        cls.lines: collections.Counter = collections.Counter()
        cls.seqnames: set[str] = set()
        for line in gtf.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            seqname, _, kind, start, end, _, strand = fields[:7]
            attributes = fields[8] if len(fields) > 8 else ""
            gene = self_gene = cls.GENE_ID.search(attributes)
            cls.seqnames.add(seqname)
            cls.lines[(kind, int(start), int(end), strand,
                       gene.group(1) if gene else None)] += 1
        cls.rows = [r for r in rows() if r["source"] == "whisnant2020-gtf"]

    @classmethod
    def _key(cls, row) -> tuple | None:
        bare = row["location"].replace("complement(", "").rstrip(")")
        span = cls.SPAN.match(bare)
        if not span:
            return None
        strand = "-" if row["location"].startswith("complement(") else "+"
        return (row["feature_type"], int(span.group(1)), int(span.group(2)), strand,
                row["name"])

    def test_the_rows_are_the_file_s_lines_as_a_multiset(self):
        """Round 16: the first version of this unioned the names at each key.

        It keyed on (type, span, strand) and collected every `gene_id` AND
        `transcript_id` at that key into a set, so 21 keys carry more than one line,
        691 of the 910 rows sat at a key accepting more than one name, and 18 rows
        accepted another GENE's name. Relabelling `IRL1` as `IRL2` and `TRL1` as
        `TRL2` over four exons, 4,996 bp inside the long inverted repeats, passed with
        all 328 tests green. Renaming all 230 `exon` rows to their `transcript_id`
        passed this gate too, and the missing-line check compared key SETS, so a line
        dropped at a duplicated key was invisible.

        The rows are the file's lines, counted: one key per line including the
        `gene_id`, compared as a multiset in both directions. Measured when written:
        910 lines, 910 rows, the two multisets equal, no union and no exemption.
        """
        counted: collections.Counter = collections.Counter()
        unparsed = []
        for row in self.rows:
            key = self._key(row)
            if key is None:
                unparsed.append((row["name"], row["location"]))
            else:
                counted[key] += 1
        self.assertEqual(unparsed, [],
                         f"{len(unparsed)} GTF-sourced rows have a location this gate "
                         f"cannot read: {unparsed[:5]}")
        missing = sorted((counted - self.lines).items())
        invented = sorted((self.lines - counted).items())
        self.assertEqual(
            (missing, invented), ([], []),
            f"the table and HSV1_Annotation.gtf are not the same multiset of "
            f"(type, span, strand, gene name): {len(missing)} in the table and not "
            f"the file {missing[:4]}, {len(invented)} in the file and not the table "
            f"{invented[:4]}")

    def test_the_row_count_is_the_file_s_line_count(self):
        """Derived from the file, not a typed floor (standing rule 47)."""
        self.assertEqual(
            len(self.rows), sum(self.lines.values()),
            f"{len(self.rows)} rows carry source `whisnant2020-gtf` and the file has "
            f"{sum(self.lines.values())} feature lines")
        self.assertGreater(sum(self.lines.values()), 0,
                           "the GTF holds no feature lines")

    def test_the_file_states_one_sequence_and_the_table_says_which(self):
        """`parse_gtf` never reads column 1, so the frame is an unasserted default.

        The GTF's coordinates are JN555585's, which is byte-identical to
        NC_001806.2, and that identity is what makes these rows comparable with the
        reference's. A second seqname would silently mix frames (root cause 6b), and
        a parser that does not read the field cannot notice.
        """
        self.assertEqual(sorted(self.seqnames), ["JN555585"],
                         f"the GTF states more than one sequence: {self.seqnames}")
        provenance = {r["provenance"] for r in self.rows}
        self.assertTrue(
            all("JN555585" in p for p in provenance),
            f"a GTF-sourced row does not name the sequence its coordinates are in: "
            f"{sorted(provenance)[:2]}")


class TheGenbankParserHandlesTheFormatSOddities(unittest.TestCase):
    """Round 15: two parser cases nothing exercised, both found by reading.

    A record none of the four committed ones resembles reaches them immediately, and
    the project's rule is that a branch nothing exercises is where the next defect
    lives. These are checked against synthetic records written here, which is the only
    way to reach a case the real inputs do not contain.
    """

    @staticmethod
    def _parse(text: str):
        import tempfile
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        path = pathlib.Path(tempfile.mkdtemp()) / "synthetic.gb"
        path.write_text(text)
        return parse_genbank(path)

    def test_a_qualifier_with_no_value_does_not_join_the_previous_one(self):
        features = self._parse(
            "LOCUS       SYNTHETIC   10 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     gene            1..10\n"
            '                     /gene="UL1"\n'
            "                     /pseudo\n"
            '                     /note="after a value-less qualifier"\n'
            "ORIGIN\n        1 acgtacgtac\n//\n")
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0]["gene"], "UL1",
                         "a value-less qualifier was appended to the gene name")
        self.assertEqual(features[0].get("pseudo"), "")
        self.assertEqual(features[0]["note"], "after a value-less qualifier")

    def test_a_repeated_qualifier_keeps_both_values(self):
        """Round 15: the second `/rpt_type` overwrote the first on 26 real features."""
        features = self._parse(
            "LOCUS       SYNTHETIC   10 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     repeat_region   1..10\n"
            '                     /rpt_type=direct\n'
            '                     /rpt_type=terminal\n'
            "ORIGIN\n        1 acgtacgtac\n//\n")
        self.assertEqual(features[0]["rpt_type"], "direct; terminal")


class EveryColumnIsTheSourceSOwnValue(unittest.TestCase):
    """Round 15: four columns were checked by nothing, or by a typed floor.

    Each of these left all 307 tests green, with the whole table rebuilt where the
    mutation was in the builder:

    - every one of the 31 reference repeat rows' `subtype` set to `inverted`, so the
      24 tandem reiteration sets and the 2 terminal repeats claim to be inverted
      repeats;
    - the `subtype` column reversed over the 72 `regulatory` rows, which preserves the
      53 polyA and 19 TATA counts exactly and mislabels 18 features, because the only
      gate was two floors (`polyA > 40`, `TATA > 10`);
    - all 599 non-empty `product` values set to `glycoprotein B`;
    - all 693 `transcript_id` values set to one constant;
    - all 72 reference `regulatory` rows reclassed `named_gene`, which is root cause
      36 in the artefact: 53 polyadenylation signals published under a class whose own
      description is "a classically named gene". Round 10 fixed that at `classify()`
      and nothing held the TABLE to it (cause 32).

    So the class is re-derived for every row by calling the classifier, and the three
    value columns are compared with the qualifier the source states for that feature.
    Measured when written: 1,779 of 1,779 rows agree, with no exemption.
    """

    @classmethod
    def setUpClass(cls):
        for name in ("ref_NC_001806.2.gb", "alt_BK012101.1.gb"):
            if not (ROOT / "data" / "raw" / name).exists():
                raise unittest.SkipTest(f"{name} is fetched, not committed")
        sys.path.insert(0, str(ROOT / "src"))
        from build_annotation import SOURCES, pre_2020_names
        from extract_names import parse_genbank, parse_gtf
        cls.known = pre_2020_names(SOURCES, ROOT / "data" / "raw")
        cls.features: dict[tuple, list[dict]] = {}
        for accession, filename, kind, _ in SOURCES:
            path = ROOT / "data" / "raw" / filename
            if not path.exists():
                raise unittest.SkipTest(f"{filename} is fetched, not committed")
            parsed = parse_genbank(path) if kind == "genbank" else parse_gtf(path)
            for feature in parsed:
                key = (accession, feature["_type"], feature["_location"])
                cls.features.setdefault(key, []).append(feature)
        cls.rows = rows()

    def test_every_row_s_class_is_what_the_classifier_gives(self):
        from build_annotation import UNNAMED_CLASSES, classify
        wrong = []
        for row in self.rows:
            name, kind = row["name"], row["feature_type"]
            if name.startswith(f"{kind}@"):
                # The builder's placeholder for a feature the record leaves unnamed.
                expected = UNNAMED_CLASSES.get(kind, ("", ""))[0]
            else:
                expected = classify(name, kind, self.known)
            if expected != row["feature_class"]:
                wrong.append((name, kind, row["feature_class"], expected))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} of {len(self.rows)} rows carry a class the classifier does "
            f"not give for that name and type (name, type, stated, classifier): "
            f"{wrong[:6]}")

    NON_GENE_TYPES = frozenset({"regulatory", "repeat_region", "stem_loop",
                                "rep_origin", "misc_feature"})
    GENE_LIKE_CLASSES = frozenset({"named_gene", "novel_orf_2020", "transcript_span",
                                   "lat_region", "origin_locus", "sORF", "uORF",
                                   "uoORF", "iORF", "dORF", "primary_mirna"})

    def test_a_feature_s_type_decides_whether_its_class_can_be_gene_like(self):
        """Round 16: the class check called the builder's own expression.

        `test_every_row_s_class_is_what_the_classifier_gives` re-derives through
        `classify`, which is character for character what the builder calls, so a
        one-word change to the table BOTH paths read passed: setting
        `UNNAMED_CLASSES["regulatory"]` to `named_gene` and rebuilding published all
        72 reference regulatory rows, 53 of them polyadenylation signals, under a
        class described as "a classically named gene", with 328 tests green. That is
        root cause 36 in the artefact and cause 53's lesson about a gate sharing the
        builder's rule.

        This is written from the artefact instead: a feature whose TYPE is not
        gene-like cannot carry a gene-like class, whatever any table says, and its
        class is its own type. The two tests are kept together deliberately, one
        agreeing with the classifier and one not depending on it.
        """
        wrong = [(r["name"], r["feature_type"], r["feature_class"]) for r in self.rows
                 if r["feature_type"] in self.NON_GENE_TYPES
                 and r["feature_class"] != r["feature_type"]]
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} rows whose feature type is not gene-like carry a class "
            f"that is not their own type: {wrong[:6]}")
        gene_like = {r["feature_class"] for r in self.rows
                     if r["feature_type"] not in self.NON_GENE_TYPES}
        self.assertEqual(
            sorted(gene_like - self.GENE_LIKE_CLASSES), [],
            f"a gene-like feature carries a class this test does not know: "
            f"{sorted(gene_like - self.GENE_LIKE_CLASSES)}. If it is a new class, "
            f"add it here and say in the commit what kind of thing it is.")
        self.assertEqual(
            sorted(self.NON_GENE_TYPES & gene_like), [],
            "a class name is being used both for a non-gene feature type and for a "
            "gene-like feature, so the column no longer separates them")

    def test_no_repeat_region_is_named_after_a_gene_it_sits_inside(self):
        """M4's rule, applied to M2's rows: round 16 found only M4 had it.

        `RL1 repeat_region 988..1040`, `RL2 repeat_region 2466..2624` and `US1
        repeat_region 132391..132516` were published under the names of the genes
        they sit inside, while `repeat_map.tsv` names the identical three features
        from their own notes. M4 has a test forbidding it (root cause 46) and M2 did
        not share the rule, in the one subject both artefacts describe.
        """
        genes = {r["name"] for r in self.rows
                 if r["feature_type"] not in self.NON_GENE_TYPES}
        offenders = []
        for row in self.rows:
            if row["feature_type"] != "repeat_region":
                continue
            if row["name"] in genes:
                offenders.append((row["name"], row["location"], row["notes"][:40]))
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} repeat regions are named after a gene that is also a "
            f"gene-like row, so a reader cannot tell the repeat from the gene inside "
            f"it: {offenders}")

    def test_each_class_description_matches_the_rows_filed_under_it(self):
        """A description that contradicts its own rows is how the reclass published.

        The summary prints one description per class, and after the mutation the
        `named_gene` row carried "A regulatory element. In this reference record 53
        are polyadenylation signals and 19 are TATA boxes".
        """
        from build_annotation import UNNAMED_CLASSES
        for kind, (label, description) in UNNAMED_CLASSES.items():
            if not label:
                # A type the builder deliberately drops, such as the whole-genome
                # `source` record, carries an empty class and no rows.
                self.assertEqual(
                    [r["name"] for r in self.rows if r["feature_type"] == kind], [],
                    f"{kind} is declared as dropped and has rows in the table")
                continue
            with self.subTest(feature_type=kind):
                self.assertEqual(
                    label, kind,
                    f"the class for feature type {kind!r} is {label!r}; an unnamed "
                    f"feature's class is its own type, so that a description cannot "
                    f"end up over rows of a different kind")
                self.assertTrue(description.strip())

    def test_subtype_product_and_transcript_are_the_source_s_qualifiers(self):
        from build_annotation import squeeze, subtype_of
        wrong, checked = [], 0
        for row in self.rows:
            stated = self.features.get(
                (row["source"], row["feature_type"], row["location"]))
            if not stated:
                wrong.append((row["name"], row["feature_type"], row["location"],
                              "the source states no such feature"))
                continue
            checked += 1
            for column, derive in (
                    ("subtype", subtype_of),
                    ("product", lambda f: squeeze(f.get("product") or "")),
                    ("transcript_id", lambda f: squeeze(f.get("transcript") or ""))):
                if not any(derive(f) == row[column] for f in stated):
                    wrong.append((row["name"], column, row[column],
                                  [derive(f) for f in stated][:2]))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} cells are not what the source states for that feature "
            f"(name, column, stated, source): {wrong[:6]}")
        self.assertEqual(
            checked, len(self.rows),
            f"{len(self.rows) - checked} rows were not compared with a source feature")


class EveryGtfRowIsALineOfTheGtf(unittest.TestCase):
    """Round 15: the 910 GTF-sourced rows, 51 per cent of the table, were ungated.

    The strand gate looked each row up by (name, type, bare location) and skipped any
    row it could not find, so a wrong name or a moved location left the gate rather
    than failing it. The cross-format check sees only the 50 rows whose triple also
    appears in `BK012101.1`, and 431 rows are `exon` or `transcript`, which that
    record never states, so they could never enter it. 848 rows spanning 2,005,286
    bases, 277 of them CDS, rested on an eight-name hand list.

    What that cost, each replayed with the whole suite and each green before this:
    swapping the names `UL23` and `UL27` over 12 rows, which puts gB's name on 1,306
    bp of thymidine kinase; moving four `UL23` CDS rows 5,000 bp; and two builder
    mutations, one shifting every exon and transcript by 137 bp (431 rows) and one
    taking each row's `gene_id` from the previous line (169 rows pointing at the
    neighbouring gene). The builder cases are the serious ones: CI rebuilds the
    artefact, so its rebuild-and-diff step reproduces the corrupted table and passes
    too.

    This reads the GTF's own columns rather than `extract_names.parse_gtf`, because a
    gate written from the parser it checks shares the parser's defects (root cause
    53). Measured when written: 910 lines, 910 rows, 0 mismatches, no exemption.
    """

    ATTRIBUTE = re.compile(r'(gene_id|transcript_id) "([^"]+)"')
    SPAN = re.compile(r"^(\d+)\.\.(\d+)$")

    @classmethod
    def setUpClass(cls):
        gtf = ROOT / "data" / "raw" / "HSV1_Annotation.gtf"
        if not gtf.exists():
            raise unittest.SkipTest("HSV1_Annotation.gtf is fetched, not committed")
        cls.lines = 0
        cls.names_at: dict[tuple, set[str]] = {}
        cls.seqnames: set[str] = set()
        for line in gtf.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            seqname, _, kind, start, end, _, strand = fields[:7]
            attributes = fields[8] if len(fields) > 8 else ""
            cls.lines += 1
            cls.seqnames.add(seqname)
            key = (kind, int(start), int(end), strand)
            cls.names_at.setdefault(key, set()).update(
                value for _, value in cls.ATTRIBUTE.findall(attributes))
        cls.rows = [r for r in rows() if r["source"] == "whisnant2020-gtf"]

    def _key(self, row) -> tuple | None:
        bare = row["location"].replace("complement(", "").rstrip(")")
        span = self.SPAN.match(bare)
        if not span:
            return None
        strand = "-" if row["location"].startswith("complement(") else "+"
        return (row["feature_type"], int(span.group(1)), int(span.group(2)), strand)

    def test_every_row_s_name_type_span_and_strand_are_one_line_of_the_file(self):
        wrong = []
        for row in self.rows:
            key = self._key(row)
            if key is None:
                wrong.append((row["name"], row["location"], "location unparsed"))
                continue
            stated = self.names_at.get(key)
            if stated is None:
                wrong.append((row["name"], row["feature_type"], row["location"],
                              "the file states no line with that type, span and "
                              "strand"))
            elif row["name"] not in stated:
                wrong.append((row["name"], row["feature_type"], row["location"],
                              f"the file calls it {sorted(stated)}"))
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} of {len(self.rows)} GTF-sourced rows are not a line of "
            f"HSV1_Annotation.gtf: {wrong[:8]}")

    def test_the_row_count_is_the_file_s_line_count(self):
        """Derived from the file, not a typed floor (standing rule 47)."""
        self.assertEqual(
            len(self.rows), self.lines,
            f"{len(self.rows)} rows carry source `whisnant2020-gtf` and the file has "
            f"{self.lines} feature lines, so the table has gained or lost features")
        self.assertGreater(self.lines, 0, "the GTF holds no feature lines")

    def test_no_line_of_the_file_is_missing_from_the_table(self):
        """The other direction: a dropped feature is as wrong as an invented one."""
        present = set()
        for row in self.rows:
            key = self._key(row)
            if key is not None:
                present.add(key)
        missing = sorted(set(self.names_at) - present)
        self.assertEqual(
            missing, [],
            f"{len(missing)} features in HSV1_Annotation.gtf have no row in the "
            f"table: {missing[:8]}")

    def test_the_file_states_one_sequence_and_the_table_says_which(self):
        """`parse_gtf` never reads column 1, so the frame is an unasserted default.

        The GTF's coordinates are JN555585's, which is byte-identical to
        NC_001806.2, and that identity is what makes these rows comparable with the
        reference's. A second seqname would silently mix frames (root cause 6b), and
        a parser that does not read the field cannot notice.
        """
        self.assertEqual(sorted(self.seqnames), ["JN555585"],
                         f"the GTF states more than one sequence: {self.seqnames}")
        provenance = {r["provenance"] for r in self.rows}
        self.assertTrue(
            all("JN555585" in p for p in provenance),
            f"a GTF-sourced row does not name the sequence its coordinates are in: "
            f"{sorted(provenance)[:2]}")


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
                # Round 15: this was `continue`, so a row whose name or location the
                # GTF does not state left the gate instead of failing it, and 848 of
                # the 910 GTF-sourced rows were tied to nothing.
                # `EveryGtfRowIsALineOfTheGtf` below is the check that reads them all;
                # this one keeps its own failure rather than passing them on.
                wrong.append((r["name"], r["feature_type"], r["location"],
                              "the GTF states no such feature under that name"))
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

    def test_every_genbank_location_matches_the_record_it_came_from(self):
        """Round 11: the strand gate covered the GTF only.

        Stripping `complement(` from every `NC_001806.2` row left the suite green, 112
        rows changed, including `RL2 gene complement(120675..124287)`. `BK012101.1` was
        protected only by the relative check between the two formats of the 2020
        annotation, which a shared strip passes, and the reference record had no
        counterpart at all: 116 of its 284 features carry `complement`. So the exact
        defect round 10 fixed was still open on the primary annotation.

        This checks each GenBank-sourced location against the record's own feature
        table, which is the only independent source for it.

        Round 12: it checked the location and not the name. A row whose `(type,
        location)` pair existed anywhere in the record passed through a blanket
        fallback, so the name beside the location was never confirmed: rotating the
        `name` column of all 77 reference CDS rows by one position left the suite green
        with 74 of 77 gene names pointing at another gene's start, a median displacement
        of 1,742 bp and 128,676 bp of CDS span relabelled. A two-row swap passed too,
        and the companion count test cannot see a permutation.

        So the name is checked as well, against the only places a name can honestly
        come from: the feature's own `/gene` qualifier, its `/note`, or the first
        `;`-delimited clause of that note (`TRL; inverted repeat flanking UL` is written
        `TRL`). That holds for 868 of the 869 GenBank rows. The one row it does not hold
        for is a TATA box the record names at all, which the builder calls
        `regulatory@complement(7592..7597)`, and that placeholder is allowed only in the
        exact form `{type}@{location}` and counted, rather than waved through by a rule
        that also admits a rotation.
        """
        import collections
        from extract_names import parse_genbank
        records = {"NC_001806.2": "ref_NC_001806.2.gb",
                   "BK012101.1": "alt_BK012101.1.gb"}
        checked, wrong, placeholders = 0, [], []
        for source, filename in sorted(records.items()):
            path = ROOT / "data" / "raw" / filename
            if not path.exists():
                continue
            names: dict[tuple[str, str], set[str]] = collections.defaultdict(set)
            for f in parse_genbank(path):
                location = (f.get("_location") or "").strip()
                if not location:
                    continue
                key = (f["_type"], location)
                names[key]  # a feature the record names at all is still a feature
                note = (f.get("note") or "").strip()
                for value in (f.get("gene"), note, note.split(";")[0]):
                    if value and value.strip():
                        names[key].add(value.strip())
            for r in rows():
                if r["source"] != source:
                    continue
                checked += 1
                key = (r["feature_type"], r["location"])
                if key not in names:
                    wrong.append((source, r["name"], r["feature_type"], r["location"],
                                  "the record states no such feature there"))
                elif r["name"] in names[key]:
                    continue
                elif r["name"] == f"{r['feature_type']}@{r['location']}":
                    placeholders.append((source, r["name"]))
                else:
                    wrong.append((source, r["name"], r["feature_type"], r["location"],
                                  f"the record calls it {sorted(names[key])}"))
        self.assertGreater(checked, 400,
                           f"only {checked} GenBank rows were checked against their "
                           f"records")
        self.assertEqual(
            wrong, [],
            f"{len(wrong)} rows state a name or a location the record does not "
            f"(source, name, type, location, what the record says): {wrong[:6]}")
        self.assertEqual(
            len(placeholders), 1,
            f"{len(placeholders)} rows carry a synthesised placeholder name rather than "
            f"a name the record states: {placeholders}. One is expected (the TATA box "
            f"at complement(7592..7597), which the record leaves unnamed); more means "
            f"names are being generated where the record has them.")

    def test_the_complement_count_matches_the_records(self):
        """A count, so a wholesale strip fails even if a row still matches by type."""
        from extract_names import parse_genbank
        for source, filename in (("NC_001806.2", "ref_NC_001806.2.gb"),
                                 ("BK012101.1", "alt_BK012101.1.gb")):
            path = ROOT / "data" / "raw" / filename
            if not path.exists():
                continue
            expected = sum(1 for f in parse_genbank(path)
                           if "complement" in (f.get("_location") or "")
                           and f["_type"] != "source")
            found = sum(1 for r in rows()
                        if r["source"] == source and "complement" in r["location"])
            with self.subTest(source=source):
                self.assertEqual(
                    found, expected,
                    f"{source}: the table has {found} complemented locations and the "
                    f"record has {expected}")
                self.assertGreater(expected, 50,
                                   f"{source} has only {expected} complemented "
                                   f"features, so this gate checks little")

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
        # Round 15: these were floors of 40 and 10 against 53 and 19, so reversing
        # the subtype column over the 72 regulatory rows preserved both counts exactly
        # and mislabelled 18 features. The counts come from the records now, and
        # `EveryColumnIsTheSourceSOwnValue` holds each row's subtype to its own
        # feature's qualifier, which is what catches a permutation.
        sys.path.insert(0, str(ROOT / "src"))
        from extract_names import parse_genbank
        from build_annotation import subtype_of
        stated = collections.Counter()
        for name in ("ref_NC_001806.2.gb", "alt_BK012101.1.gb"):
            path = ROOT / "data" / "raw" / name
            if not path.exists():
                raise unittest.SkipTest(f"{name} is fetched, not committed")
            for f in parse_genbank(path):
                if f["_type"] == "regulatory":
                    stated[subtype_of(f)] += 1
        for subtype in ("polyA_signal_sequence", "TATA_box"):
            with self.subTest(subtype=subtype):
                self.assertEqual(
                    subtypes[subtype], stated[subtype],
                    f"the records state {stated[subtype]} {subtype} features and the "
                    f"table carries {subtypes[subtype]}")

    def test_the_gtf_comments_are_kept(self):
        """386 GTF lines carry the functional names and all 910 notes were empty."""
        noted = [r for r in rows() if r["source"] == "whisnant2020-gtf" and r["notes"]]
        self.assertGreater(len(noted), 200,
                           f"only {len(noted)} GTF rows carry their comment")
        self.assertTrue(any("ICP4" in r["notes"] for r in noted),
                        "the comment naming ICP4 is not in the table")


if __name__ == "__main__":
    unittest.main(verbosity=2)
