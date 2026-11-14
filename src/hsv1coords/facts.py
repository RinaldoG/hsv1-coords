"""Every load-bearing figure, derived from the artefacts rather than typed.

Six review rounds found stale numbers in the documents. Each round's fix corrected the
figures that round had quoted, and the next round found more, including in lines the
previous round had already named. Round 6 changed a single figure in a single document
ninety times and only five turned the suite red.

The failure is structural, not inattention: a number written into prose has no link to
the thing it describes. So every figure that matters is derived here, each with the
phrasings the documents use for it, and `tests/test_docs_match_artefacts.py` fails when
a document states a different value for any of them.

Adding a figure to a document without adding it here leaves it ungated, which is the
remaining hole and is stated as such in the roadmap rather than papered over.
"""
from __future__ import annotations

import collections
import csv
import pathlib
import statistics
from dataclasses import dataclass, field

ROOT = pathlib.Path(__file__).resolve().parents[2]
DERIVED = ROOT / "data" / "derived"
RAW = ROOT / "data" / "raw"


@dataclass(frozen=True)
class Fact:
    """One derived figure, and how the documents refer to it."""

    key: str
    value: int | float | str
    nouns: tuple[str, ...]
    """Phrases that FOLLOW the figure in prose, e.g. "diploid locus pairs".

    A noun must not contain a number. One that did ("is printed as 109,564 in") made
    the pattern capture the preceding figure instead, so the check compared the wrong
    number and failed on correct text.
    """
    note: str = ""
    prefixes: tuple[str, ...] = ()
    """Phrases that PRECEDE the figure, e.g. "Median mappable fraction".

    Round 7: five phrasings matched nothing because the documents write the words
    before the number, and the matcher only looked after it. They were printed as
    spares carried for future wording while the figure they named was ungated in that
    document. A prefix must not contain a number, for the same reason a noun must not.

    It sits after `note` deliberately: putting it before shifted every positional
    `note` argument into this field, and a string iterates character by character, so
    each fact silently gained a handful of one-character phrasings.
    """
    words: dict[int, str] = field(default_factory=dict)

    @property
    def forms(self) -> set[str]:
        out = {str(self.value)}
        if isinstance(self.value, int):
            out.add(f"{self.value:,}")
            if self.value in SPELLED:
                out.add(SPELLED[self.value])
        return out


SPELLED = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
           8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve",
           13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen",
           17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}


def _rows(path: pathlib.Path) -> list[dict[str, str]]:
    """Rows of a TSV, skipping any `#` preamble.

    `base_agreement.tsv` carries five comment lines before its header, and without this
    the first comment became the header: every row then parsed into one column, the
    figures came out empty, and `derive` raised on an empty sequence instead of
    reporting a missing artefact. A silent version of that would have published three
    facts with no values.
    """
    with path.open() as fh:
        lines = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def missing_inputs() -> list[str]:
    """Which artefacts are absent, so a smaller gate cannot look like a bigger one.

    Round 7: a fresh clone printed "17 of 18 derived figures gated" against "19 of 26"
    with the inputs present, because `derive` drops every fact whose artefact is
    missing. The ratio improved while the gate shrank, with no notice. The test prints
    this list so the number is read in the right context.
    """
    needed = {
        "data/raw/MANIFEST.tsv": RAW / "MANIFEST.tsv",
        "data/raw/ref_NC_001806.2.gb": RAW / "ref_NC_001806.2.gb",
        "data/raw/ref_NC_001806.1.gb": RAW / "ref_NC_001806.1.gb",
        "data/derived/annotation_v0.1.0.tsv": DERIVED / "annotation_v0.1.0.tsv",
        "data/derived/aliases_literature.tsv": DERIVED / "aliases_literature.tsv",
        "data/derived/repeat_map.tsv": DERIVED / "repeat_map.tsv",
        "data/derived/diploid_pairs.tsv": DERIVED / "diploid_pairs.tsv",
        "data/derived/version_crosswalk.tsv": DERIVED / "version_crosswalk.tsv",
        "data/derived/crosswalk_coverage.tsv": DERIVED / "crosswalk_coverage.tsv",
        "data/derived/base_agreement.tsv": DERIVED / "base_agreement.tsv",
    }
    return sorted(name for name, path in needed.items() if not path.exists())


def derive() -> list[Fact]:
    """Every gated figure. Anything unavailable is omitted rather than guessed."""
    facts: list[Fact] = []

    manifest = RAW / "MANIFEST.tsv"
    if manifest.exists():
        facts.append(Fact("inputs", len(_rows(manifest)),
                          ("inputs fetched and checksummed",),
                          "rows in data/raw/MANIFEST.tsv"))

    ann = DERIVED / "annotation_v0.1.0.tsv"
    if ann.exists():
        rows = _rows(ann)
        facts.append(Fact("annotation_rows", len(rows),
                          ("feature rows", "rows across three sources",
                           "rows, up from")))
        # Distinct names, not rows. Round 10: this counted rows, so the figure in
        # ROADMAP.md described how many gene, CDS and mRNA features the spans are
        # stated across rather than how many spans there are, and the rule that
        # classed them matched any hyphenated name, so AL-RNA, ORF-O, ORF-P and Ori-S
        # were counted among them.
        span_rows = [r for r in rows if r["source"] == "BK012101.1"
                     and r["feature_class"] == "transcript_span"]
        facts.append(Fact("transcript_spans", len({r["name"] for r in span_rows}),
                          ("named transcript units spanning",)))
        facts.append(Fact("transcript_span_rows", len(span_rows),
                          ("rows, counting each gene",)))
        nested = {"sORF", "uORF", "uoORF", "iORF", "dORF"}
        facts.append(Fact(
            "nested_orfs",
            len({r["name"] for r in rows if r["source"] == "BK012101.1"
                 and r["feature_class"] in nested}),
            ("nested reading frames",)))

    aliases = DERIVED / "aliases_literature.tsv"
    if aliases.exists():
        rows = _rows(aliases)
        facts.append(Fact("alias_rows", len(rows), ("cited alias rows",)))
        facts.append(Fact("alias_genes",
                          len({r["canonical_systematic_name"] for r in rows}),
                          ("genes in the", "genes in the dictionary")))
        # Identifiers, not citation strings: there are 104 distinct citation strings
        # and 103 distinct identifiers, and the documents mean the latter. Gating this
        # is what surfaced the difference.
        facts.append(Fact("alias_citations",
                          len({r["citation_identifier"] for r in rows
                               if r.get("citation_identifier")}),
                          ("distinct citation identifiers",)))

    repeat_map = DERIVED / "repeat_map.tsv"
    pairs_file = DERIVED / "diploid_pairs.tsv"
    if repeat_map.exists() and pairs_file.exists():
        regions = _rows(repeat_map)
        pairs = _rows(pairs_file)
        partial = sum(1 for p in pairs
                      if str(p.get("offset_delta", "")).startswith("partial"))
        facts.append(Fact("regions", len(regions), ("regions covering",)))
        # Round 7's second sweep caught these two as misses that are ours and
        # derivable, as against the dates, page numbers and DOIs that are not. The
        # union is computed here rather than read from the report, so the figure and
        # the prose cannot drift apart.
        covered: set[int] = set()
        for r in regions:
            covered.update(range(int(r["start"]), int(r["end"]) + 1))
        facts.append(Fact("repeat_covered_bp", len(covered),
                          ("bp of repeated or ambiguous sequence",)))
        facts.append(Fact("repeat_covered_percent",
                          f"{100 * len(covered) / 152222:.1f}",
                          ("per cent of the genome",)))
        facts.append(Fact("diploid_pairs", len(pairs), ("diploid locus pairs",)))
        facts.append(Fact("partial_pairs", partial,
                          ("of them marked partial",)))

    version = DERIVED / "version_crosswalk.tsv"
    if version.exists():
        blocks = _rows(version)
        offsets = {int(b["target_start"]) - int(b["query_start"]) for b in blocks}
        facts.append(Fact("version_blocks", len(blocks),
                          ("aligned blocks", "blocks, from")))
        facts.append(Fact("version_offsets", len(offsets), ("distinct offsets",)))

    # The base-agreement figures. Round 9: README.md stated "97 per cent or better over
    # 717,688 scored positions" as the evidence that the lift-over defect was fixed, and
    # both numbers were the pre-fix measurement: the widened gate scores 1,676,839
    # positions and four genomes fall below 97. Neither figure was derivable from any
    # artefact, so the gate could not see them, while README.md claimed every figure the
    # repository derives is gated. They are derivable now because the measurement is
    # pinned in a committed file.
    agreement = DERIVED / "base_agreement.tsv"
    if agreement.exists():
        rows = [r for r in _rows(agreement) if r.get("scored")]
        scored = sum(int(r["scored"]) for r in rows)
        matching = sum(int(r["matching"]) for r in rows)
        worst = min(100 * int(r["matching"]) / int(r["scored"]) for r in rows)
        facts.append(Fact("scored_positions", scored,
                          ("scored positions",),
                          note="positions compared against the reference"))
        facts.append(Fact("worst_agreement", f"{worst:.1f}",
                          ("per cent or better",),
                          note="the weakest genome's base agreement",
                          prefixes=("agrees with the reference at",)))
        facts.append(Fact("pooled_agreement", f"{100 * matching / scored:.2f}",
                          ("per cent across the corpus",),
                          note="pooled base agreement",
                          prefixes=("pooled agreement of",)))

    # The feature counts the documents compare, and the LAT overhang. Round 9 showed
    # all three could be changed with the suite green: they are the arithmetic the
    # "77 against 300 CDS" section rests on, and 778 bp is the measurement that makes
    # the LAT pair partial rather than exact.
    annotation = DERIVED / "annotation_v0.1.0.tsv"
    if annotation.exists():
        rows = _rows(annotation)
        counts = collections.Counter(
            (r["source"], r["feature_type"]) for r in rows)
        # Phrasings have to tell these four apart: "300 CDS features and 83 gene
        # names" and "77 CDS features and 79 gene features" sit in the same paragraph,
        # so a noun of "CDS features" alone matched the wrong number and reported a
        # contradiction in correct text. A pattern that matches both figures is worse
        # than no pattern: it fails on the truth.
        facts.append(Fact("ref_cds", counts[("NC_001806.2", "CDS")],
                          (),
                          note="CDS features in the reference record",
                          prefixes=("`NC_001806.2` carries",
                                    "against the reference's")))
        facts.append(Fact("ref_genes", counts[("NC_001806.2", "gene")],
                          ("gene features",),
                          note="gene features in the reference record"))
        facts.append(Fact("alt_cds", counts[("BK012101.1", "CDS")],
                          (),
                          note="CDS features in the re-annotation",
                          prefixes=("`BK012101.1` carries", "bp of sequence with")))
        facts.append(Fact("alt_genes", counts[("BK012101.1", "gene")],
                          ("gene names",),
                          note="gene names in the re-annotation"))

    pairs_file = DERIVED / "diploid_pairs.tsv"
    if pairs_file.exists():
        overhangs = []
        for r in _rows(pairs_file):
            if r["locus"] != "LAT" or r["feature_type"] != "ncRNA":
                continue
            len_1 = int(r["copy_1_end"]) - int(r["copy_1_start"])
            len_2 = int(r["copy_2_end"]) - int(r["copy_2_start"])
            overhangs.append(abs(len_2 - len_1))
        if overhangs:
            facts.append(Fact("lat_overhang_bp", max(overhangs),
                              ("bp past `IRL`", "bp into `IRS`",
                               "bp beyond the end of"),
                              note="how far the LAT copy runs past IRL",
                              prefixes=("it extends", "which runs")))

    coverage = DERIVED / "crosswalk_coverage.tsv"
    if coverage.exists():
        rows = _rows(coverage)
        facts.append(Fact("genomes", len(rows),
                          ("complete genomes", "genomes aligned to the reference"),
                          prefixes=("**Built.**",)))
        facts.append(Fact(
            "median_mappable",
            f"{statistics.median(float(r['mapped_percent']) for r in rows):.2f}",
            ("per cent, after review",),
            # The documents write these words BEFORE the number, which is why the two
            # phrasings that named them matched nothing and were printed as spares.
            note="median mappable fraction",
            prefixes=("Median mappable fraction",
                      "median mappable fraction to")))
        rc = [r for r in rows if r.get("orientation") == "reverse_complemented"]
        facts.append(Fact("reverse_complemented", len(rc),
                          ("genomes were previously unmappable",)))
        for r in rc:
            # A distinctive phrasing per accession, so a percentage cannot match any
            # other percentage in the documents. Round 7: the three were stated as a
            # bare list ("79.40, 79.67 and 69.87 per cent"), where no phrasing can pin
            # an individual figure, so the documents now name each accession.
            facts.append(Fact(f"revcomp_{r['accession']}", r["mapped_percent"],
                              (f"per cent for {r['accession']}",),
                              note=f"{r['accession']} mappable fraction",
                              prefixes=(f"`{r['accession']}` now maps at",)))

    # Structural constants of the genome. These are the most load-bearing numbers in
    # the documents and were previously ungated entirely: a mutation sweep changed one
    # figure in one document twenty times and caught none of them, because the gate
    # only knew about counts of derived rows.
    ref = RAW / "ref_NC_001806.2.gb"
    if ref.exists():
        import re as _re
        body = ref.read_text()
        seq = "".join(_re.findall(r"[acgtnACGTN]",
                                  body.split("ORIGIN")[1].split("//")[0]
                                  .replace("\n", "")))
        # One phrasing per place a document states it. Round 7: this figure is the most
        # load-bearing in the repository, was stated five times, and was reachable at
        # exactly one of them, because the coverage check pooled the four documents.
        facts.append(Fact("genome_bp", len(seq),
                          ("bp, and", "bp of sequence",
                           "bp. The offset between them"),
                          prefixes=("A sweep of all", "runs to",
                                    "`NC_001806.2` is")))
        repeats = {}
        for loc, quals in _re.findall(
                r"^     repeat_region\s+(\S+)\n((?:\s+/.*\n)+)", body, _re.M):
            note = _re.search(r'/note="([^"]+)"', quals)
            m = _re.match(r"(\d+)\.\.(\d+)$", loc)
            if note and m and "inverted repeat flanking" in note.group(1):
                repeats[note.group(1).split(";")[0]] = (int(m.group(1)),
                                                        int(m.group(2)))
        # The four repeat lengths are two distinct values, each already derived below as
        # `long_repeat_bp` and `short_repeat_bp`. Round 7: they were carried here as
        # four more facts with an empty phrasing tuple behind `if False`, which made
        # them permanently ungateable and inflated the coverage denominator by four.
        if "TRL" in repeats and "IRS" in repeats:
            facts.append(Fact("long_repeat_bp",
                              repeats["TRL"][1] - repeats["TRL"][0] + 1,
                              ("bases, as are",)))
            facts.append(Fact("short_repeat_bp",
                              repeats["IRS"][1] - repeats["IRS"][0] + 1,
                              ("bp each. They are not identical",)))

    v1 = RAW / "ref_NC_001806.1.gb"
    if v1.exists():
        import re as _re2
        b1 = v1.read_text()
        s1 = "".join(_re2.findall(r"[acgtnACGTN]",
                                  b1.split("ORIGIN")[1].split("//")[0]
                                  .replace("\n", "")))
        # The backticks are part of the text once it is flattened, so a phrasing that
        # strips them matches nothing. Round 7: all three phrasings here did, which left
        # the one figure that distinguishes the two reference versions, the whole
        # subject of root cause 6b, ungated at all four places it is stated.
        facts.append(Fact("genome_v1_bp", len(s1),
                          ("bp and `NC_001806.2` is", "bp; `NC_001806.2` is",
                           "bp, which is `.1`"),
                          prefixes=("`NC_001806.1` is",)))

    try:
        from .conventions import ConventionSet
        cs = ConventionSet.load()
        # The worked example's own numbers, derived from behaviour rather than typed.
        lat = cs["szpara2014_trimmed"].to_convention(118777)
        if lat is not None:
            facts.append(Fact("lat_in_szpara", lat,
                              ("in a Szpara-trimmed analysis",
                               "is computed here from the rule")))
        facts.append(Fact("conventions", len(cs.conventions),
                          ("handlings, all", "handlings, each", "handlings,")))
    except Exception:                                   # artefacts not built
        pass

    return facts


# `superseded_values` was removed in round 8. It generated plausible wrong values for
# an absence check, it had zero callers from the day it was written, and it was the
# helper standing rule 14 named ("a documentation gate asserts the absence of
# superseded values, per document"). Round 7's own root cause 18 lists "a helper was
# added with zero callers, which reads as coverage and is not" as one of three
# compounding causes, and then round 7 did not remove this one.
#
# It is deleted rather than wired because the check it was for is now done better and
# in the other direction: `tests/test_docs_match_artefacts.py` asserts that every
# occurrence of a derived value in every document has one of its fact's phrasings
# beside it, and that every number carrying such a phrasing equals the derived value.
# That catches any wrong value, not a guessed list of plausible ones, so an absence
# check over invented candidates would add nothing. Standing rule 14 is amended in
# PITFALLS.md to say so.
