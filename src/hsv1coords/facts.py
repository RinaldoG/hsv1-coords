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
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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
        spans = sum(1 for r in rows if r["source"] == "BK012101.1"
                    and r["feature_class"] == "transcript_span")
        facts.append(Fact("transcript_spans", spans,
                          ("of its named features are transcript",
                           "of that record's features are transcript")))

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

    coverage = DERIVED / "crosswalk_coverage.tsv"
    if coverage.exists():
        rows = _rows(coverage)
        facts.append(Fact("genomes", len(rows), ("complete genomes",)))
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
                           "bp. The offset between them")))
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
                           "bp, which is `.1`")))

    try:
        from .conventions import ConventionSet
        cs = ConventionSet.load()
        # The worked example's own numbers, derived from behaviour rather than typed.
        lat = cs["szpara2014_trimmed"].to_convention(118777)
        if lat is not None:
            facts.append(Fact("lat_in_szpara", lat,
                              ("in a Szpara-trimmed analysis",
                               "in a Szpara-trimmed")))
        facts.append(Fact("conventions", len(cs.conventions),
                          ("handlings, all", "handlings, each", "handlings,")))
    except Exception:                                   # artefacts not built
        pass

    return facts


def superseded_values(fact: Fact, span: int = 60) -> set[str]:
    """Plausible wrong values for this figure, for an absence check.

    Only generated for small integers, where an off-by-a-few is the realistic drift.
    Percentages and large counts are checked by presence of the correct value and by
    the absence of any OTHER figure carrying the same noun, which `tests` does.
    """
    if not isinstance(fact.value, int) or fact.value > span:
        return set()
    out = set()
    for n in range(1, span + 1):
        if n == fact.value:
            continue
        out.add(str(n))
        if n in SPELLED:
            out.add(SPELLED[n])
    return out
