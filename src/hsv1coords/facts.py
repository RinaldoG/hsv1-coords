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
                          ("inputs fetched and checksummed", "inputs, checksummed",
                           "inputs and the skeleton"),
                          "rows in data/raw/MANIFEST.tsv"))

    ann = DERIVED / "annotation_v0.1.0.tsv"
    if ann.exists():
        rows = _rows(ann)
        facts.append(Fact("annotation_rows", len(rows),
                          ("feature rows", "rows across three sources")))
        spans = sum(1 for r in rows if r["source"] == "BK012101.1"
                    and r["feature_class"] == "transcript_span")
        facts.append(Fact("transcript_spans", spans,
                          ("of its named features are transcript",
                           "features are multi-gene transcript spans",
                           "of the 2020 record's features are transcript")))

    aliases = DERIVED / "aliases_literature.tsv"
    if aliases.exists():
        rows = _rows(aliases)
        facts.append(Fact("alias_rows", len(rows),
                          ("cited aliases", "cited rows", "alias rows")))
        facts.append(Fact("alias_genes",
                          len({r["canonical_systematic_name"] for r in rows}),
                          ("genes, 103 citations", "genes in the dictionary")))

    repeat_map = DERIVED / "repeat_map.tsv"
    pairs_file = DERIVED / "diploid_pairs.tsv"
    if repeat_map.exists() and pairs_file.exists():
        regions = _rows(repeat_map)
        pairs = _rows(pairs_file)
        partial = sum(1 for p in pairs
                      if str(p.get("offset_delta", "")).startswith("partial"))
        facts.append(Fact("regions", len(regions),
                          ("regions covering", "repeated regions,")))
        facts.append(Fact("diploid_pairs", len(pairs),
                          ("diploid locus pairs", "diploid pairs")))
        facts.append(Fact("partial_pairs", partial,
                          ("of them marked partial", "partial pairs",
                           "of them are marked partial")))

    version = DERIVED / "version_crosswalk.tsv"
    if version.exists():
        blocks = _rows(version)
        offsets = {int(b["target_start"]) - int(b["query_start"]) for b in blocks}
        facts.append(Fact("version_blocks", len(blocks),
                          ("aligned blocks", "blocks, from", "blocks there are")))
        facts.append(Fact("version_offsets", len(offsets), ("distinct offsets",)))

    coverage = DERIVED / "crosswalk_coverage.tsv"
    if coverage.exists():
        rows = _rows(coverage)
        facts.append(Fact("genomes", len(rows),
                          ("complete genomes", "complete HSV-1 genomes")))
        facts.append(Fact(
            "median_mappable",
            f"{statistics.median(float(r['mapped_percent']) for r in rows):.2f}",
            ("per cent, after review", "Median mappable fraction",
             "median mappable fraction"),
            "median mappable fraction"))
        rc = [r for r in rows if r.get("orientation") == "reverse_complemented"]
        facts.append(Fact("reverse_complemented", len(rc),
                          ("genomes deposited on the opposite strand",
                           "were previously at 0 per cent")))
        for r in rc:
            # A distinctive noun per accession, so a percentage cannot match any
            # other percentage in the documents.
            facts.append(Fact(f"revcomp_{r['accession']}", r["mapped_percent"],
                              (f"per cent for {r['accession']}",),
                              f"{r['accession']} mappable fraction"))

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
        facts.append(Fact("genome_bp", len(seq),
                          ("bp, and", "bases of DNA", "bp of the reference",
                           "bp; NC_001806.1", "identical 152,222 bases")))
        repeats = {}
        for loc, quals in _re.findall(
                r"^     repeat_region\s+(\S+)\n((?:\s+/.*\n)+)", body, _re.M):
            note = _re.search(r'/note="([^"]+)"', quals)
            m = _re.match(r"(\d+)\.\.(\d+)$", loc)
            if note and m and "inverted repeat flanking" in note.group(1):
                repeats[note.group(1).split(";")[0]] = (int(m.group(1)),
                                                        int(m.group(2)))
        for name, (a, b) in repeats.items():
            facts.append(Fact(f"{name}_length", b - a + 1,
                              (f"bp each for {name}",) if False else (),
                              f"{name} length"))
        if "TRL" in repeats and "IRS" in repeats:
            facts.append(Fact("long_repeat_bp",
                              repeats["TRL"][1] - repeats["TRL"][0] + 1,
                              ("bp each, identical", "bp copies",
                               "bp of sequence in two copies",
                               "bases, as are")))
            facts.append(Fact("short_repeat_bp",
                              repeats["IRS"][1] - repeats["IRS"][0] + 1,
                              ("bp each. A position", "across 6,633")))

    v1 = RAW / "ref_NC_001806.1.gb"
    if v1.exists():
        import re as _re2
        b1 = v1.read_text()
        s1 = "".join(_re2.findall(r"[acgtnACGTN]",
                                  b1.split("ORIGIN")[1].split("//")[0]
                                  .replace("\n", "")))
        facts.append(Fact("genome_v1_bp", len(s1),
                          ("bp and NC_001806.2", "bp, which is NC_001806.1",
                           "bp against")))

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
