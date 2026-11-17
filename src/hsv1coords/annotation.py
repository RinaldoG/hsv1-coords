"""What is at a reference position. The question every hand check needed and none could ask.

Three independent reviewers resolved published coordinates through this library on 17
November 2026 and all three put the same gap first. Every one of their hand checks
needed the answer to "what feature sits at `NC_001806.2` 112,361", and the library had
no call for it: the annotation existed as a TSV and as GFF3, and reaching it meant
grepping `annotation_v0.1.0.tsv` and parsing `ref_NC_001806.2.gb` by hand, including
writing the minus-strand arithmetic. One of them wrote it as "resolving a coordinate
without this is only half an answer: you get a number and still cannot say what it
means". `codon()` is the half that existed, and it goes the other way, from a gene and
a codon number to positions.

Two things this is careful about, because getting either wrong would make it worse than
nothing.

**An intron is not a feature.** A spliced feature covers its exons, not the gaps
between them. `NC_001806.2` annotates LAT as a gene over `118777..127151` and as an
ncRNA over `join(118805..119464,121421..127151)`, and 120,505 is inside the gene and
inside neither exon, which is precisely how the reviewers established that Amelio et
al.'s CTCCC cluster is in the LAT intron. So `features_at` reports what COVERS the
position, `introns_at` reports features whose outer extent contains it while none of
their spans do, and the two are separate answers rather than one blurred one. A single
list would have made the intron case indistinguishable from the exon case, which is
the distinction the question was being asked for.

**A position in a repeat names two places.** The annotation would happily report the
IRL copy of a feature and say nothing about the TRL copy carrying the same bases, which
is the error this whole repository exists to prevent. `describe` consults the repeat
map and leads with the ambiguity when there is one.

What this does not do: it does not translate, it does not lift, and it answers about
the reference only. A coordinate from another strain is lifted first, and a coordinate
in `NC_001806.1` is translated first; both refusals name the call to make.
"""
from __future__ import annotations

import csv
import functools
import re
from dataclasses import dataclass

from ._coordinates import as_coordinate
from ._paths import data_path, require

ANNOTATION = data_path("derived", "annotation_v0.1.0.tsv")
REFERENCE = "NC_001806.2"
GENOME_BP = 152222
SPAN = re.compile(r"(\d+)\.\.[<>]?(\d+)")


@dataclass(frozen=True)
class AnnotatedFeature:
    """One annotated feature of one record, with its spans already parsed."""

    name: str
    feature_type: str
    feature_class: str
    location: str
    source: str
    product: str
    notes: str
    provenance: str
    spans: tuple[tuple[int, int], ...]
    strand: str
    """"forward" or "complement", from the location's own `complement(...)`.

    Round 10 shipped 423 minus-strand features as plain ascending ranges because a
    parser never read the strand column, so the absence of `complement` read as a
    positive assertion of the plus strand for 83 genes. The strand is taken from the
    location text here, which is the only place this table states it.
    """

    @property
    def extent(self) -> tuple[int, int]:
        """The outermost bounds, which for a spliced feature include its introns."""
        return (min(a for a, _ in self.spans), max(b for _, b in self.spans))

    def covers(self, position: int) -> bool:
        """True when the position is inside a span, so inside an exon."""
        return any(a <= position <= b for a, b in self.spans)

    def within_extent(self, position: int) -> bool:
        start, end = self.extent
        return start <= position <= end

    @property
    def spliced(self) -> bool:
        return len(self.spans) > 1

    def __str__(self) -> str:
        name = self.name or "(unnamed)"
        product = f", {self.product}" if self.product else ""
        strand = " on the complement strand" if self.strand == "complement" else ""
        return (f"{self.feature_type} {name}{product} at {self.location}{strand} "
                f"[{self.source}]")


def _parse(location: str) -> tuple[tuple[tuple[int, int], ...], str]:
    strand = "complement" if "complement(" in location else "forward"
    spans = tuple((int(a), int(b)) for a, b in SPAN.findall(location))
    return spans, strand


@functools.lru_cache(maxsize=1)
def _rows() -> tuple[AnnotatedFeature, ...]:
    out = []
    with require(ANNOTATION).open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            spans, strand = _parse(row["location"])
            if not spans:
                # A single-base location such as `1234` states no span. None exist in
                # the shipped table; a row that gained one would otherwise be dropped
                # silently, which is root cause 6, filtering on a derived field.
                continue
            out.append(AnnotatedFeature(
                name=row["name"], feature_type=row["feature_type"],
                feature_class=row["feature_class"], location=row["location"],
                source=row["source"], product=row.get("product", ""),
                notes=row.get("notes", ""), provenance=row.get("provenance", ""),
                spans=spans, strand=strand))
    if not out:
        raise ValueError(
            f"{ANNOTATION} states no features, so every question about a position "
            f"would be answered with an empty list, which reads like 'nothing is "
            f"there'. Run `make annotation`.")
    return tuple(out)


class NotAReferencePosition(ValueError):
    """A position that cannot be a reference position, or a frame that is not this one."""


def _check(position: object) -> int:
    pos = as_coordinate(position, "a reference position")
    if not 1 <= pos <= GENOME_BP:
        raise NotAReferencePosition(
            f"{pos} is outside {REFERENCE}, which is 1 to {GENOME_BP}. A coordinate "
            f"from another strain is lifted first:\n"
            f"    from hsv1coords import lift\n"
            f"    lift(position, 'JQ673480.1', '{REFERENCE}')\n"
            f"and a coordinate in NC_001806.1, which is 152,261 bp, is translated "
            f"first:\n"
            f"    from hsv1coords.versions import VersionCrosswalk\n"
            f"    VersionCrosswalk.load().convert(position, 'NC_001806.1', "
            f"'{REFERENCE}')")
    return pos


def sources() -> tuple[str, ...]:
    """Every record the annotation table describes.

    Three of them describe the identical 152,222 bases under different names and at
    different scopes: 77 CDS in `NC_001806.2` against about 300 in `BK012101.1`. So a
    question about a position has to say which annotation it is asking about, or get
    all of them and be told which is which.
    """
    return tuple(sorted({f.source for f in _rows()}))


def features_at(position: int, source: str | None = REFERENCE
                ) -> tuple[AnnotatedFeature, ...]:
    """Features whose spans cover this reference position.

    `source=None` asks every annotation rather than the reference's own, which is
    usually what a reader wants second: the 2020 re-annotation names things the 2018
    record does not.

    A spliced feature is reported only when the position is in one of its exons.
    `introns_at` is the other half of that answer, deliberately separate.

    >>> [f.name for f in features_at(112361) if f.feature_type == "CDS"]
    ['UL53']
    """
    pos = _check(position)
    return tuple(f for f in _rows()
                 if (source is None or f.source == source) and f.covers(pos))


def introns_at(position: int, source: str | None = REFERENCE
               ) -> tuple[AnnotatedFeature, ...]:
    """Spliced features that contain this position between their exons, not in one.

    This is a different answer from `features_at`, not a longer one. A caller slicing
    the reference between a spliced feature's outer bounds picks up intron bases
    without noticing, and a caller told only "LAT" about position 120,505 would
    reasonably conclude the position is transcribed, which it is not.
    """
    pos = _check(position)
    return tuple(f for f in _rows()
                 if (source is None or f.source == source)
                 and f.spliced and f.within_extent(pos) and not f.covers(pos))


def describe(position: int, source: str | None = REFERENCE) -> str:
    """What is at this position, in the form a caller can print.

    Leads with the repeat ambiguity when there is one, because a position that names
    two places has no single answer to this question either, and a list of features
    that did not say so would be the most confident wrong answer this library could
    give.
    """
    pos = _check(position)
    lines = [f"{REFERENCE}:{pos:,}"]
    try:
        from .repeats import RepeatMap
        status = RepeatMap.load().classify(pos)
    except Exception:                                   # artefacts not built
        status = None
    if status is not None and not status.well_defined:
        alts = ", ".join(f"{a:,}" for a in status.alternatives) or "none computable"
        lines.append(f"  NOT well defined: this position is inside a repeat present "
                     f"in more than one copy, and the same bases are at {alts}. "
                     f"Whatever is listed below is listed for this copy only.")
    covering = features_at(pos, source)
    if covering:
        lines.append("  covered by:")
        lines += [f"    {f}" for f in covering]
    else:
        lines.append("  covered by: nothing in "
                     + (source or "any annotation here"))
    inside = introns_at(pos, source)
    if inside:
        lines.append("  between the exons of (so NOT transcribed as part of):")
        lines += [f"    {f}" for f in inside]
    return "\n".join(lines)
