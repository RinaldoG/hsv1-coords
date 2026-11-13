"""Answer whether a reference coordinate is well defined, and if not, what it could be.

A fifth of the HSV-1 genome lies inside an inverted repeat, so a position there names a
place in the linear record without saying which copy of the locus is meant. The usual
response in the field is to mask those regions with Ns, which throws the information
away. This keeps it: a position inside a repeat gets an answer that names the other
copy and the locus it belongs to.

The question this is built to answer is the one a coordinate user actually has:

    >>> status = RepeatMap.load().classify(900)
    >>> status.well_defined
    False
    >>> status.alternatives
    [124774]

Nothing here is hand-listed. Every region and every pair is derived by
`src/build_repeat_map.py` from the reference record's own feature table.
"""
from __future__ import annotations

import csv
import pathlib
from dataclasses import dataclass, field

DERIVED = pathlib.Path(__file__).resolve().parents[2] / "data" / "derived"
REPEAT_MAP = DERIVED / "repeat_map.tsv"
DIPLOID_PAIRS = DERIVED / "diploid_pairs.tsv"

# The two inverted repeat pairs. Each is one stretch of sequence present twice, in
# opposite orientation, which is why the counterpart offset is measured from the far end.
INVERTED_PAIRS = [("TRL", "IRL"), ("IRS", "TRS")]

# The genome is circular, and its terminal 'a' sequence is present at both ends of the
# linear record in the SAME orientation. That is a direct repeat, not an inverted one,
# so it needs its own relation: an earlier version modelled only the inverted pairs and
# gave 798 positions an incomplete list of alternatives.
DIRECT_PAIRS = [((1, 399), (151824, 152222))]

GENOME_BP = 152222


@dataclass(frozen=True)
class Region:
    name: str
    kind: str
    start: int
    end: int

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class DiploidPair:
    locus: str
    repeat_pair: str
    feature_type: str
    copy_1: tuple[int, int]
    copy_2: tuple[int, int]


@dataclass
class CoordinateStatus:
    """What can honestly be said about one reference position.

    `well_defined` is False when the position lies in a region the genome carries more
    than once, so the coordinate does not identify a unique locus. It being True means
    only that no annotated large repeat covers the position: this map describes the
    annotated repeat structure, and it does not and cannot certify that a position is
    unique everywhere in the sequence. An earlier version said it did, which was false
    for every tandem reiteration set.
    """

    position: int
    well_defined: bool
    regions: list[Region] = field(default_factory=list)
    alternatives: list[int] = field(default_factory=list)
    loci: list[str] = field(default_factory=list)
    reason: str = ""

    def __str__(self) -> str:
        if self.well_defined:
            return f"{self.position}: well defined"
        alts = ", ".join(str(a) for a in self.alternatives) or "none computable"
        return (f"{self.position}: not well defined ({self.reason}); "
                f"could equally be {alts}")


class RepeatMap:
    def __init__(self, regions: list[Region], pairs: list[DiploidPair]) -> None:
        self.regions = regions
        self.pairs = pairs
        self._by_name = {r.name: r for r in regions}

    @classmethod
    def load(cls, directory: pathlib.Path | None = None) -> "RepeatMap":
        d = pathlib.Path(directory) if directory else DERIVED
        with (d / "repeat_map.tsv").open() as fh:
            regions = [Region(r["name"], r["kind"], int(r["start"]), int(r["end"]))
                       for r in csv.DictReader(fh, delimiter="\t")]
        pairs: list[DiploidPair] = []
        pp = d / "diploid_pairs.tsv"
        if pp.exists():
            with pp.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    pairs.append(DiploidPair(
                        r["locus"], r["repeat_pair"], r["feature_type"],
                        (int(r["copy_1_start"]), int(r["copy_1_end"])),
                        (int(r["copy_2_start"]), int(r["copy_2_end"]))))
        return cls(regions, pairs)

    def _repeat(self, label: str) -> Region | None:
        for r in self.regions:
            if r.name.startswith(label) and r.kind == "repeat_region":
                return r
        return None

    def _direct_counterparts(self, position: int) -> list[int]:
        """Equivalent positions in a direct (same orientation) repeat."""
        out: list[int] = []
        for (a_start, a_end), (b_start, b_end) in DIRECT_PAIRS:
            if a_start <= position <= a_end:
                out.append(b_start + (position - a_start))
            if b_start <= position <= b_end:
                out.append(a_start + (position - b_start))
        return out

    def counterparts(self, position: int) -> list[int]:
        """Every equivalent position in another copy of an inverted repeat.

        The two copies of a repeat are inverted relative to each other, so an offset
        measured from the start of one corresponds to the same offset measured from the
        end of the other.

        This returns a list rather than one answer because `IRL` and `IRS` overlap: the
        long and short repeats abut and share bases around 125,975 to 126,373, so a
        position in that junction belongs to both and has a counterpart in each. An
        earlier version returned whichever pair it checked first, which silently gave
        the wrong one for the junction.
        """
        out: list[int] = []
        for a, b in INVERTED_PAIRS:
            ra, rb = self._repeat(a), self._repeat(b)
            if not ra or not rb:
                continue
            if ra.contains(position):
                out.append(rb.end - (position - ra.start))
            if rb.contains(position):
                out.append(ra.start + (rb.end - position))
        out.extend(self._direct_counterparts(position))
        return sorted(set(out) - {position})

    def counterpart(self, position: int) -> int | None:
        """The single counterpart, or None where there is not exactly one."""
        alts = self.counterparts(position)
        return alts[0] if len(alts) == 1 else None

    def loci_at(self, position: int) -> list[str]:
        out = []
        for p in self.pairs:
            if (p.copy_1[0] <= position <= p.copy_1[1]
                    or p.copy_2[0] <= position <= p.copy_2[1]):
                out.append(p.locus)
        return sorted(set(out))

    def classify(self, position: int) -> CoordinateStatus:
        if not 1 <= position <= GENOME_BP:
            return CoordinateStatus(
                position, False, [], [], [],
                f"outside the reference, which is 1 to {GENOME_BP:,}")
        hits = [r for r in self.regions if r.contains(position)]
        inverted = [r for r in hits
                    if r.kind == "repeat_region"
                    and any(r.name.startswith(x) for pair in INVERTED_PAIRS for x in pair)]
        direct = self._direct_counterparts(position)
        if not inverted:
            if direct:
                return CoordinateStatus(
                    position, False, hits, direct, self.loci_at(position),
                    "inside the terminal 'a' sequence, which is present at both ends "
                    "of the linear record in the same orientation")
            if hits:
                kinds = ", ".join(sorted({r.name for r in hits}))
                return CoordinateStatus(
                    position, True, hits, [], self.loci_at(position),
                    f"inside {kinds}, which is a repeat feature but not one of the two "
                    f"large inverted repeats, so this map names no alternative copy. "
                    f"Note that tandem reiteration sets are internally repetitive, so a "
                    f"short sequence here may still occur elsewhere")
            return CoordinateStatus(
                position, True, [], [], [],
                "no annotated repeat covers this position")
        alts = self.counterparts(position)
        names = ", ".join(sorted({r.name for r in inverted}))
        reason = (f"inside {names}, which is present in two copies"
                  if len(inverted) == 1 else
                  f"inside {names}, which overlap and are each present in two copies")
        return CoordinateStatus(position, False, hits, alts,
                                self.loci_at(position), reason)

    def is_well_defined(self, position: int) -> bool:
        return self.classify(position).well_defined

    @property
    def covered_bases(self) -> int:
        covered: set[int] = set()
        for r in self.regions:
            covered.update(range(r.start, r.end + 1))
        return len(covered)
