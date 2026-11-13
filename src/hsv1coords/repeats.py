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
    [125474]

Every region, every pair and every copy's orientation is derived by
`src/build_repeat_map.py` from the reference record's own feature table and sequence,
and read from `data/derived/repeat_map.tsv` here. Round 7 found this paragraph was not
true: the three copies of the 'a' sequence and their orientations were a literal list
in this module, and the inverted pairs were a literal list of names, so a rebuild that
moved either would have left them stale and silent. Both are now read from the
artefact, and `tests/test_m4_repeat_map.py` fails if this module reintroduces either.
"""
from __future__ import annotations

import csv
import pathlib
import re
from dataclasses import dataclass, field

DERIVED = pathlib.Path(__file__).resolve().parents[2] / "data" / "derived"
REPEAT_MAP = DERIVED / "repeat_map.tsv"
DIPLOID_PAIRS = DERIVED / "diploid_pairs.tsv"

GENOME_BP = 152222


@dataclass(frozen=True)
class Region:
    name: str
    kind: str
    start: int
    end: int
    orientation: str = "single"
    notes: str = ""

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
    note: str = ""

    @property
    def exact(self) -> bool:
        """True when both copies have the same length and the pairing is positional.

        A partial pair is one where the same locus is named in both repeat copies but
        the two differ in length, which happens when a copy is truncated by the end of
        the linear record or extends past the repeat. LAT is both: its long-repeat copy
        runs 778 bp beyond IRL, and its other copy is cut off at position 1.
        """
        return not self.note.startswith("partial")


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
        self.inverted_pairs = self._derive_inverted_pairs()
        self.multi_copy = self._derive_multi_copy()

    def _derive_inverted_pairs(self) -> list[tuple[str, str]]:
        """The inverted repeat pairs, from what the record says each region flanks.

        The reference annotates each as "inverted repeat flanking UL" or "... US", so
        the pairing is in the record and does not need naming here. A segment that does
        not have exactly two flanking repeats is skipped rather than guessed at.
        """
        flanking: dict[str, list[str]] = {}
        for r in self.regions:
            m = re.search(r"inverted repeat flanking (\w+)", r.notes or "")
            if m and r.kind == "repeat_region":
                flanking.setdefault(m.group(1), []).append(r.name)
        out = []
        for segment in sorted(flanking):
            names = sorted(set(flanking[segment]),
                           key=lambda n: self._by_name[n].start)
            if len(names) == 2:
                out.append((names[0], names[1]))
        return out

    def _derive_multi_copy(self) -> dict[str, list[Region]]:
        """Regions the record names more than once, keyed by name.

        The 'a' sequence is the live case: three copies, not two, and not all in the
        same orientation. The orientations come from the artefact's own column, which
        `src/build_repeat_map.py` measures against the sequence.
        """
        groups: dict[str, list[Region]] = {}
        for r in self.regions:
            groups.setdefault(r.name, []).append(r)
        return {n: sorted(c, key=lambda r: r.start)
                for n, c in groups.items() if len(c) > 1}

    @classmethod
    def load(cls, directory: pathlib.Path | None = None) -> "RepeatMap":
        d = pathlib.Path(directory) if directory else DERIVED
        with (d / "repeat_map.tsv").open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
            if rows and "orientation" not in rows[0]:
                raise ValueError(
                    f"{d / 'repeat_map.tsv'} has no orientation column. Each copy's "
                    f"orientation is measured against the sequence by "
                    f"`src/build_repeat_map.py`; without it a counterpart cannot be "
                    f"computed. Rebuild with `make repeatmap`.")
            regions = [Region(r["name"], r["kind"], int(r["start"]), int(r["end"]),
                              r.get("orientation", "single"), r.get("notes", ""))
                       for r in rows]
        pairs: list[DiploidPair] = []
        pp = d / "diploid_pairs.tsv"
        if pp.exists():
            with pp.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    pairs.append(DiploidPair(
                        r["locus"], r["repeat_pair"], r["feature_type"],
                        (int(r["copy_1_start"]), int(r["copy_1_end"])),
                        (int(r["copy_2_start"]), int(r["copy_2_end"])),
                        r.get("offset_delta", "")))
        return cls(regions, pairs)

    def _repeat(self, label: str) -> Region | None:
        for r in self.regions:
            if r.name.startswith(label) and r.kind == "repeat_region":
                return r
        return None

    def _multi_copy_counterparts(self, position: int) -> list[int]:
        """Every other copy of a region the record names more than once.

        The live case is the 'a' sequence, which the reference carries three times and
        not all in the same orientation, so a position in one has two counterparts and
        the sign of the offset differs between them. Copies of unequal length, or whose
        orientation could not be measured, are skipped: the honest answer there is no
        counterpart rather than a computed one.
        """
        out: list[int] = []
        for copies in self.multi_copy.values():
            home = next((r for r in copies if r.contains(position)), None)
            if home is None or home.orientation not in ("forward", "reverse"):
                continue
            if len({r.length for r in copies}) != 1:
                continue
            offset = (position - home.start if home.orientation == "forward"
                      else home.end - position)
            for r in copies:
                if r is home:
                    continue
                if r.orientation not in ("forward", "reverse"):
                    continue
                out.append(r.start + offset if r.orientation == "forward"
                           else r.end - offset)
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
        for a, b in self.inverted_pairs:
            ra, rb = self._repeat(a), self._repeat(b)
            if not ra or not rb:
                continue
            if ra.contains(position):
                out.append(rb.end - (position - ra.start))
            if rb.contains(position):
                out.append(ra.start + (rb.end - position))
        out.extend(self._multi_copy_counterparts(position))
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
                    and any(r.name.startswith(x)
                            for pair in self.inverted_pairs for x in pair)]
        a_seq = self._multi_copy_counterparts(position)
        if not inverted:
            if a_seq:
                return CoordinateStatus(
                    position, False, hits, a_seq, self.loci_at(position),
                    "inside the 'a' sequence, which the reference carries three times: "
                    "at both ends of the linear record in the same orientation, and "
                    "once inverted between the segments")
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
