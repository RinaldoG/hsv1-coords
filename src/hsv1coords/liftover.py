"""Lift a coordinate between an HSV-1 genome and the reference, or decline.

Every genome in the crosswalk is aligned to `NC_001806.2` and reduced to ungapped blocks.
A lift is a binary search over those blocks. Two genomes are related through the reference.

The contract is that this **declines** rather than guesses. A position in an insertion
relative to the reference has no reference coordinate, and `lift` returns `None` for it.
A position that lands in a region the alignment could not place uniquely likewise returns
`None`. Returning a nearby coordinate instead would be the silent error this whole
repository exists to prevent, and it is exactly how a position inside the inverted repeats
ends up pointing at the wrong copy of LAT.

Coordinates are 1-based and inclusive on the public interface, matching how papers and
GenBank cite positions. The stored blocks are 0-based half-open, as the aligner emits them.
"""
from __future__ import annotations

import bisect
import csv
import functools
import pathlib
from dataclasses import dataclass, field

DEFAULT_DIR = (pathlib.Path(__file__).resolve().parents[2]
               / "data" / "derived" / "crosswalk")
REFERENCE = "NC_001806.2"


@dataclass(frozen=True)
class Block:
    query_start: int   # 0-based
    target_start: int  # 0-based
    length: int


class Crosswalk:
    """The ungapped blocks relating one genome to the reference.

    A genome deposited on the opposite strand from the reference is reverse
    complemented before alignment, so its block offsets are positions in the reverse
    complement, not in the deposited record. That frame is declared on the first line
    of the block file and applied here, at the one boundary where a deposited
    coordinate enters. Round 7 found it declared in a different file and applied
    nowhere, which made every lift against those three genomes wrong by a median of 30
    to 50 kb while looking perfectly ordinary.
    """

    def __init__(self, accession: str, blocks: list[Block],
                 orientation: str = "forward", deposited_length: int | None = None,
                 reverse_strand_bp: int = 0) -> None:
        if orientation not in ("forward", "reverse_complemented"):
            raise ValueError(f"unknown orientation {orientation!r} for {accession}")
        if orientation == "reverse_complemented" and not deposited_length:
            raise ValueError(
                f"{accession} is reverse complemented but its block file states no "
                f"deposited_length, so a deposited coordinate cannot be placed")
        self.accession = accession
        self.orientation = orientation
        self.deposited_length = deposited_length
        self.reverse_strand_bp = reverse_strand_bp
        """How much of this genome aligned on the other strand and was dropped.

        For most genomes this is a few hundred base pairs of repeat artefact. For a
        genome deposited as a segmental isomer it is a whole segment: `ON783215.1` has
        its long segment forward and its short segment reversed, so about 30 kb has no
        reference coordinate. Round 9: a lift there returned None with nothing saying
        why, and nothing in the documents told a reader a whole segment was missing.
        """
        # Tuples, not lists, because `load` is cached and hands the same object to
        # every caller. Round 9: one caller doing `cw.blocks.clear()` changed what every
        # later `load` returned, taking `to_reference(60000)` from 61474 to None for a
        # caller that did nothing wrong. The docstrings said the objects were "treated
        # as immutable within a run", which is a request, not a guarantee.
        self.blocks = tuple(sorted(blocks, key=lambda b: b.query_start))
        self._q_starts = [b.query_start for b in self.blocks]
        self._by_target = tuple(sorted(blocks, key=lambda b: b.target_start))
        self._t_starts = [b.target_start for b in self._by_target]

    @classmethod
    @functools.lru_cache(maxsize=256)   # > 183 genomes: at 64 a corpus-wide
    def load(cls, accession: str, directory: pathlib.Path | None = None) -> "Crosswalk":
        """Cached: a genome-to-genome lift loads two crosswalks, twice over.

        Round 8 counted four `load` calls per cross-genome lift, each re-reading and
        re-sorting up to 1088 blocks.

        Round 9 measured the cache thrashing at `maxsize=64` against 183 genomes: a
        sweep over the corpus gave 6 hits and 185 misses, so the speed-up held only for
        a working set smaller than the corpus. The blocks are tuples, so the shared
        object cannot be mutated by one caller and seen by another.
        """
        directory = directory or DEFAULT_DIR
        path = pathlib.Path(directory) / f"{accession}.blocks.tsv"
        if not path.exists():
            raise FileNotFoundError(f"no crosswalk for {accession} at {path}")
        lines = path.read_text().splitlines(keepends=True)
        if not lines or not lines[0].startswith("#"):
            raise ValueError(
                f"{path} states no frame on its first line. A block file without a "
                f"frame cannot be read safely: a reverse-complemented genome would be "
                f"silently treated as forward. Rebuild with `make crosswalk`.")
        meta = dict(kv.split("=", 1) for kv in lines[0].lstrip("#").split())
        with path.open() as fh:
            next(fh)
            blocks = [Block(int(r["query_start"]), int(r["target_start"]),
                            int(r["length"]))
                      for r in csv.DictReader(fh, delimiter="\t")]
        length = meta.get("deposited_length")
        return cls(accession, blocks, meta.get("orientation", ""),
                   int(length) if length else None,
                   int(meta.get("reverse_strand_bp") or 0))

    def _to_aligned(self, position: int) -> int:
        """A 1-based deposited position as a 1-based position in the aligned frame."""
        if self.orientation == "forward":
            return position
        return self.deposited_length + 1 - position

    def _to_deposited(self, position: int) -> int:
        """The inverse. The reverse complement is its own inverse, hence one formula."""
        return self._to_aligned(position)

    @staticmethod
    def available(directory: pathlib.Path | None = None) -> list[str]:
        directory = pathlib.Path(directory or DEFAULT_DIR)
        return sorted(p.name.replace(".blocks.tsv", "")
                      for p in directory.glob("*.blocks.tsv"))

    def to_reference(self, position: int) -> int | None:
        """Lift a 1-based deposited position in this genome to the reference, or None."""
        if not isinstance(position, int) or isinstance(position, bool):
            raise TypeError(f"a coordinate must be an int, got {position!r}")
        if position < 1:
            return None
        if self.deposited_length and position > self.deposited_length:
            return None
        p = self._to_aligned(position) - 1
        i = bisect.bisect_right(self._q_starts, p) - 1
        if i < 0:
            return None
        b = self.blocks[i]
        if p >= b.query_start + b.length:
            return None
        return b.target_start + (p - b.query_start) + 1

    def from_reference(self, position: int) -> int | None:
        """Lift a 1-based reference position into this genome's deposited frame."""
        if not isinstance(position, int) or isinstance(position, bool):
            raise TypeError(f"a coordinate must be an int, got {position!r}")
        if position < 1:
            return None
        p = position - 1
        i = bisect.bisect_right(self._t_starts, p) - 1
        if i < 0:
            return None
        b = self._by_target[i]
        if p >= b.target_start + b.length:
            return None
        return self._to_deposited(b.query_start + (p - b.target_start) + 1)

    @property
    def mapped_bases(self) -> int:
        return sum(b.length for b in self.blocks)

    def __repr__(self) -> str:
        return (f"Crosswalk({self.accession!r}, {len(self.blocks)} blocks, "
                f"{self.mapped_bases:,} bp mapped)")


@dataclass
class LiftResult:
    """A lifted coordinate together with what can honestly be said about it."""

    position: int | None
    source: str
    target: str
    reference_position: int | None = None
    well_defined: bool = True
    alternatives: list[int] = field(default_factory=list)
    loci: list[str] = field(default_factory=list)
    reason: str = ""
    complemented: bool = False
    """True when the two coordinates are on opposite strands.

    Three of the 183 deposited genomes are on the opposite strand from the reference, so
    a lift into or out of one of them is positionally correct while the base at the
    answer is the complement of the base at the question. Round 8: the frame handling was
    right and nothing told the user, so someone lifting variant positions out of
    `ON783216.1` would have got reference positions whose alleles were the complement of
    theirs, with no warning. That is the shape of error this repository exists to
    prevent, so it is a field and it is in `__str__`.
    """

    inverted_alternatives: tuple[int, ...] = ()
    """Which of `alternatives` are on the opposite strand from this answer.

    Round 9 found the alternatives printed bare: "could equally be 7597", where 7597 is
    the inverted copy and its allele is the complement. Fed from
    `CoordinateStatus.inverted_alternatives`, which is verified base by base against the
    reference sequence in `tests/test_m4_repeat_map.py`. Last in the dataclass, as
    standing rule 22 requires.
    """

    def __bool__(self) -> bool:
        return self.position is not None

    def __str__(self) -> str:
        strand = (", ON THE OPPOSITE STRAND, so the base here is the complement of "
                  "the base there" if self.complemented else "")
        if self.position is None:
            # The reason, which this dropped. Round 9 added a decline that names a
            # missing segment and it did not reach the printed form, which is the same
            # failure as the strand note in round 8: a correct value nobody is told.
            because = f": {self.reason}" if self.reason else ""
            return (f"{self.source}:{self.position} does not lift to {self.target}"
                    f"{because}")
        if self.well_defined:
            return f"{self.source} -> {self.target}: {self.position}{strand}"
        alts = ", ".join(
            f"{a} (opposite strand, so complement the allele)"
            if a in self.inverted_alternatives else str(a)
            for a in self.alternatives) or "none computable"
        loci = f" ({', '.join(self.loci)})" if self.loci else ""
        return (f"{self.source} -> {self.target}: {self.position}{strand}, but NOT "
                f"well defined{loci}: {self.reason}. The reference position could "
                f"equally be {alts}")


def lift_detailed(position: int, source: str, target: str,
                  directory: pathlib.Path | None = None) -> LiftResult:
    """Lift a coordinate and consult the repeat map about the result.

    `lift` gives a number. This gives the number plus whether it identifies a unique
    locus, which it does not when the position falls inside one of the inverted
    repeats: roughly a fifth of this genome. Reporting such a coordinate as though it
    were unambiguous is the error the repeat map exists to prevent, so anything acting
    on a coordinate should use this rather than `lift`.
    """
    from .repeats import RepeatMap

    def orientation_of(accession: str) -> str:
        if accession == REFERENCE:
            return "forward"
        return Crosswalk.load(accession, directory).orientation

    if source == REFERENCE:
        ref_pos: int | None = position
    else:
        ref_pos = Crosswalk.load(source, directory).to_reference(position)
    # Opposite strands iff exactly one of the two is a reverse-complemented deposit.
    complemented = ((orientation_of(source) == "reverse_complemented")
                    != (orientation_of(target) == "reverse_complemented"))
    if ref_pos is None:
        # No answer, so no strand: `complemented` described a comparison between two
        # coordinates and one of them does not exist. Round 9 found it set to True on a
        # declined lift, which reads as a fact about an answer that was never given.
        return LiftResult(None, source, target, None, True, [], [],
                          _why_no_coordinate(source, target, position, directory),
                          complemented=False)

    out = position if source == target else (
        ref_pos if target == REFERENCE
        else Crosswalk.load(target, directory).from_reference(ref_pos))

    status = RepeatMap.load().classify(ref_pos)
    return LiftResult(out, source, target, ref_pos, status.well_defined,
                      status.alternatives, status.loci, status.reason,
                      complemented=complemented,
                      inverted_alternatives=status.inverted_alternatives)


SEGMENT_BP = 5000
"""Above this much dropped reverse-strand alignment, a genome is a segmental isomer.

A few hundred base pairs on the other strand is repeat artefact. Thirty thousand is a
segment. The threshold only decides the wording of a decline, never a coordinate.
"""


def _why_no_coordinate(source: str, target: str, position: int,
                       directory: pathlib.Path | None) -> str:
    """The reason a lift declined, naming a missing segment where that is the cause.

    Round 9: this said "no reference coordinate for this position" whether the position
    sat in a small insertion or inside the 30 kb of `ON783215.1` that the crosswalk does
    not represent at all, because that genome's short segment is deposited inverted and
    alignments on the other strand are dropped. A user lifting a variant in gD or ICP4
    out of it got None with no way to tell the difference.
    """
    plain = "no reference coordinate for this position"
    if source == REFERENCE:
        return plain
    try:
        cw = Crosswalk.load(source, directory)
    except FileNotFoundError:          # pragma: no cover
        return plain
    if cw.reverse_strand_bp < SEGMENT_BP:
        return plain
    return (f"no reference coordinate for this position. {source} has "
            f"{cw.reverse_strand_bp:,} bp that aligned to the reference on the opposite "
            f"strand, which this crosswalk does not represent: it is deposited as a "
            f"segmental isomer, with one segment inverted relative to the reference. "
            f"Positions in that segment have no coordinate here at all, rather than "
            f"having one that could not be placed")


def lift(position: int, source: str, target: str,
         directory: pathlib.Path | None = None) -> int | None:
    """Lift a 1-based position from one genome to another, through the reference.

    Returns None where there is no single honest answer. That is three cases: either
    leg has no coordinate, the position is not in the genome at all, or the reference
    position sits in a region present in more than one copy, where the alignment's
    single answer is one of two or three equally valid ones.

    The third case is the one round 7 found open. The aligner runs with
    `--secondary=no`, so it reports one placement inside the inverted repeats and the
    ambiguity filter never sees the other: roughly 92 per cent of positions in those
    repeats were receiving a confident answer while the documents said the lift-over
    declined there. The repeat map already knows the counterparts, so this asks it.

    Use `lift_detailed` to get the number together with its alternatives rather than a
    refusal.
    """
    detailed = lift_detailed(position, source, target, directory)
    if detailed.position is None or not detailed.well_defined:
        return None
    return detailed.position
