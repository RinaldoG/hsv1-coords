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
                 orientation: str = "forward", deposited_length: int | None = None
                 ) -> None:
        if orientation not in ("forward", "reverse_complemented"):
            raise ValueError(f"unknown orientation {orientation!r} for {accession}")
        if orientation == "reverse_complemented" and not deposited_length:
            raise ValueError(
                f"{accession} is reverse complemented but its block file states no "
                f"deposited_length, so a deposited coordinate cannot be placed")
        self.accession = accession
        self.orientation = orientation
        self.deposited_length = deposited_length
        self.blocks = sorted(blocks, key=lambda b: b.query_start)
        self._q_starts = [b.query_start for b in self.blocks]
        self._by_target = sorted(blocks, key=lambda b: b.target_start)
        self._t_starts = [b.target_start for b in self._by_target]

    @classmethod
    @functools.lru_cache(maxsize=64)
    def load(cls, accession: str, directory: pathlib.Path | None = None) -> "Crosswalk":
        """Cached: a genome-to-genome lift loads two crosswalks, twice over.

        Round 8 counted four `load` calls per cross-genome lift, each re-reading and
        re-sorting up to 1088 blocks. Same immutability assumption as `RepeatMap.load`.
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
                   int(length) if length else None)

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

    def __bool__(self) -> bool:
        return self.position is not None

    def __str__(self) -> str:
        strand = (", ON THE OPPOSITE STRAND, so the base here is the complement of "
                  "the base there" if self.complemented else "")
        if self.position is None:
            return f"{self.source}:{self.position} does not lift to {self.target}"
        if self.well_defined:
            return f"{self.source} -> {self.target}: {self.position}{strand}"
        alts = ", ".join(str(a) for a in self.alternatives) or "none computable"
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
        return LiftResult(None, source, target, None, True, [], [],
                          "no reference coordinate for this position",
                          complemented=complemented)

    out = position if source == target else (
        ref_pos if target == REFERENCE
        else Crosswalk.load(target, directory).from_reference(ref_pos))

    status = RepeatMap.load().classify(ref_pos)
    return LiftResult(out, source, target, ref_pos, status.well_defined,
                      status.alternatives, status.loci, status.reason,
                      complemented=complemented)


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
