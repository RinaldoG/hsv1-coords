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
import pathlib
from dataclasses import dataclass

DEFAULT_DIR = (pathlib.Path(__file__).resolve().parents[2]
               / "data" / "derived" / "crosswalk")
REFERENCE = "NC_001806.2"


@dataclass(frozen=True)
class Block:
    query_start: int   # 0-based
    target_start: int  # 0-based
    length: int


class Crosswalk:
    """The ungapped blocks relating one genome to the reference."""

    def __init__(self, accession: str, blocks: list[Block]) -> None:
        self.accession = accession
        self.blocks = sorted(blocks, key=lambda b: b.query_start)
        self._q_starts = [b.query_start for b in self.blocks]
        self._by_target = sorted(blocks, key=lambda b: b.target_start)
        self._t_starts = [b.target_start for b in self._by_target]

    @classmethod
    def load(cls, accession: str, directory: pathlib.Path | None = None) -> "Crosswalk":
        directory = directory or DEFAULT_DIR
        path = pathlib.Path(directory) / f"{accession}.blocks.tsv"
        if not path.exists():
            raise FileNotFoundError(f"no crosswalk for {accession} at {path}")
        with path.open() as fh:
            blocks = [Block(int(r["query_start"]), int(r["target_start"]),
                            int(r["length"]))
                      for r in csv.DictReader(fh, delimiter="\t")]
        return cls(accession, blocks)

    @staticmethod
    def available(directory: pathlib.Path | None = None) -> list[str]:
        directory = pathlib.Path(directory or DEFAULT_DIR)
        return sorted(p.name.replace(".blocks.tsv", "")
                      for p in directory.glob("*.blocks.tsv"))

    def to_reference(self, position: int) -> int | None:
        """Lift a 1-based position in this genome to the reference, or None."""
        p = position - 1
        i = bisect.bisect_right(self._q_starts, p) - 1
        if i < 0:
            return None
        b = self.blocks[i]
        if p >= b.query_start + b.length:
            return None
        return b.target_start + (p - b.query_start) + 1

    def from_reference(self, position: int) -> int | None:
        """Lift a 1-based reference position into this genome, or None."""
        p = position - 1
        i = bisect.bisect_right(self._t_starts, p) - 1
        if i < 0:
            return None
        b = self._by_target[i]
        if p >= b.target_start + b.length:
            return None
        return b.query_start + (p - b.target_start) + 1

    @property
    def mapped_bases(self) -> int:
        return sum(b.length for b in self.blocks)

    def __repr__(self) -> str:
        return (f"Crosswalk({self.accession!r}, {len(self.blocks)} blocks, "
                f"{self.mapped_bases:,} bp mapped)")


def lift(position: int, source: str, target: str,
         directory: pathlib.Path | None = None) -> int | None:
    """Lift a 1-based position from one genome to another, through the reference.

    Returns None if either leg has no coordinate, which is the honest answer rather
    than an approximation.
    """
    if source == target:
        # Even the identity case must be a real position in a genome we hold, or the
        # caller gets a confident answer about a coordinate nobody has checked.
        if source != REFERENCE:
            cw = Crosswalk.load(source, directory)
            return position if cw.to_reference(position) is not None else None
        return position
    if source == REFERENCE:
        ref_pos: int | None = position
    else:
        ref_pos = Crosswalk.load(source, directory).to_reference(position)
    if ref_pos is None:
        return None
    if target == REFERENCE:
        return ref_pos
    return Crosswalk.load(target, directory).from_reference(ref_pos)
