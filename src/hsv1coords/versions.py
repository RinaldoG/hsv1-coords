"""Convert a coordinate between the two versions of the reference accession.

`NC_001806.1` is 152,261 bp and `NC_001806.2` is 152,222 bp. The difference is not a
constant shift: across the 22 aligned blocks there are **nine distinct offsets**,
ranging from -39 to +5. So a coordinate quoted against one version and read against the
other is wrong by whichever offset applies at that position, and neither record says
which version a paper used unless the paper does.

Published work exists against both. Morse et al. 2017 states 152,261 bp, which is `.1`.
Everything else in this repository is `.2`.
"""
from __future__ import annotations

import bisect
import csv
import pathlib
from dataclasses import dataclass

DEFAULT = (pathlib.Path(__file__).resolve().parents[2]
           / "data" / "derived" / "version_crosswalk.tsv")
OLD, NEW = "NC_001806.1", "NC_001806.2"


@dataclass(frozen=True)
class VersionBlock:
    old_start: int     # 0-based in NC_001806.1
    new_start: int     # 0-based in NC_001806.2
    length: int

    @property
    def offset(self) -> int:
        return self.new_start - self.old_start


class VersionCrosswalk:
    """Ungapped blocks relating NC_001806.1 to NC_001806.2."""

    def __init__(self, blocks: list[VersionBlock]) -> None:
        self.blocks = sorted(blocks, key=lambda b: b.old_start)
        self._old = [b.old_start for b in self.blocks]
        self._by_new = sorted(blocks, key=lambda b: b.new_start)
        self._new = [b.new_start for b in self._by_new]

    @classmethod
    def load(cls, path: pathlib.Path | None = None) -> "VersionCrosswalk":
        with (pathlib.Path(path) if path else DEFAULT).open() as fh:
            return cls([VersionBlock(int(r["query_start"]), int(r["target_start"]),
                                     int(r["length"]))
                        for r in csv.DictReader(fh, delimiter="\t")])

    def to_new(self, position: int) -> int | None:
        """A 1-based NC_001806.1 position, in NC_001806.2, or None."""
        p = position - 1
        i = bisect.bisect_right(self._old, p) - 1
        if i < 0:
            return None
        b = self.blocks[i]
        if p >= b.old_start + b.length:
            return None
        return b.new_start + (p - b.old_start) + 1

    def to_old(self, position: int) -> int | None:
        p = position - 1
        i = bisect.bisect_right(self._new, p) - 1
        if i < 0:
            return None
        b = self._by_new[i]
        if p >= b.new_start + b.length:
            return None
        return b.old_start + (p - b.new_start) + 1

    def offset_at(self, position: int) -> int | None:
        """The offset applying at this NC_001806.1 position, or None if unmapped."""
        new = self.to_new(position)
        return None if new is None else new - position

    @property
    def distinct_offsets(self) -> list[int]:
        return sorted({b.offset for b in self.blocks})

    def convert(self, position: int, source: str, target: str) -> int | None:
        if source == target:
            return position
        if source == OLD and target == NEW:
            return self.to_new(position)
        if source == NEW and target == OLD:
            return self.to_old(position)
        raise ValueError(f"unknown version pair {source} to {target}")
