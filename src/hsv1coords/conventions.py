"""Reconcile HSV-1 coordinates between the conventions different groups actually use.

Knowing that a region is repeated is available from several places. The thing that is
not available anywhere, and the reason a printed HSV-1 coordinate cannot be taken at
face value, is that groups handle those regions differently. One excludes them from
analysis by hard-coded coordinates. One masks them during mapping and restores them
afterwards. One trims the genome to a single copy of each repeat, which moves every
coordinate downstream. The same printed number means different things in each.

This holds each convention as data, with its definition quoted from its own paper, and
answers two separate questions about a position:

    is it usable under that convention?   -> `usable`
    what is the equivalent position?      -> `convert`

Those are deliberately separate. A convention that only drops data leaves coordinates
unchanged, so "unusable here" and "moves to a different number" are different answers,
and conflating them is the error this resolver exists to prevent.

A convention whose definition could not be read from its own source is refused, and the
refusal names what was tried. That rule stays because this project has twice encoded
something from a second-hand description that turned out to be wrong. It is not an
excuse to stop looking: `szpara2014_trimmed` was marked unverified on the strength of
one failed API call, and review found the paper free to read on the publisher's own
site. A failed request is not a paywall.

**The two collapsed references disagree with each other, and that is the real cost.**
Szpara trims to one copy of each inverted repeat, keeping `IRL` and `IRS` and removing
`TRL` and `TRS`. Morse removes `IRL` and `TRS` instead. Both are 1-based coordinates on
a shortened genome, both are published without an offset table, and in the region
between them the same printed number means two different places.
"""
from __future__ import annotations

import csv
import pathlib
from dataclasses import dataclass, field

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "conventions.tsv"
GENOME_BP = 152222


class UnverifiedConvention(RuntimeError):
    """Raised when asked to apply a convention whose definition was never read."""


@dataclass(frozen=True)
class Convention:
    name: str
    kind: str                  # none, exclude, mask_then_restore, trim
    intervals: list[tuple[int, int]]
    reference: str
    verified: bool
    source: str
    quote: str

    @property
    def changes_coordinates(self) -> bool:
        """Only removing sequence moves the numbers. Masking and excluding do not."""
        return self.kind == "trim"

    @property
    def affected_bases(self) -> int:
        return sum(b - a + 1 for a, b in self.intervals)

    def retained_blocks(self, genome_bp: int = GENOME_BP) -> list[tuple[int, int]]:
        """The reference intervals this convention keeps, in order."""
        if not self.changes_coordinates:
            return [(1, genome_bp)]
        blocks: list[tuple[int, int]] = []
        cursor = 1
        for a, b in sorted(self.intervals):
            if a > cursor:
                blocks.append((cursor, a - 1))
            cursor = max(cursor, b + 1)
        if cursor <= genome_bp:
            blocks.append((cursor, genome_bp))
        return blocks

    def length(self, genome_bp: int = GENOME_BP) -> int:
        return sum(b - a + 1 for a, b in self.retained_blocks(genome_bp))

    def to_convention(self, reference_position: int,
                      genome_bp: int = GENOME_BP) -> int | None:
        """Reference coordinate to this convention's coordinate, or None if removed."""
        if not self.changes_coordinates:
            return reference_position
        offset = 0
        for a, b in self.retained_blocks(genome_bp):
            if a <= reference_position <= b:
                return offset + (reference_position - a) + 1
            offset += b - a + 1
        return None

    def to_reference(self, position: int, genome_bp: int = GENOME_BP) -> int | None:
        """This convention's coordinate back to the reference, or None if out of range."""
        if not self.changes_coordinates:
            return position
        remaining = position
        for a, b in self.retained_blocks(genome_bp):
            size = b - a + 1
            if remaining <= size:
                return a + remaining - 1
            remaining -= size
        return None


@dataclass
class PositionReport:
    """What one reference position means under every convention."""

    position: int
    usable_under: list[str] = field(default_factory=list)
    unusable_under: list[tuple[str, str]] = field(default_factory=list)
    refused: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)

    @property
    def comparable_everywhere(self) -> bool:
        """Usable under every convention this resolver can evaluate.

        Deliberately ignores the refused ones: counting them would make this False at
        every position in the genome whenever any convention is unverified, which tells
        a reader nothing. `refused` is reported separately so it is still visible.
        """
        return not self.unusable_under

    def __str__(self) -> str:
        lines = [f"reference position {self.position}:"]
        for name in self.usable_under:
            lines.append(f"  usable under {name}")
        for name, why in self.unusable_under:
            lines.append(f"  NOT usable under {name}: {why}")
        for name in self.refused:
            lines.append(f"  cannot say for {name}: its definition is unverified")
        for c in self.caveats:
            lines.append(f"  caveat: {c}")
        return "\n".join(lines)


def _parse_intervals(text: str) -> list[tuple[int, int]]:
    if not text.strip() or text.strip() == "UNVERIFIED":
        return []
    out = []
    for part in text.split(";"):
        a, b = part.split("-")
        out.append((int(a), int(b)))
    return out


class ConventionSet:
    def __init__(self, conventions: list[Convention]) -> None:
        self.conventions = conventions
        self._by_name = {c.name: c for c in conventions}

    @classmethod
    def load(cls, path: pathlib.Path | None = None) -> "ConventionSet":
        with (pathlib.Path(path) if path else DATA).open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return cls([Convention(
            name=r["convention"],
            kind=r["kind"],
            intervals=_parse_intervals(r["intervals"]),
            reference=r["reference"],
            verified=r["verified"].strip().lower() == "yes",
            source=r["source"],
            quote=r["quote"],
        ) for r in rows])

    def __getitem__(self, name: str) -> Convention:
        return self._by_name[name]

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.conventions]

    @property
    def verified(self) -> list[Convention]:
        return [c for c in self.conventions if c.verified]

    def _require_verified(self, c: Convention) -> None:
        if not c.verified:
            raise UnverifiedConvention(
                f"{c.name} was never verified against its source, so this resolver "
                f"will not apply it. {c.source}")

    def usable(self, position: int, convention: str) -> tuple[bool, str]:
        """Is this reference position present in an analysis under that convention?"""
        c = self[convention]
        self._require_verified(c)
        if not 1 <= position <= GENOME_BP:
            return False, f"outside the reference, which is 1 to {GENOME_BP:,}"
        for a, b in c.intervals:
            if a <= position <= b:
                if c.kind == "exclude":
                    return False, f"inside {a:,}-{b:,}, excluded from the analysis"
                if c.kind == "mask_then_restore":
                    return True, (f"inside {a:,}-{b:,}, which is masked during mapping "
                                  f"and restored afterwards, so the position is "
                                  f"retained")
                return False, f"inside {a:,}-{b:,}"
        return True, "not in any affected region"

    def convert(self, position: int, source: str, target: str) -> int | None:
        """The equivalent position under another convention, or None if it is not there.

        Conventions that only drop data do not move coordinates, so this returns the
        same number and `usable` is what tells you whether it means anything.
        """
        cs, ct = self[source], self[target]
        self._require_verified(cs)
        self._require_verified(ct)
        ref = cs.to_reference(position)
        if ref is None:
            return None
        if not ct.changes_coordinates:
            ok, _ = self.usable(ref, target)
            return ref if ok else None
        return ct.to_convention(ref)

    def offsets_between(self, a: str, b: str,
                        step: int = 1) -> list[tuple[int, int, int]]:
        """Where two conventions disagree, as (reference position, a, b) triples.

        Two collapsed references that remove different blocks do not differ by a
        constant, so an offset quoted as a single number is wrong somewhere.
        """
        ca, cb = self[a], self[b]
        self._require_verified(ca)
        self._require_verified(cb)
        out = []
        for pos in range(1, GENOME_BP + 1, step):
            pa, pb = ca.to_convention(pos), cb.to_convention(pos)
            if pa is not None and pb is not None and pa != pb:
                out.append((pos, pa, pb))
        return out

    def report(self, position: int, repeat_map=None) -> PositionReport:
        """Which conventions retain this position, with caveats from the repeat map.

        A convention retaining a position does not make the position meaningful. An
        exclusion list written by hand can miss part of a repeat: Lassalle's stops at
        150,962 while TRS runs to 152,222, so positions in that tail are retained by
        the convention and still sit in a duplicated region. Reporting the first
        without the second is how a reader concludes two coordinates are comparable
        when they are not.
        """
        rep = PositionReport(position)
        for c in self.conventions:
            if not c.verified:
                rep.refused.append(c.name)
                continue
            ok, why = self.usable(position, c.name)
            if ok:
                rep.usable_under.append(c.name)
            else:
                rep.unusable_under.append((c.name, why))

        if repeat_map is None:
            try:
                from .repeats import RepeatMap
                repeat_map = RepeatMap.load()
            except Exception:                      # the map is optional here
                repeat_map = None
        if repeat_map is not None and rep.usable_under:
            status = repeat_map.classify(position)
            if not status.well_defined:
                alts = ", ".join(str(a) for a in status.alternatives) or "unknown"
                rep.caveats.append(
                    f"retained by the conventions above, but the position is not well "
                    f"defined in the reference: {status.reason}. It could equally be "
                    f"{alts}")
        return rep
