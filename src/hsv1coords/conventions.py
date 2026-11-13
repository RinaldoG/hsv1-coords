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

**A convention whose definition could not be read from its own source is refused.**
`szpara2014_trimmed` is declared and not applied, because the paper is paywalled and
the trimming intervals were never verified. This project has twice encoded something
from a second-hand description that turned out to be wrong, so the rule is now explicit
in the code rather than left to judgement.
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


@dataclass
class PositionReport:
    """What one reference position means under every convention."""

    position: int
    usable_under: list[str] = field(default_factory=list)
    unusable_under: list[tuple[str, str]] = field(default_factory=list)
    refused: list[str] = field(default_factory=list)

    @property
    def comparable_everywhere(self) -> bool:
        return not self.unusable_under and not self.refused

    def __str__(self) -> str:
        lines = [f"reference position {self.position}:"]
        for name in self.usable_under:
            lines.append(f"  usable under {name}")
        for name, why in self.unusable_under:
            lines.append(f"  NOT usable under {name}: {why}")
        for name in self.refused:
            lines.append(f"  cannot say for {name}: its definition is unverified")
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
        if cs.changes_coordinates or ct.changes_coordinates:
            raise UnverifiedConvention(
                "coordinate-shifting conventions are not implemented, because the only "
                "one encountered so far could not be verified from its source")
        ok, _ = self.usable(position, target)
        return position if ok else None

    def report(self, position: int) -> PositionReport:
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
        return rep
