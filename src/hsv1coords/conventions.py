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
# Each reference version has its own length. A bounds check must use the length of the
# version the coordinate is actually in. An earlier version of this module used a
# `GENOME_BP + 100` slack so that NC_001806.1 positions would pass a NC_001806.2 bounds
# check, which was the visible symptom of coordinates being applied in the wrong frame.
VERSION_LENGTHS = {"NC_001806.2": 152222, "NC_001806.1": 152261}


class UnverifiedConvention(RuntimeError):
    """Raised when asked to apply a convention this resolver cannot apply safely.

    Two cases. Either the definition was never read from its source, or it was read
    but the coordinates cannot be applied here, which is what `definition_only` means.
    `morse2017_nr` is the second case: it is defined against NC_001806.1 and everything
    else here is .2, and the offset between those two versions is not constant.
    """


class ReferenceVersionMismatch(UnverifiedConvention):
    """Raised when two conventions are defined against different reference versions.

    NC_001806.1 is 152,261 bp and NC_001806.2 is 152,222 bp, and the difference is not
    a constant shift. Reading a `.1` position in the `.2` frame, there are nine distinct
    offsets across 22 aligned blocks. Measured at four positions, one per segment:

        offset at 60,000 is +1
        offset at 120,000 is +2
        offset at 130,000 is +3
        offset at 150,000 is -39

    Converting between conventions across that boundary without a crosswalk would
    produce confident, wrong coordinates, which is the error this repository exists to
    prevent.

    Round 8: every sign in this paragraph was inverted and the magnitudes mis-located,
    against an artefact the documentation gate was already pinning, in the text a user
    reads to decide whether to trust a conversion. The gate reads the four top-level
    documents only, so nothing under `src/` was covered. Each line above is now checked
    against `VersionCrosswalk.offset_at` at that position, which is a positional check:
    a first attempt only tested that each stated offset existed somewhere in the set,
    and flipping +1 to -1 passed it, because -1 is a real offset elsewhere in the
    genome.
    """


@dataclass(frozen=True)
class Convention:
    name: str
    kind: str                  # none, exclude, mask_then_restore, trim
    intervals: list[tuple[int, int]]
    reference: str
    verified: bool
    source: str
    quote: str
    reference_version: str = "NC_001806.2"
    applicable: bool = True      # False when the definition is known but unusable here
    records: list[tuple[str, int, int]] = field(default_factory=list)
    """Named records, for a convention whose artefact is more than one sequence.

    Morse's NR-reference is four FASTA records (RL, UL, RS, US), confirmed from the
    CHROM column of its own deposited VCFs, so a Morse coordinate is a pair like
    `UL:50788` and not an integer. An earlier version of this module returned a single
    integer for it, copied from the model used for Szpara's single-sequence trim. That
    integer appears in no Morse file.
    """
    boundary_uncertainty_bp: int = 0
    """How far an interval boundary might be out.

    Non-zero when the boundaries were reconstructed rather than read. An answer for a
    position within this distance of a boundary is not reliable, and saying so is the
    difference between a usable record and a confidently wrong one.
    """

    boundaries: str = "read"
    """Whether the intervals were read from the source or reconstructed from the record.

    `read` means the source prints the coordinates, so the only question about them is
    transcription. `inferred` means the source prints none and the interval comes from
    the reference record's own annotation together with what the source says in words.

    Round 9 found `boundary_uncertainty_bp` set to 0 on all five conventions, including
    three whose own source field says the boundaries were not printed, so `near_a_boundary`
    was False everywhere and the UNRELIABLE path was unreachable on the shipped data. The
    field was the wrong instrument: Whisnant and Szpara print no coordinates, but the
    interval taken for them is the record's own TRL and TRS extent, which is exact to the
    base. What is uncertain there is the identification, not the edge, and a tolerance in
    base pairs cannot express "this might be the wrong repeat". So the uncertainty that
    applies is stated as a kind rather than a distance, and it reaches the user through
    `usable` and `report`.

    Last in the dataclass, as standing rule 22 requires.
    """

    boundaries_evidence: str = ""
    """The words in this row's own source that justify the `boundaries` flag.

    Round 10: the gate on that flag matched four phrases of this repository's prose, so
    rewording a source field honestly turned the suite red on correct data, and a
    reworded row could also have been mislabelled without the gate noticing. The
    evidence is quoted per row now, and the suite checks that the quote appears in the
    source column, which ties the flag to the source rather than to a regex over our own
    writing.

    Last in the dataclass, as standing rule 22 requires.
    """

    @property
    def boundaries_are_inferred(self) -> bool:
        return self.boundaries == "inferred"

    def near_a_boundary(self, position: int) -> bool:
        if not self.boundary_uncertainty_bp:
            return False
        return any(abs(position - edge) <= self.boundary_uncertainty_bp
                   for a, b in self.intervals for edge in (a, b))

    @property
    def changes_coordinates(self) -> bool:
        """Only removing sequence moves the numbers. Masking and excluding do not."""
        return self.kind in ("trim", "split_records")

    @property
    def is_split(self) -> bool:
        return self.kind == "split_records"

    def locate(self, reference_position: int) -> tuple[str, int] | None:
        """For a split convention: the record and 1-based offset within it, or None."""
        if not self.is_split:
            raise ValueError(f"{self.name} is not a split-record convention")
        for name, a, b in self.records:
            if a <= reference_position <= b:
                return name, reference_position - a + 1
            
        return None

    @property
    def affected_bases(self) -> int:
        return sum(b - a + 1 for a, b in self.intervals)

    @property
    def genome_length(self) -> int:
        """The length of the reference version THIS convention is defined against."""
        return VERSION_LENGTHS.get(self.reference_version, GENOME_BP)

    def retained_blocks(self, genome_bp: int | None = None) -> list[tuple[int, int]]:
        """The reference intervals this convention keeps, in order."""
        genome_bp = genome_bp if genome_bp is not None else self.genome_length
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

    def length(self, genome_bp: int | None = None) -> int:
        return sum(b - a + 1 for a, b in self.retained_blocks(genome_bp))

    def coordinate(self, reference_position: int) -> str | None:
        """A coordinate as this convention's own files would print it."""
        if self.is_split:
            hit = self.locate(reference_position)
            return None if hit is None else f"{hit[0]}:{hit[1]}"
        v = self.to_convention(reference_position)
        return None if v is None else str(v)

    def to_convention(self, reference_position: int,
                      genome_bp: int | None = None) -> int | None:
        """A coordinate IN THIS CONVENTION'S OWN REFERENCE VERSION, collapsed.

        The input must already be in `self.reference_version`. Callers holding a
        base-version coordinate translate it first; `ConventionSet` does this.
        """
        if self.is_split:
            raise ValueError(
                f"{self.name} is a four-record reference, so a position in it is a "
                f"pair like 'UL:50788', not an integer. Use `coordinate` or `locate`.")
        if not self.changes_coordinates:
            return reference_position
        offset = 0
        for a, b in self.retained_blocks(genome_bp):
            if a <= reference_position <= b:
                return offset + (reference_position - a) + 1
            offset += b - a + 1
        return None

    def from_record(self, record: str, offset: int) -> int | None:
        """For a split convention: a record and 1-based offset back to the reference.

        The inverse of `locate`, and the route a user holding `UL:50788` from the
        deposited VCFs actually needs. Round 8 found there was none, and that
        `convert`'s error message sent them to `coordinate_in`, which goes the other
        way.
        """
        if not self.is_split:
            raise ValueError(f"{self.name} is not a split-record convention")
        if not isinstance(offset, int) or isinstance(offset, bool) or offset < 1:
            raise ValueError(f"an offset within a record is a positive int, got "
                             f"{offset!r}")
        for name, a, b in self.records:
            if name == record:
                return a + offset - 1 if offset <= b - a + 1 else None
        raise KeyError(
            f"{self.name} has no record {record!r}; it has "
            f"{[n for n, _, _ in self.records]}")

    def to_reference(self, position: int, genome_bp: int | None = None) -> int | None:
        """This convention's coordinate back to the reference, or None if out of range.

        Refuses for a split convention, as `to_convention` already did. Round 8 found
        the guard on only one of the pair, so `morse2017_nr.to_reference(50788)`
        returned 50788: a four-record reference treated as one collapsed sequence,
        wrong by roughly 9 kb, and an integer that appears in no Morse file. Use
        `from_record` instead.
        """
        if self.is_split:
            raise ValueError(
                f"{self.name} is a four-record reference, so a position in it is a "
                f"pair like 'UL:50788', not an integer. Use `from_record('UL', 50788)`.")
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
    refused: list[tuple[str, str]] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    inferred_intervals: list[str] = field(default_factory=list)
    """Conventions in this report whose intervals are not printed by their source.

    Round 9: `boundary_uncertainty_bp` was 0 on all five conventions, three of which take
    their intervals from the reference record's own annotation because the source prints
    none, and no answer said so. It is listed once per report rather than as a caveat per
    position, because it is a fact about the convention and attaching it to every
    position buried the caveats that are about the position.

    Last in the dataclass, as standing rule 22 requires.
    """

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
        for name, why in self.refused:
            lines.append(f"  cannot say for {name}: {why}")
        for c in self.caveats:
            lines.append(f"  caveat: {c}")
        if self.inferred_intervals:
            lines.append(
                f"  provenance: {', '.join(self.inferred_intervals)} print no "
                f"coordinates in their own text, so the intervals used are the "
                f"reference record's annotation read together with what the source says "
                f"in words. The edges are exact to the base; the identification is an "
                f"inference. See the source column of data/conventions.tsv")
        return "\n".join(lines)


def _parse_records(text: str) -> list[tuple[str, int, int]]:
    if not text.strip():
        return []
    out = []
    for part in text.split(";"):
        name, span = part.split(":")
        a, b = span.split("-")
        out.append((name.strip(), int(a), int(b)))
    return out


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
        out = []
        for r in rows:
            flag = r["verified"].strip().lower()
            # yes            verified and in this resolver's base version
            # other_version  verified, but defined against another reference version;
            #                usable only once a version crosswalk exists
            # no             definition never read from its source
            out.append(Convention(
                name=r["convention"],
                kind=r["kind"],
                intervals=_parse_intervals(r["intervals"]),
                reference=r["reference"],
                verified=flag in ("yes", "other_version"),
                source=r["source"],
                quote=r["quote"],
                reference_version=(r.get("reference_version") or "NC_001806.2").strip(),
                applicable=flag == "yes",
                boundary_uncertainty_bp=int(r.get("boundary_uncertainty_bp") or 0),
                boundaries=(r.get("boundaries") or "read").strip(),
                boundaries_evidence=(r.get("boundaries_evidence") or "").strip(),
                records=_parse_records(r.get("records", "")),
            ))
        return cls(out)

    def __getitem__(self, name: str) -> Convention:
        return self._by_name[name]

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.conventions]

    @property
    def verified(self) -> list[Convention]:
        """Verified from source AND applicable in this resolver's coordinate space."""
        return [c for c in self.conventions if c.verified and c.applicable]

    @property
    def other_version(self) -> list[Convention]:
        """Verified, but defined against a different version of the reference."""
        return [c for c in self.conventions if c.verified and not c.applicable]

    BASE_VERSION = "NC_001806.2"

    def _require_verified(self, c: Convention) -> None:
        if not c.verified:
            raise UnverifiedConvention(
                f"{c.name} was never verified against its source, so this resolver "
                f"will not apply it. {c.source}")
        if not c.applicable and not self._version_crosswalk_available():
            raise ReferenceVersionMismatch(
                f"{c.name} is defined against {c.reference_version}, and this resolver "
                f"works in {self.BASE_VERSION}. The two differ by a non-constant "
                f"offset, so applying it would give confident wrong answers. "
                f"{c.source}")

    @staticmethod
    def _version_crosswalk_available() -> bool:
        from .versions import DEFAULT
        return DEFAULT.exists()

    def _to_convention_version(self, position: int, c: Convention) -> int | None:
        """Translate a base-version position into the version this convention uses."""
        if c.reference_version == self.BASE_VERSION:
            return position
        from .versions import VersionCrosswalk
        return VersionCrosswalk.load().convert(position, self.BASE_VERSION,
                                               c.reference_version)

    def usable(self, position: int, convention: str) -> tuple[bool, str]:
        """Is this position present in an analysis under that convention?

        Carries the interval provenance, because this is the answer a caller branches
        on. `report` uses `usable_plain` instead and states the provenance once for the
        whole report rather than on every line.
        """
        ok, why = self.usable_plain(position, convention)
        return ok, self._with_provenance(self[convention], why)

    def usable_plain(self, position: int, convention: str) -> tuple[bool, str]:
        """Is this reference position present in an analysis under that convention?

        `position` is always in this resolver's base version. Where a convention is
        defined against another version, the position is translated first, and the
        answer says so, because the translation is not a constant shift.
        """
        c = self[convention]
        self._require_verified(c)
        if c.reference_version != self.BASE_VERSION:
            native = self._to_convention_version(position, c)
            if native is None:
                return False, (f"has no coordinate in {c.reference_version}, which is "
                               f"the version this convention is defined against")
            ok, why = self._usable_native(native, c)
            answer = (f"{why} (as {c.reference_version} position {native:,}, "
                      f"offset {native - position:+d})")
            if c.near_a_boundary(native):
                answer += (f" -- UNRELIABLE: within {c.boundary_uncertainty_bp} bp of "
                           f"an interval boundary that was reconstructed, not read")
            return ok, answer
        return self._usable_native(position, c)

    def _with_provenance(self, c: Convention, answer: str) -> str:
        """Append the interval provenance to an answer that rests on it.

        Round 10: `Convention.boundaries`' docstring said the inference "reaches the user
        through `usable` and `report`". It reached them through `report` only, so
        `usable(118777, "morse2017_nr")` named an interval that is not printed in the
        source it cites and said nothing about that. A docstring claiming a path exists
        is not the path.
        """
        if not c.boundaries_are_inferred:
            return answer
        return (f"{answer} [interval inferred: {c.name} prints no coordinates in its "
                f"own text, so this interval is the reference record's annotation read "
                f"with what the source says in words]")

    def _warning_for(self, position: int, c: Convention) -> str | None:
        """Any caveat attaching to this convention's answer, usable or not.

        `report` previously dropped these whenever the answer was "usable", which
        discarded every warning on exactly the positions where a reader would act on
        the answer.
        """
        native = (self._to_convention_version(position, c)
                  if c.reference_version != self.BASE_VERSION else position)
        if native is None:
            return None
        if c.near_a_boundary(native):
            return (f"{c.name}: within {c.boundary_uncertainty_bp} bp of a "
                    f"reconstructed interval boundary, so this answer is unreliable")
        return None

    def _usable_native(self, position: int, c: Convention) -> tuple[bool, str]:
        limit = VERSION_LENGTHS.get(c.reference_version, GENOME_BP)
        if not 1 <= position <= limit:
            return False, (f"outside {c.reference_version}, which is 1 to {limit:,}")
        if c.is_split:
            hit = c.locate(position)
            if hit is None:
                removed = next((f"{a:,}-{b:,}" for a, b in c.intervals
                                if a <= position <= b), "a removed block")
                return False, f"inside {removed}, which is not in the NR-reference"
            return True, f"present as {hit[0]}:{hit[1]:,}"
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

        **This says nothing about repeat ambiguity.** A reference position inside an
        inverted repeat has two or three equally valid counterparts, and this returns
        the conversion of the position as given, not a judgement about whether that
        position identifies one locus. `report` attaches that caveat and `lift`
        declines on it; round 7 recorded the asymmetry rather than papering over it,
        because converting a coordinate between two printed conventions is a different
        question from knowing which copy of a locus it names.
        """
        cs, ct = self[source], self[target]
        self._require_verified(cs)
        self._require_verified(ct)
        if cs.is_split or ct.is_split:
            split = cs if cs.is_split else ct
            raise ValueError(
                f"{split.name} is a multi-record reference; its coordinates are pairs "
                f"like 'UL:50788'. To go from a reference position INTO it use "
                f"`coordinate_in(position, '{split.name}')`; to come OUT of it use "
                f"`self['{split.name}'].from_record(record, offset)` and convert the "
                f"reference position from there.")
        # Translate into the source's own version, out to the base version, then into
        # the target's version. Omitting these two steps is how root cause 6b survived
        # a fix that was supposed to remove it: `usable` translated and `convert` did
        # not, so the two contradicted each other.
        native_ref = cs.to_reference(position)
        if native_ref is None:
            return None
        base = self._to_base_version(native_ref, cs)
        if base is None:
            return None
        native_target = self._to_convention_version(base, ct)
        if native_target is None:
            return None
        if not ct.changes_coordinates:
            ok, _ = self._usable_native(native_target, ct)
            return native_target if ok else None
        return ct.to_convention(native_target)

    def coordinate_in(self, position: int, convention: str) -> str | None:
        """A base-version position, as the named convention's own files would print it."""
        c = self[convention]
        self._require_verified(c)
        native = self._to_convention_version(position, c)
        return None if native is None else c.coordinate(native)

    def _to_base_version(self, position: int, c: Convention) -> int | None:
        if c.reference_version == self.BASE_VERSION:
            return position
        from .versions import VersionCrosswalk
        return VersionCrosswalk.load().convert(position, c.reference_version,
                                               self.BASE_VERSION)

    # `offsets_between` was removed in round 7. It had no caller and no test, and it
    # could not do what its docstring claimed: it was written to compare two collapsed
    # references, there is only one collapsed non-split convention in the set, and the
    # other collapsed one is split, so every pair it existed for raised before
    # comparing anything. It also reported "no disagreement" for an exclude convention,
    # because it had no `usable` filter, and dropped the trimmed positions that were
    # the interesting part. Use `report` for a position, or `convert` for a pair.

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
                rep.refused.append(
                    (c.name, "its definition has not been read from its source"))
                continue
            if not c.applicable and not self._version_crosswalk_available():
                rep.refused.append(
                    (c.name, f"it is defined against {c.reference_version}, which "
                             f"differs from {self.BASE_VERSION} by a non-constant "
                             f"offset, and no version crosswalk is built"))
                continue
            ok, why = self.usable_plain(position, c.name)
            if c.boundaries_are_inferred:
                rep.inferred_intervals.append(c.name)
            warning = self._warning_for(position, c)
            if warning:
                rep.caveats.append(warning)
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
                # `describe_alternative`, not `str`: every alternative in the large
                # inverted repeats is on the opposite strand, and this report printed
                # the bare number. Round 9.
                alts = ", ".join(status.describe_alternative(a)
                                 for a in status.alternatives) or "unknown"
                rep.caveats.append(
                    f"retained by the conventions above, but the position is not well "
                    f"defined in the reference: {status.reason}. It could equally be "
                    f"{alts}")
        return rep
