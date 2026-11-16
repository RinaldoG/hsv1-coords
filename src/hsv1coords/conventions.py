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
`TRL` and `TRS`. Morse removes `IRL` and `TRS` instead, and splits what remains into four
named records, so a Morse coordinate is a pair like `UL:50788` rather than an integer.

Measured, for the reference positions where both can be expressed: in the long unique
region a Morse offset is one base higher than the Szpara number for the same place, and
in the short region the two differ by 116,760 bases (`RS`) or 123,391 bases (`US`),
because Morse restarts numbering at each record and Szpara does not. So a bare number
carried from one into the other is wrong by one base over two thirds of the genome and by
more than a hundred kilobases in the rest.

An earlier version of this paragraph said both references are "published without an
offset table" and described Morse with the single-sequence model this module raises on.
Round 10 corrected both: Morse's Table 1 names a BED file of the NR-reference
coordinates, which this repository has not retrieved and which would be the better source
than the stated rule it derives from, and the module has refused to return a bare integer
for Morse since round 7. What `data/conventions.tsv` records for that row is what was
available and what was derived.
"""
from __future__ import annotations

import csv
import functools
import pathlib
import types

from ._coordinates import as_coordinate
from ._paths import data_path, require
from dataclasses import dataclass, field
from typing import Sequence

DATA = data_path("conventions.tsv")
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


_INTERVAL_BASES = types.MappingProxyType({
    # basis -> (whether the intervals were read, what a reader is owed)
    "deposited_record": (True, "the convention is the deposited record itself"),
    "paper_transcribed": (
        True, "transcribed by hand from a paper this repository does not hold, "
              "so the quotation in data/conventions.tsv is the evidence"),
    # Round 18: these two caveats are read by nothing today, because `basis_caveat`
    # has one caller and only for `paper_transcribed` (cause 34), and the first
    # repeated the absolute "prints no coordinates" that round 17 narrowed in the
    # artefact and in the documents. Narrowed here too, so the constant does not
    # outlive the correction.
    "inferred_from_record": (
        False, "the source prints no coordinates for these intervals, so they are "
               "the reference record's own annotation read with what the source says "
               "in words"),
    "derived_from_stated_rule": (
        False, "derived from the source's stated rule, because the coordinate "
               "file it names has not been retrieved"),
})
"""What each `interval_basis` value means, and whether it counts as read.

Round 15: `boundaries` was a two-valued flag whose `read` half rested on two
columns this repository writes, `quote` and `source`. A reviewer rewrote one
clause of `source` to dodge the contradiction check, appended one fabricated
sentence to `quote`, flipped the flag and dropped the inference caveat out of
`usable()` with the suite green but for two README quotations, which the same
reviewer then updated as any developer would. Cause 49's own lesson says a burden
checked against our own prose is not checked.

A free-text flag cannot carry that burden, so the basis is a controlled value, the
flag is derived from it, and what a reader is owed is stated per basis and
printed. `paper_transcribed` is the honest name for the one row whose numbers we
typed from a paper we do not hold, and it says so in the answer.
"""


@dataclass(frozen=True)
class Convention:
    name: str
    kind: str                  # none, exclude, mask_then_restore, trim
    intervals: tuple[tuple[int, int], ...]
    """The reference intervals this convention keeps.

    A tuple, not a list. Round 9 hardened `Crosswalk.blocks`, `RepeatMap.regions` and
    `RepeatMap.pairs` against a caller mutating a shared object, and round 12 found the
    two siblings it missed: `Convention` is frozen and this field was a mutable list, so
    `cs["szpara2014_trimmed"].intervals.append((50000, 70000))` changed
    `to_convention(60000)` from 50787 to None for every holder of that object.

    Round 13: the round 12 docstring, its test and its commit message all said
    `ConventionSet.load` is cached, and it is not, which made the test's second
    assertion vacuous (it passed because the file was re-read) and overstated the
    blast radius of the defect. One holder of a set is still enough: a set is passed
    between callers, and `facts.py` and `report` both walk it. `__post_init__` coerces
    this field whatever a constructing caller passed, so the guarantee belongs to the
    type rather than to the loader.
    """
    reference: str
    verified: bool
    source: str
    quote: str
    reference_version: str = "NC_001806.2"
    applicable: bool = True      # False when the definition is known but unusable here
    records: tuple[tuple[str, int, int], ...] = ()
    """Named records, for a convention whose artefact is more than one sequence.

    A tuple for the same reason `intervals` is one, and round 13's finding: this was
    the third mutable field in the dataclass round 12 hardened, and it carries
    Morse's four record extents, so `cs["morse2017_nr"].records.clear()` turned
    `usable(140000, "morse2017_nr")` from "present as US:7,392" into "inside a removed
    block, which is not in the NR-reference" for every later holder of that set.

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

    _basis: tuple = ()
    """`(whether the intervals were read, what a reader is owed)`, resolved once.

    Set by `__post_init__` from `INTERVAL_BASES`, so the answer belongs to the
    instance. Round 17: rebinding the class attribute used to change it.
    """

    interval_basis: str = "deposited_record"
    """Where this convention's intervals came from, as a controlled value.

    See `INTERVAL_BASES`. `boundaries` stays in the artefact as the two-valued
    summary, and the suite requires the two to agree, so neither can be edited alone.

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

    def __post_init__(self) -> None:
        """Freeze the two sequence fields whoever constructed the row.

        Round 13: `intervals` was typed as a tuple and the synthetic rows in the test
        suite pass a list, which the dataclass stores as it is given. So the
        immutability round 12 established held for rows built by `ConventionSet.load`
        and for nobody else. Coercing here makes it a property of the type.
        """
        object.__setattr__(self, "intervals",
                           tuple(tuple(i) for i in self.intervals))
        object.__setattr__(self, "records",
                           tuple(tuple(r) for r in self.records))
        # Round 17: `MappingProxyType` blocks item assignment and not rebinding, and
        # the class itself is not frozen, so `Convention.INTERVAL_BASES = {...}` in
        # one statement turned the inference caveat off for every later caller: the
        # outcome causes 49 and 68 are about, in the fix for them. The pair is
        # resolved here and stored on the frozen instance, so an object that already
        # answers with the caveat keeps answering with it whatever happens to the
        # class afterwards.
        basis = _INTERVAL_BASES.get(self.interval_basis)
        if basis is None:
            raise ValueError(
                f"{self.name}'s interval_basis is {self.interval_basis!r}, which is "
                f"not one of {sorted(_INTERVAL_BASES)}")
        object.__setattr__(self, "_basis", tuple(basis))

    INTERVAL_BASES = _INTERVAL_BASES
    """The controlled values, exposed for a caller that wants to read them.

    Round 18: resolving the pair onto the instance was not enough, because
    `ConventionSet.load` is not cached, so every public call builds fresh objects and
    rebinding this attribute still turned the inference caveat off for all of them.
    `__post_init__` reads the module-level `_INTERVAL_BASES`, which the class does not
    expose, so this name is documentation and not the mechanism. The test that was
    cited for the round 17 fix asserted only that item assignment on a mapping proxy
    raises, which was true before that fix as well: it had no power over the defect
    (standing rule 10).
    """

    def __post_init_basis__(self) -> None:          # documentation only
        """`__post_init__` does the work; this name is here to be searched for."""

    @property
    def boundaries_are_inferred(self) -> bool:
        return not self._basis[0]

    @property
    def basis_caveat(self) -> str:
        return self._basis[1]

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

    def _require_native(self, position: int, native: bool) -> None:
        """Refuse a base-version position on a convention defined in another version.

        Round 16: `locate`, `coordinate` and `from_record` took a parameter named
        `reference_position` and applied it straight to Morse's NC_001806.1 records,
        while `ConventionSet.coordinate_in` translated first, so the two disagreed at
        109,639 of the 152,222 positions, four of them naming a DIFFERENT record:
        `cs['morse2017_nr'].coordinate(60000)` gave 'UL:50788' against the set's
        'UL:50787', and at 9,213 'UL:1' against 'RL:9212'. 72 per cent is the same
        figure root cause 6b records for `convert` before its fix; the fix went to
        `ConventionSet` and not to these three methods next door.

        They are reachable and useful for a caller who already holds a native
        coordinate, which is why they take `native=True` rather than being made
        private.
        """
        if native or self.reference_version == ConventionSet.BASE_VERSION:
            return
        raise ReferenceVersionMismatch(
            f"{self.name} is defined against {self.reference_version} and this "
            f"position looks like a {ConventionSet.BASE_VERSION} one. Use "
            f"`ConventionSet.coordinate_in(position, {self.name!r})`, which "
            f"translates first, or pass `native=True` if the position is already a "
            f"{self.reference_version} coordinate.")

    def locate(self, reference_position: int,
               native: bool = False) -> tuple[str, int] | None:
        """For a split convention: the record and 1-based offset within it, or None.

        `reference_position` is a coordinate in THIS convention's own reference
        version; see `_require_native`.
        """
        if not self.is_split:
            raise ValueError(f"{self.name} is not a split-record convention")
        reference_position = as_coordinate(reference_position)
        self._require_native(reference_position, native)
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

    def coordinate(self, reference_position: int,
                   native: bool = False) -> str | None:
        """A coordinate as this convention's own files would print it.

        `reference_position` is a coordinate in THIS convention's own reference
        version; see `_require_native`.
        """
        if self.is_split:
            hit = self.locate(reference_position, native=native)
            return None if hit is None else f"{hit[0]}:{hit[1]}"
        self._require_native(reference_position, native)
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
        reference_position = as_coordinate(reference_position)
        limit = VERSION_LENGTHS.get(self.reference_version, GENOME_BP)
        if not 1 <= reference_position <= limit:
            # Round 16: there was no bounds check, so `full_length.to_convention(
            # 1000000)` answered 1,000,000 and `to_convention(-5)` answered -5,
            # which is the unchecked identity case round 9 closed in
            # `VersionCrosswalk.convert` and round 11 in `lift_detailed`.
            return None
        if not self.changes_coordinates:
            return reference_position
        offset = 0
        for a, b in self.retained_blocks(genome_bp):
            if a <= reference_position <= b:
                return offset + (reference_position - a) + 1
            offset += b - a + 1
        return None

    def from_record(self, record: str, offset: int,
                    native: bool = False) -> int | None:
        """For a split convention: a record and 1-based offset back to the reference.

        The inverse of `locate`, and the route a user holding `UL:50788` from the
        deposited VCFs actually needs. Round 8 found there was none, and that
        `convert`'s error message sent them to `coordinate_in`, which goes the other
        way.

        The answer is in this library's base version, NC_001806.2, like every other
        coordinate it hands out. Round 16: it returned the convention's own
        NC_001806.1 position under a docstring reading "back to the reference", so a
        user holding `UL:50788` got a number in a frame nothing told them about
        (standing rule 19). `native=True` returns the untranslated one, which is what
        `locate` round-trips against.
        """
        if not self.is_split:
            raise ValueError(f"{self.name} is not a split-record convention")
        try:
            offset = as_coordinate(offset, "an offset within a record")
        except TypeError as wrong:
            raise ValueError(str(wrong)) from wrong
        if offset < 1:
            raise ValueError(f"an offset within a record is a positive int, got "
                             f"{offset!r}")
        for name, a, b in self.records:
            if name == record:
                if offset > b - a + 1:
                    # Round 17: this returned the same None as a position that HAS no
                    # NC_001806.2 equivalent, and 7 offsets inside the four records
                    # are in that second case (RL:8776, RL:8777, UL:894, UL:24223,
                    # UL:53026, UL:77111 and one more), so a user holding `UL:894`
                    # from a deposited VCF could not tell "your coordinate does not
                    # exist" from "that base is not in the newer reference". The same
                    # fault `unmapped_stretch` was fixed for in round 16.
                    raise ValueError(
                        f"{record} holds {b - a + 1:,} bases, so {record}:{offset} is "
                        f"not a position in it. A None from this method means the "
                        f"position exists and has no {ConventionSet.BASE_VERSION} "
                        f"equivalent.")
                own = a + offset - 1
                if native:
                    return own
                if self.reference_version != ConventionSet.BASE_VERSION:
                    # Round 16: this returned an NC_001806.1 position under a
                    # docstring reading "back to the reference", and `convert`'s own
                    # error message sends users here. The frame travels with the
                    # value: it is translated to the base version, which is what
                    # every other coordinate this library hands out is in.
                    from .versions import VersionCrosswalk
                    return VersionCrosswalk.load().convert(
                        own, self.reference_version, ConventionSet.BASE_VERSION)
                return own
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
        limit = self.length(genome_bp)
        position = as_coordinate(position)
        if not 1 <= position <= limit:
            # Round 16: there was no bounds check, so `to_reference(0)` gave 9,213,
            # inside the TRL that `szpara2014_trimmed` removes and which
            # `to_convention(9213)` answers None for, and
            # `to_reference(-20)` gave 9,193. That is the identity-case hole round 9
            # closed in `VersionCrosswalk.convert` and round 11 in `lift_detailed`,
            # still open in the exported `Convention`.
            return None
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
    transcribed_intervals: list[tuple[str, str]] = field(default_factory=list)
    """Conventions here whose intervals ARE the source's, typed in from a paper.

    Round 15: `read` was taken to need no disclosure, and its evidence was a quotation
    in a column this repository writes. One row is in that position,
    `lassalle2020_excluded`, and the answer says so now.
    """

    checks_not_run: list[str] = field(default_factory=list)
    """Checks that could not be made, so silence is not read as a clean answer.

    Round 18: with no repeat map at all, `lift` raised and named the file while this
    report answered with `caveats == []` and no notice, at the repository's own
    worked example: a caller could not tell "checked and well defined" from "not
    checked". Root cause 20's outcome through the one door round 17 left open.

    Last in the dataclass, as standing rule 22 requires.
    """

    interval_provenance: list[tuple[str, str]] = field(default_factory=list)
    """Each convention here whose intervals its source does not print, and why.

    Round 9: `boundary_uncertainty_bp` was 0 on all five conventions, three of which take
    their intervals from the reference record's own annotation because the source prints
    none, and no answer said so. It is listed once per report rather than as a caveat per
    position, because it is a fact about the convention and attaching it to every
    position buried the caveats that are about the position.

    Round 13: this held the names only, and the sentence `__str__` printed beside them
    was hand-written prose that counted the data: "Two of those sources print no
    coordinates at all ... for morse2017_nr the paper names a coordinate file that has
    not been retrieved". Both halves are statements about rows of
    `data/conventions.tsv`, pinned in library code, so a fourth inferred row or a
    retrieved file would have made the library say something false, and
    `boundaries_evidence` sat in the dataclass read by nothing while the suite gated it
    (root cause 34). Each row's own evidence travels with its name now and is what the
    report prints.

    Last in the dataclass, as standing rule 22 requires.
    """

    @property
    def inferred_intervals(self) -> list[str]:
        """The names only, for a caller that branches rather than prints."""
        return [name for name, _ in self.interval_provenance]

    @property
    def comparable_everywhere(self) -> bool:
        """Usable under every convention this resolver can evaluate.

        Deliberately ignores the refused ones: counting them would make this False at
        every position in the genome whenever any convention is unverified, which tells
        a reader nothing. `refused` is reported separately so it is still visible.
        """
        if self.checks_not_run:
            # Round 18: this returned the strongest claim the type can make while the
            # ambiguity check had not run at all.
            return False
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
        for missing in self.checks_not_run:
            lines.append(f"  not checked: {missing}")
        if self.interval_provenance:
            each = "; ".join(f"{name} ({why})" for name, why in
                             self.interval_provenance)
            lines.append(
                f"  provenance: the intervals used here were not read from the source "
                f"for {each}. The edges are exact to the base; which features they are "
                f"is the inference. See the source column of data/conventions.tsv")
        if self.transcribed_intervals:
            each = "; ".join(f"{name} ({why})" for name, why in
                             self.transcribed_intervals)
            lines.append(
                f"  transcription: the intervals used here are the source's own "
                f"numbers, typed in by hand, for {each}")
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


def _parse_intervals(text: str) -> tuple[tuple[int, int], ...]:
    if not text.strip() or text.strip() == "UNVERIFIED":
        return ()
    out = []
    for part in text.split(";"):
        a, b = part.split("-")
        out.append((int(a), int(b)))
    return tuple(out)


class ConventionSet:
    def __init__(self, conventions: "Sequence[Convention]") -> None:
        # A tuple, round 13: this was the fourth mutable container holding the same
        # data, and `conventions` is what `facts.py`, `report` and every gate in the
        # suite walk, so one caller sorting or clearing it changes what the next one
        # sees. Callers that add a row build a new set, which is what they already do.
        self.conventions = tuple(conventions)
        self._by_name = {c.name: c for c in self.conventions}

    @classmethod
    def load(cls, path: pathlib.Path | None = None) -> "ConventionSet":
        with require(pathlib.Path(path) if path else DATA).open() as fh:
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
                interval_basis=(r.get("interval_basis")
                                or "deposited_record").strip(),
                boundaries_evidence=(r.get("boundaries_evidence") or "").strip(),
                records=_parse_records(r.get("records", "")),
            ))
        if not out:
            # Round 16: an empty table gave `comparable_everywhere` True, because
            # `not self.unusable_under` is vacuously true over zero conventions, so
            # the strongest claim the type can make came out of the weakest possible
            # data (cause 55 inside library code).
            raise ValueError(
                f"{require(pathlib.Path(path) if path else DATA)} states no "
                f"conventions, so every question about a position would be answered "
                f"by silence. Rebuild or restore the file.")
        return cls(out)

    def __getitem__(self, name: str) -> Convention:
        try:
            return self._by_name[name]
        except KeyError:
            # Round 11: a bare `KeyError: 'hg38'` with no list of what is valid, while
            # `Convention.from_record` names its four records on a bad one.
            raise KeyError(
                f"unknown convention {name!r}; this resolver holds "
                f"{sorted(self._by_name)}") from None

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
        # Loaded once. Round 12: this called `VersionCrosswalk.load()` on every
        # invocation, so every Morse query re-read and re-parsed the TSV, once per
        # position. `Crosswalk` and `RepeatMap` are cached for exactly this reason.
        return self._version_crosswalk().convert(position, self.BASE_VERSION,
                                                 c.reference_version)

    @staticmethod
    @functools.lru_cache(maxsize=1)
    def _version_crosswalk():
        from .versions import VersionCrosswalk
        return VersionCrosswalk.load()

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
        return (f"{answer} [interval not read from the source: see the source column of "
                f"data/conventions.tsv for what was available and what was derived]")

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
        # Round 11: `usable(5.5, ...)` answered while `report(5.5)` raised, which is the
        # same package disagreeing with itself about what a coordinate is. Every other
        # entry point in the library refuses a non-int.
        position = as_coordinate(position)
        limit = VERSION_LENGTHS.get(c.reference_version, GENOME_BP)
        if not 1 <= position <= limit:
            return False, (f"outside {c.reference_version}, which is 1 to {limit:,}")
        if c.is_split:
            hit = c.locate(position, native=True)
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

        "Not there" covers both ends: a position either convention drops has no
        equivalent under the pair, and `usable` is what says which end dropped it and
        why. Conventions that only exclude data do not MOVE the coordinate, so where
        both ends retain the position the number is unchanged; `coordinate_in` is the
        way to ask how one convention would print a position without asking whether a
        second one keeps it.

        Round 16 corrected this docstring and the asymmetry under it: the text said
        the number came back with `usable` carrying the caveat, the code returned None
        for an excluded target, and nothing checked the source end at all.

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
        # Round 16: the exclusion check was applied to the TARGET only, so
        # `convert(118777, 'full_length', 'lassalle2020_excluded')` was None while
        # `convert(118777, 'lassalle2020_excluded', 'full_length')` answered 118,777:
        # the same pair of conventions, opposite answers, depending on which way the
        # question was asked (cause 56). Both ends are checked now. The contract that
        # several rounds built and tested is kept: a position a convention excludes
        # has no equivalent under it, and `usable` says why.
        source_ok, _ = self._usable_native(native_ref, cs)
        if not source_ok:
            return None
        if not ct.changes_coordinates:
            ok, _ = self._usable_native(native_target, ct)
            return native_target if ok else None
        return ct.to_convention(native_target)

    def coordinate_in(self, position: int, convention: str) -> str | None:
        """A base-version position, as the named convention's own files would print it."""
        position = as_coordinate(position)
        c = self[convention]
        self._require_verified(c)
        native = self._to_convention_version(position, c)
        # `native=True`: this has already translated into the convention's own
        # version, which is the whole point of going through the set.
        return None if native is None else c.coordinate(native, native=True)

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
                rep.interval_provenance.append((c.name, c.boundaries_evidence))
            elif c.interval_basis == "paper_transcribed":
                rep.transcribed_intervals.append((c.name, c.basis_caveat))
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
            except FileNotFoundError as absent:
                # Round 16: this caught EVERY exception, so a `repeat_map.tsv` whose
                # orientation column had been dropped made `lift`, `lift_detailed`
                # and `copies_in_the_reference` all raise and name the file, while
                # `report` answered with no caveats and no notice. Round 14 narrowed
                # the identical clause in `copies_in_the_reference` for this reason
                # and the sibling here was missed (standing rules 54 and 61). Only
                # the map's ABSENCE is a reason to answer without it.
                repeat_map = None
                rep.checks_not_run.append(
                    f"the repeat map is not built, so nothing here says whether this "
                    f"position is well defined: {absent}")
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
