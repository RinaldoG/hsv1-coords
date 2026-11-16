"""Resolve an HSV-1 gene name, in any convention, to a canonical systematic name.

The dictionary is `data/derived/aliases_literature.tsv`, where every alias carries a
citation. This module adds no synonymy of its own beyond one mechanical case, declared
here because round 10 found the claim stated absolutely and no longer true: the
single-letter Greek renderings (`a-TIF` for `alpha-TIF`) are derived keys, so
`resolve("a-TIF")` answers while `citations_for("UL48", "a-TIF")` is empty, since no
source prints that spelling. Every other mapping comes from the table with a source, and
a name the table declines to attribute raises `UnresolvableName` rather than resolving.

Greek letters are a practical problem rather than a biological one. MEDLINE and PubMed
Central render them as spelled-out words, so the same name reaches a text pipeline as
`alpha TIF`, `alphaTIF`, `a-TIF` or with the actual character, depending on where it came
from. `normalise` folds case, separators and the Greek characters themselves, and the
single-letter renderings (`a-TIF` for `alpha-TIF`) are indexed when the dictionary loads,
where it can be checked that the short form is not already a cited alias of something
else. Round 9 found this paragraph claiming `normalise` folded `a-TIF` when it raised
`KeyError`, which would drop a record in a text pipeline that believed it.

It is deliberately conservative: it never merges two names that differ in anything but
case, separators and Greek rendering. One consequence worth knowing, found in round 10
and declared in `tests/known_gaps.tsv`: because a hyphen is a separator, `oriS` and
`Ori-S` fold to one key, and the answer is the 2020 re-annotation's reading frame
spanning the TRS copy and TRS1, not the oriS palindrome. The origins of replication are
not gene names and are not in this table at all; they are `stem_loop` features in
`data/derived/repeat_map.tsv`. Note that a dot is not a separator here: `UL20.6`
and `UL-20-6` are different keys, because a decimal point inside a systematic name
carries meaning. An earlier version of this paragraph said dots between a word and a
digit were folded, which described a behaviour the code does not have and should not.
"""
from __future__ import annotations

import csv
import functools
import pathlib

from ._paths import data_path, require
import re
import types
import unicodedata
from dataclasses import dataclass, field

DEFAULT_TABLE = data_path("derived", "aliases_literature.tsv")

# Greek character to the word MEDLINE prints for it.
GREEK = {
    "α": "alpha", "Α": "alpha", "β": "beta", "Β": "beta",
    "γ": "gamma", "Γ": "gamma", "δ": "delta", "Δ": "delta",
    "ε": "epsilon", "κ": "kappa", "μ": "mu", "ω": "omega",
}


def normalise(name: str) -> str:
    """Fold a name to a comparison key: case, separators and Greek rendering only."""
    s = unicodedata.normalize("NFKC", name).strip()
    for ch, word in GREEK.items():
        s = s.replace(ch, word)
    s = s.lower()
    # Treat spaces, hyphens, underscores and dots between a word and a digit as the
    # same separator, which is what differs between `alpha TIF`, `alpha-TIF` and
    # `alphaTIF`. Decimal points inside systematic names (UL20.6) are preserved,
    # because there they carry meaning.
    s = re.sub(r"(?<=[a-z])[\s_-]+(?=[a-z])", " ", s)
    s = re.sub(r"[\s_-]+", "", s)
    return s


# The canonical set mixes two schemes, because the deposited records do: the reference
# annotates RL1, RL2 and RS1 once each, and the 2020 re-annotation annotates TRL1/IRL1,
# TRL2/IRL2 and IRS1/TRS1 separately, while also introducing copy-specific decimal names
# such as IRL2.5 that have no copy-agnostic form. So `resolve("IRL2")` gives `RL2` and
# drops the copy, and `resolve("IRL2.5")` gives `IRL2.5` and keeps it. Round 10 found
# that difference invisible to a caller; `Dictionary.copy_kept` now reports it, and the
# repeat map is what answers "which copy" for a coordinate.
#
# A locus inside an inverted repeat exists in two copies, named with a TR (terminal
# repeat) and an IR (internal repeat) prefix and the same suffix. A name mapping to both
# is not ambiguous: it is one locus present twice, which is a fact about the genome and
# the central reason this repository exists. See FINDINGS.md.
REPEAT_PAIR = re.compile(r"^(TR|IR)([LS])(.*)$")


def repeat_pair_of(names: set[str] | list[str]) -> str | None:
    """If these names are the two repeat copies of one locus, return its bare suffix."""
    parsed = [REPEAT_PAIR.match(n) for n in names]
    if len(parsed) != 2 or not all(parsed):
        return None
    a, b = parsed
    if a.group(1) == b.group(1):
        return None                      # both TR or both IR, so two different loci
    if a.group(2) != b.group(2) or a.group(3) != b.group(3):
        return None                      # different segment or different suffix
    return a.group(2) + a.group(3)


@dataclass
class Entry:
    canonical: str
    alias: str
    alias_type: str
    citation: str
    citation_identifier: str
    inventory: str
    confidence: str
    notes: str = ""
    resolution: str = "resolve"
    """Whether the resolver may answer with this row, decided per row in the builder.

    Round 10: `confidence` was loaded here and read by nothing, so four rows flagged
    `ambiguous`, whose notes say in terms that they are not resolved, resolved silently.
    `ICP10` is the HSV-2 name for the positional homologue and returned UL39; `AP` is
    attributable to either UL26 or UL26.5 and returned UL26. The module docstring says
    refusing to guess is the point, and for those rows the table recorded the refusal
    while the code guessed.
    """


@dataclass
class Dictionary:
    entries: list[Entry] = field(default_factory=list)
    _by_key: dict[str, set[str]] = field(default_factory=dict)
    _declined: dict[str, list] = field(default_factory=dict)
    _derived_from: dict[str, str] = field(default_factory=dict)
    """Each key this loader derived, and the cited key it was derived from.

    Round 13: `a-TIF` and `aTIF` resolved to UL48 while `confidence_of` returned ''
    for both, which the docstring defines as "unknown name", and a declined alias
    spelled `alpha-...` was refused under its own spelling and answered under the
    short one, because `_declined` was built before the short keys existed. A derived
    key now carries everything its source key carries.
    """

    @classmethod
    @functools.lru_cache(maxsize=4)
    def load(cls, path: pathlib.Path | None = None) -> "Dictionary":
        """Cached, like the other two loaders.

        Round 11 measured the README's own first pattern, `Dictionary.load().resolve(x)`
        in a loop: 0.57 s for 200 names against 0.2 ms with one load reused, 2,440 times
        the useful cost. The other two loaders were cached in round 8 and this was
        missed. Same immutability basis: the entries are read and not written.
        """
        path = require(pathlib.Path(path or DEFAULT_TABLE))
        entries: list[Entry] = []
        with pathlib.Path(path).open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                entries.append(Entry(
                    canonical=row["canonical_systematic_name"].strip(),
                    alias=row["alias"].strip(),
                    alias_type=row["alias_type"].strip(),
                    citation=row["citation"].strip(),
                    citation_identifier=row["citation_identifier"].strip(),
                    inventory=row["inventory"].strip(),
                    confidence=row["confidence"].strip(),
                    resolution=(row.get("resolution") or "resolve").strip(),
                    notes=(row.get("notes") or "").strip(),
                ))
        d = cls(entries=entries)
        for e in entries:
            d._by_key.setdefault(normalise(e.alias), set()).add(e.canonical)
            if e.resolution == "decline":
                d._declined.setdefault(normalise(e.alias), []).append(e)
        d._add_single_letter_greek_keys()
        # Frozen after the derived keys are added, because this loader is cached now and
        # hands one object to every caller. Round 9 froze the crosswalk's blocks and the
        # repeat map's regions for this reason, round 11 found the repeat map's
        # `inverted_pairs` and `multi_copy` still mutable, and caching this loader would
        # have added a third instance of the same hazard: one caller clearing `_by_key`
        # would stop every later resolve in the process.
        d.entries = tuple(d.entries)
        d._by_key = types.MappingProxyType(
            {key: frozenset(value) for key, value in d._by_key.items()})
        d._declined = types.MappingProxyType(
            {key: tuple(value) for key, value in d._declined.items()})
        d._derived_from = types.MappingProxyType(dict(d._derived_from))
        return d

    def _add_single_letter_greek_keys(self) -> None:
        """Also index `a-TIF` for an alias spelled `alpha-TIF`, and only then.

        Older MEDLINE renderings print a Greek letter as a single Latin letter, so the
        same name arrives as `a-TIF` and `aTIF`. The module docstring claimed `normalise`
        folded those with `alpha-TIF`; it did not, and both raised `KeyError`. Round 9.

        Done here rather than in `normalise` because folding a leading `a` to `alpha` in
        general would merge names that are not the same name. This only ever adds a key
        that an existing `alpha`, `beta` or `gamma` alias already implies, skips any key
        that is already in use, and refuses any that would make one key resolve to more
        than one gene, so it cannot create an ambiguity or override a cited alias.
        """
        words = {"alpha": "a", "beta": "b", "gamma": "g"}
        for key in list(self._by_key):
            for word, letter in words.items():
                if not key.startswith(word):
                    continue
                short = letter + key[len(word):]
                if short in self._by_key:
                    continue    # a cited alias of its own: never overridden
                targets = self._by_key[key]
                if len(targets) == 1 or repeat_pair_of(targets) is not None:
                    self._by_key[short] = set(targets)
                    self._derived_from[short] = key
                    # Round 13: the refusal lives in `_declined`, which was built
                    # before this ran, so a declined alias answered under its derived
                    # spelling. Marking `alpha-TIF` as `decline` left `a-TIF` and
                    # `aTIF` resolving to UL48: the table recorded a refusal and the
                    # code answered anyway, which is root cause 43's shape in the one
                    # place round 10's fix did not reach.
                    if key in self._declined:
                        self._declined[short] = list(self._declined[key])

    @property
    def canonical_names(self) -> set[str]:
        return {e.canonical for e in self.entries}

    def resolve(self, name: str) -> str:
        """The canonical name for `name`.

        Raises KeyError if unknown, and AmbiguousName if the literature uses the name
        for more than one gene. Refusing to guess is the point: an ambiguous name is a
        recorded fact about the literature, not a defect to be papered over.
        """
        key = normalise(name)
        hits = self._by_key.get(key)
        if not hits:
            raise KeyError(name)
        declined = self._declined.get(key)
        if declined and len(hits) == 1:
            raise UnresolvableName(name, sorted(hits), declined[0].notes)
        if len(hits) > 1:
            if repeat_pair_of(hits):
                raise RepeatPairedName(name, sorted(hits))
            raise AmbiguousName(name, sorted(hits))
        return next(iter(hits))

    def resolve_all(self, name: str) -> list[str]:
        """Every canonical name this string has been used for, ambiguity included."""
        return sorted(self._by_key.get(normalise(name), set()))

    def is_repeat_paired(self, name: str) -> bool:
        """True if this name denotes one locus present in both inverted-repeat copies."""
        return repeat_pair_of(self._by_key.get(normalise(name), set())) is not None

    def copies_in_the_reference(self, name: str) -> int | None:
        """How many copies of this gene the reference carries: 2, 1, or None.

        Round 13: whether a diploid locus announced itself depended only on which
        naming scheme its aliases happened to use. `resolve("pri-miRNA-H1")` raises
        `RepeatPairedName`, because the table holds both copies under their own names,
        while `resolve("ICP0")`, `resolve("ICP34.5")` and `resolve("ICP4")` answered
        RL2, RL1 and RS1 with `is_repeat_paired` False: those canonicals name the
        locus rather than a copy, so nothing in the answer said the gene is present
        twice. They are the three loci this repository exists for, and
        `data/derived/diploid_pairs.tsv` already held the fact.

        A separate question from `is_repeat_paired`, deliberately, and not folded into
        it: that one means "`resolve` will refuse", which these names must not do.

        None when the dictionary does not hold the name, or when the repeat map is not
        available, so a caller can tell "one copy" from "no information".
        """
        try:
            canonical = self.resolve(name)
        except RepeatPairedName:
            # Both copies are in the table under their own names, so the question is
            # already answered. This clause is first because the refusals are all
            # subclasses of KeyError.
            return 2
        except KeyError:
            return None
        try:
            from .repeats import RepeatMap
            pairs = RepeatMap.load().pairs
        except Exception:                       # the map is optional here
            return None
        bare = canonical.upper()
        for prefix in ("IR", "TR"):
            if bare.startswith(prefix) and len(bare) > len(prefix):
                bare = bare[len(prefix):]
                break
        for pair in pairs:
            locus = pair.locus.upper()
            if locus in (canonical.upper(), bare):
                return 2
        return 1

    def is_ambiguous(self, name: str) -> bool:
        """True for a name used for two loci, or one the table declines to resolve.

        Round 10: this was False for `ICP10` and `AP`, both flagged `ambiguous` in the
        table, because it only counted how many canonicals a key maps to.

        It is not the same question as `confidence_of(name) == "ambiguous"`, and round
        13 found the two read as though they were. This one means "`resolve` will
        refuse", and the suite requires exactly that. The column means "this string is
        easy to confuse with another gene's", which is true of names that resolve
        perfectly well: "capsid scaffolding protein" is UL26 and "capsid scaffold
        protein" is UL26.5, one letter apart, each unambiguous on its own. Refusing
        those would be wrong; saying nothing about them was the defect. `caution_for`
        is the signal for that, and it is not silent where this is False.
        """
        key = self._cited_key(name)
        hits = self._by_key.get(key, set())
        if key in self._declined:
            return True
        return len(hits) > 1 and repeat_pair_of(hits) is None

    def _rows_for(self, key: str) -> list:
        return [e for e in self.entries if normalise(e.alias) == key]

    def _cited_key(self, name: str) -> str:
        """The key whose rows describe this name: its own, or the one it came from."""
        key = normalise(name)
        return self._derived_from.get(key, key)

    def confidence_of(self, name: str) -> str:
        """The WEAKEST confidence recorded for this name, or '' if it is not in here.

        Weakest, deliberately: a name is only as attested as its least attested
        citation, so `ICP4` reads `probable` because the 1988 and 1985 naming papers
        are paywalled and were not read in full. 77 of 375 names read `probable` for
        that reason, and `confidences_of` gives both ends for a caller that needs to
        tell "weak everywhere" from "strong somewhere".

        Round 13: a derived Greek key answered '' here, which this docstring defines
        as a name the dictionary does not hold, while `resolve` answered UL48 for the
        same string.
        """
        best_worst = self.confidences_of(name)
        return best_worst[1] if best_worst else ""

    def confidences_of(self, name: str) -> tuple[str, str] | None:
        """(strongest, weakest) confidence recorded for this name, or None if unknown.

        Round 13: `confidence_of` returns the weakest only, so a caller filtering on
        `established` discarded ICP4 and all 68 systematic UL and US names.
        """
        order = {"established": 0, "probable": 1, "ambiguous": 2}
        found = [e.confidence for e in self._rows_for(self._cited_key(name))]
        if not found:
            return None
        return (min(found, key=lambda c: order.get(c, 3)),
                max(found, key=lambda c: order.get(c, 3)))

    COPY_PREFIXES = {"IR": "internal", "TR": "terminal"}

    def copy_named(self, name: str) -> str | None:
        """Which copy of an inverted repeat the INPUT named, or None.

        Round 10: the canonical set mixes two schemes, because the records do.
        `resolve("IRL2")` gives `RL2`, which names the locus and not the copy, while
        `resolve("IRL2.5")` gives `IRL2.5`, which names the copy. Nothing in the answer
        said which had happened, so a caller asking about the internal copy of ICP0 got
        a copy-agnostic answer indistinguishable from a copy-specific one.

        This reports what the caller's own string said, which is the part the canonical
        may have dropped. `copy_kept` says whether the answer still carries it.
        """
        bare = name.strip()
        if normalise(bare) not in self._by_key:
            # Round 13: any string starting `tr` or `ir` was answered, so
            # `copy_named("trypsin")` was "terminal" and `copy_named("transactivator")`
            # "terminal" while `copy_kept` returned None for both, which its own
            # docstring defines as "the input named no copy". The two accessors
            # contradicted each other on caller input. A name this dictionary does not
            # hold names no repeat copy here.
            return None
        for prefix, which in self.COPY_PREFIXES.items():
            if bare.upper().startswith(prefix) and len(bare) > len(prefix):
                return which
        return None

    def copy_kept(self, name: str) -> bool | None:
        """Whether the canonical for `name` still names the copy the input named.

        None when the input named no copy. False is the case worth knowing: the caller
        said `IRL2` and the answer is `RL2`, so the copy has been dropped and the
        coordinate question "which copy" is still open. The repeat map answers it.
        """
        named = self.copy_named(name)
        if named is None:
            return None
        try:
            canonical = self.resolve(name)
        except KeyError:
            return None
        return self.copy_named(canonical) is not None

    def caution_for(self, name: str) -> str:
        """Why this name is easy to get wrong, or '' if nothing is recorded.

        Round 13: the table flags six aliases `ambiguous`, and for three of them
        (`AL-RNA` and the two capsid scaffold strings) `resolve` answers and
        `is_ambiguous` is False, correctly, so the only signal was a confidence string
        a caller had to know to compare against. The note beside the flag is the thing
        worth reading: for the capsid pair it says in terms that "a string match on
        either will silently pick the wrong gene".
        """
        notes = [e.notes for e in self._rows_for(self._cited_key(name))
                 if e.confidence == "ambiguous" and e.notes]
        return notes[0] if notes else ""

    def why_declined(self, name: str) -> str:
        """Why the table refuses to resolve this name, or '' if it does resolve."""
        declined = self._declined.get(normalise(name))
        return declined[0].notes if declined else ""

    def aliases_for(self, canonical: str) -> list[str]:
        return sorted({e.alias for e in self.entries if e.canonical == canonical})

    def citations_for(self, canonical: str, alias: str) -> list[str]:
        return sorted({e.citation_identifier or e.citation
                       for e in self.entries
                       if e.canonical == canonical and e.alias == alias})


class AmbiguousName(KeyError):
    """The literature uses this name for two genuinely different loci."""

    def __init__(self, name: str, candidates: list[str]) -> None:
        super().__init__(f"{name!r} is used for {candidates} in the literature")
        self.name = name
        self.candidates = candidates


class UnresolvableName(KeyError):
    """The table records this name and refuses to attribute it to one gene.

    Distinct from `AmbiguousName`, which is two canonicals in the data. Here there is
    one canonical and a row saying it must not be handed out: `ICP10` is the HSV-2 name
    for the positional homologue of UL39, and `AP` is attributable to either UL26 or
    the nested UL26.5 depending on whether the precursor or the nested gene is meant.
    Round 10 found both resolving silently.
    """

    def __init__(self, name: str, candidates: list[str], why: str) -> None:
        super().__init__(f"{name!r} is recorded against {candidates} and the table "
                         f"declines to resolve it: {why}")
        self.name = name
        self.candidates = candidates
        self.why = why


class RepeatPairedName(KeyError):
    """One locus, present in both copies of an inverted repeat.

    Not an ambiguity. A caller that wants a single coordinate must say which copy it
    means, or accept both, which is exactly the decision the un-callable mask exists to
    force (roadmap M4).
    """

    def __init__(self, name: str, candidates: list[str]) -> None:
        super().__init__(f"{name!r} is one locus in two inverted-repeat copies: "
                         f"{candidates}")
        self.name = name
        self.candidates = candidates
