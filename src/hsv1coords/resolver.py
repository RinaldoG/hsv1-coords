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
case, separators and Greek rendering. Note that a dot is not a separator here: `UL20.6`
and `UL-20-6` are different keys, because a decimal point inside a systematic name
carries meaning. An earlier version of this paragraph said dots between a word and a
digit were folded, which described a behaviour the code does not have and should not.
"""
from __future__ import annotations

import csv
import pathlib
import re
import unicodedata
from dataclasses import dataclass, field

DEFAULT_TABLE = (pathlib.Path(__file__).resolve().parents[2]
                 / "data" / "derived" / "aliases_literature.tsv")

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

    @classmethod
    def load(cls, path: pathlib.Path | None = None) -> "Dictionary":
        path = path or DEFAULT_TABLE
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

    def is_ambiguous(self, name: str) -> bool:
        """True for a name used for two loci, or one the table declines to resolve.

        Round 10: this was False for `ICP10` and `AP`, both flagged `ambiguous` in the
        table, because it only counted how many canonicals a key maps to.
        """
        key = normalise(name)
        hits = self._by_key.get(key, set())
        if key in self._declined:
            return True
        return len(hits) > 1 and repeat_pair_of(hits) is None

    def confidence_of(self, name: str) -> str:
        """The weakest confidence recorded for this name, or '' if unknown.

        The column existed and no caller could see it.
        """
        key = normalise(name)
        order = {"established": 0, "probable": 1, "ambiguous": 2}
        found = [e.confidence for e in self.entries if normalise(e.alias) == key]
        if not found:
            return ""
        return max(found, key=lambda c: order.get(c, 3))

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
