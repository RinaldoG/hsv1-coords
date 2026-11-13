"""Resolve an HSV-1 gene name, in any convention, to a canonical systematic name.

The dictionary is `data/derived/aliases_literature.tsv`, where every alias carries a
citation. This module adds no synonymy of its own: if a mapping is not in the table with a
source, this module will not produce it.

Greek letters are a practical problem rather than a biological one. MEDLINE and PubMed
Central render them as spelled-out words, so the same name reaches a text pipeline as
`alpha TIF`, `alphaTIF`, `a-TIF` or with the actual character, depending on where it came
from. `normalise` folds those together so that a lookup works on whichever form arrives.
It is deliberately conservative: it never merges two names that differ in anything but
case, separators and Greek rendering.
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


@dataclass
class Dictionary:
    entries: list[Entry] = field(default_factory=list)
    _by_key: dict[str, set[str]] = field(default_factory=dict)

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
                    notes=(row.get("notes") or "").strip(),
                ))
        d = cls(entries=entries)
        for e in entries:
            d._by_key.setdefault(normalise(e.alias), set()).add(e.canonical)
        return d

    @property
    def canonical_names(self) -> set[str]:
        return {e.canonical for e in self.entries}

    def resolve(self, name: str) -> str:
        """The canonical name for `name`.

        Raises KeyError if unknown, and AmbiguousName if the literature uses the name
        for more than one gene. Refusing to guess is the point: an ambiguous name is a
        recorded fact about the literature, not a defect to be papered over.
        """
        hits = self._by_key.get(normalise(name))
        if not hits:
            raise KeyError(name)
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
        """True only for a name genuinely used for two different loci."""
        hits = self._by_key.get(normalise(name), set())
        return len(hits) > 1 and repeat_pair_of(hits) is None

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
