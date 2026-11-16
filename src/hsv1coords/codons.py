"""Codon numbers to reference coordinates, because that is how the literature cites.

Round 18: a reviewer worked a real clinical paper through this library and could not.
Brunnemann et al. (2016) Sci Rep 6:29903, PMC4947914, reports aciclovir resistance as
"G61A and P84L mediate broad cross-resistance against aciclovir, penciclovir,
brivudin", measured in single substitutions introduced into a BAC-cloned strain 17+ and
assayed by cytopathic-effect reduction in cell culture. Those are amino acid positions
in UL23's thymidine kinase, and the whole resistance literature is written that way.
This repository held `UL23 CDS complement(46673..47803)` and no way to get from codon
61 to a coordinate, so the reviewer wrote a GenBank parser, a codon table and the
minus-strand arithmetic themselves. That is the gap this closes.

What it does NOT do: translate. It gives the three reference positions of a codon, in
the order the coding strand reads them, and says which strand that is. A caller wanting
the amino acid reads those bases from the record. Nothing here asserts that a published
substitution is correct: it answers "which positions does that paper's codon number
mean in NC_001806.2".
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass

from ._coordinates import as_coordinate
from ._paths import data_path, require

ANNOTATION = data_path("derived", "annotation_v0.1.0.tsv")
REFERENCE = "NC_001806.2"
SPAN = re.compile(r"(\d+)\.\.>?(\d+)")


@dataclass(frozen=True)
class Codon:
    """Where one codon of one gene sits in the reference."""

    gene: str
    codon: int
    positions: tuple[int, int, int]
    """The three reference positions, in the order the coding strand reads them.

    For a gene on the complement these DESCEND, which is the point: a caller taking
    `min(positions)` and adding 2 would read the wrong three bases for about half the
    genes in this record.
    """
    strand: str          # "forward" or "complement"
    codons_in_gene: int

    @property
    def spans_an_intron(self) -> bool:
        """True when the three positions are not contiguous.

        39 CDS rows in this record are spliced, RL2 among them, so a codon can sit
        across a junction. A caller slicing the reference between min and max would
        pick up intron bases without noticing.
        """
        ordered = sorted(self.positions)
        return ordered[2] - ordered[0] != 2


def _cds_rows() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    with require(ANNOTATION).open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["source"] != REFERENCE or row["feature_type"] != "CDS":
                continue
            out.setdefault(row["name"], []).append(row)
    return out


def coding_positions(gene: str) -> tuple[list[int], str]:
    """Every coding position of `gene`, in reading order, and its strand.

    Raises KeyError when the reference states no CDS under that name, and ValueError
    when it states more than one, which is the two-copy case: `RL2` is annotated in
    both repeat copies, and a codon number means a position in each. The caller picks
    a copy by asking for `RL2` and reading the error, which names both locations.
    """
    rows = _cds_rows()
    if gene not in rows:
        raise KeyError(
            f"{REFERENCE} states no CDS named {gene!r}. Names come from the record's "
            f"own `/gene` qualifier; `hsv1coords.Dictionary.resolve` turns a "
            f"published name such as 'thymidine kinase' into one of them.")
    if len(rows[gene]) > 1:
        raise ValueError(
            f"{gene} has {len(rows[gene])} CDS features in {REFERENCE}, because it "
            f"lies in a repeat and is present in two copies: "
            f"{[r['location'] for r in rows[gene]]}. A codon number means a position "
            f"in each, and which one a paper meant is not recoverable from the codon "
            f"number alone, so this refuses rather than choosing.")
    location = rows[gene][0]["location"]
    strand = "complement" if location.startswith("complement(") else "forward"
    spans = [(int(a), int(b)) for a, b in SPAN.findall(location)]
    positions: list[int] = []
    for start, end in spans:
        positions.extend(range(start, end + 1))
    if strand == "complement":
        positions.reverse()
    return positions, strand


def codon(gene: str, number: int) -> Codon:
    """The reference positions of codon `number` of `gene`, 1-based.

    >>> codon("UL23", 61).positions
    (46620, 46619, 46618)

    The published labels this answers are the ones in the resistance literature:
    `codon("UL23", 61)` for G61A and `codon("UL23", 84)` for P84L.
    """
    number = as_coordinate(number, "a codon number")
    positions, strand = coding_positions(gene)
    whole = len(positions) // 3
    if not 1 <= number <= whole:
        raise ValueError(
            f"{gene}'s CDS is {len(positions):,} bases, which is {whole:,} codons "
            f"including the stop, so codon {number} is not in it.")
    start = (number - 1) * 3
    three = positions[start:start + 3]
    return Codon(gene=gene, codon=number, positions=(three[0], three[1], three[2]),
                 strand=strand, codons_in_gene=whole)
