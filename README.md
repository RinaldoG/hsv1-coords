# hsv1-coords

A coordinate and annotation layer for herpes simplex virus 1, so that a position in the
HSV-1 genome means one thing.

> ## Review round 7, and what it changed, 13 November 2026
>
> An internal review found three faults that made parts of this repository wrong. They
> were published here before they were fixed, and this note is kept rather than deleted
> so the record of what was wrong survives the fix.
>
> 1. **`lift()` returned wrong coordinates for `ON783216.1`, `ON783217.1` and
>    `PQ569920.1`**, by a median of 30 to 50 kb. These three genomes are deposited on the
>    opposite strand, so their block files held offsets into the reverse complement, and
>    the only record of that sat in a column of a different file which nothing read.
>    **Fixed:** the frame is now stated on the first line of each block file and applied
>    where a deposited coordinate enters. Every one of the 183 genomes now agrees with the
>    reference at 97 per cent or better over 717,688 scored positions; these three were at
>    28 per cent. Re-introducing the fault fails a named test.
> 2. **The lift-over did not decline inside the inverted repeats, although this file said
>    it did.** The aligner runs with `--secondary=no`, so the second placement is
>    discarded before the ambiguity filter can see it, and roughly 92 per cent of
>    positions in those repeats still received one confident answer. **Fixed:** `lift`
>    consults the repeat map, which already held verified counterparts, and declines where
>    a position has more than one valid copy. `lift_detailed` returns the number together
>    with its alternatives, and `Crosswalk.to_reference` still gives the alignment's own
>    answer for anyone who wants it.
> 3. **The documentation gate was smaller than its own description, and two of the tests
>    reporting its coverage could not fail.** **Partly fixed, and measured:** the catch
>    rate over 90 single-figure mutations went from 34 to 45. Every figure this repository
>    derives is gated and the suite asserts that; the remaining misses are dates, page
>    numbers, DOIs and counts quoted from other people's records, which this approach
>    cannot pin. That limit is stated in `ROADMAP.md` rather than closed.
>
> The convention resolver, the repeat map, the version crosswalk and the gene dictionary
> were re-derived from primary sources during the same review and were unaffected.
> `FINDINGS.md` carries the full list, and the private working log carries the root
> causes.

**Status: early, and usable in part.** Seven components are built and gated: the gene
dictionary, the annotation, the strain crosswalk, the repeat map, the convention
resolver, the version crosswalk and the documentation gate that keeps the figures in
this file true. Two acceptance criteria are outstanding: the annotation is still a
classed feature table rather than GFF3, and no published coordinate has yet been
resolved from the literature by hand across several papers and strains. See `ROADMAP.md`
for the criteria
each component must clear, and `FINDINGS.md` for what building it established, including
corrections to this repository's own earlier claims.

```
pip install -e .          # or run anything below with PYTHONPATH=src
make inputs genomes all   # fetches the records and builds the artefacts
```

Every output below is copied from a real run, not written by hand.

```python
from hsv1coords import ConventionSet, Dictionary, RepeatMap, lift
from hsv1coords.liftover import lift_detailed

Dictionary.load().resolve("VP16")                        # 'UL48'
Dictionary.load().resolve("alpha-TIF")                   # 'UL48'

# A position outside the repeats lifts to a number.
lift(60000, "JQ673480.1", "JQ780693.1")                  # 59722

# One inside a repeat returns None, because there is no single honest answer.
lift(118777, "NC_001806.2", "JQ780693.1")                # None
print(lift_detailed(118777, "NC_001806.2", "JQ780693.1"))
# NC_001806.2 -> JQ780693.1: 118185, but NOT well defined (LAT,
# regulatory@complement(7592..7597) / LAT): inside IRL, which is present in two
# copies. The reference position could equally be 7597

# Three deposits are on the opposite strand, and a lift says so.

print(lift_detailed(50000, "ON783216.1", "NC_001806.2"))
# ON783216.1 -> NC_001806.2: 75929, ON THE OPPOSITE STRAND, so the base here is
# the complement of the base there

print(RepeatMap.load().classify(118777))                 # the LAT TATA box
# 118777: not well defined (inside IRL, which is present in two copies);
# could equally be 7597

print(ConventionSet.load().report(118777))
# reference position 118777:
#   usable under full_length
#   usable under whisnant2020_masked
#   usable under szpara2014_trimmed
#   NOT usable under lassalle2020_excluded: inside 116,926-132,647, excluded from
#   the analysis
#   NOT usable under morse2017_nr: inside 117,160-125,973, which is not in the
#   NR-reference (as NC_001806.1 position 118,776, offset -1)
```

So one position in the gene that dominates latent transcription becomes 109,564 in a
Szpara-trimmed analysis, is excluded from a second, absent from a third, and retained by
a fourth, while also naming no unique locus in the reference. The 109,564 is computed
here from the rule Szpara et al. state, not read off their record: that paper prints no
coordinates for the trimmed genome. The Lassalle intervals come from a figure legend
rather than Methods and look like gene boundaries rather than the annotated repeat
extents, which `data/conventions.tsv` records against that row.

The lift-over declines rather than guessing. Between the two deposited genomes labelled
strain KOS the offset is not constant: measured across the alignment it ranges from -79 to
-887 bp depending on position, so a fixed correction would be wrong almost everywhere.

## The problem this addresses

The HSV-1 genome is well sequenced and poorly described.

- The RefSeq record `NC_001806.2` is still flagged `PROVISIONAL REFSEQ`, was last modified
  on 13 August 2018, and names `Maq v. 0.7.1` as its assembler.
- `BK012101.1` annotates the identical 152,222 bp of sequence with 300 CDS features
  against the
  reference's 77.
- Four deposited genomes carry "KOS" in the strain label, under three different label
  strings (`KOS`, `isolate KOS`, `isolate KOS 1.1`), at 152,011, 151,974, 151,794 and
  151,024 bp. Two of those are a deliberate pair: Colgrove et al., *Virology* 2016;
  487:215-221 sequenced KOS and KOS1.1 and published provenance markers distinguishing
  them. The problem is not that they differ, it is that nothing maps between them. Across
  the two records both labelled plainly `KOS`, the coordinate offset measured here ranges
  from -79 to -887 bp depending on where in the genome you are, so a position cited
  against "KOS" cannot be transferred by a constant correction.
- `UL48`, `VP16` and `alpha-TIF` are one gene. UniProt records that synonymy per protein
  (`P06492`), but there is no genome-wide, cited, machine-readable mapping that covers the
  systematic, protein, functional and molecular-weight conventions together with the names
  introduced by the 2020 re-annotation.
- Both published HSV-1 assembly pipelines are unmaintained (VirAmp last pushed 8 April
  2015, VirGA 25 July 2018), and VirGA is still cited as the assembly method in later
  GenBank deposits.

This matters most exactly where it is least convenient. LAT, which dominates latent
transcription, sits in the long inverted repeat and extends past it: the reference
annotates one copy at `118777..127151`, which runs 778 bp beyond the end of `IRL`
(`117161..126373`) into `IRS`, and the other at `complement(<1..7569)`, truncated because
the genome is circular and the remainder lies in `TRS` at the other end of the linear
record. So the gene that dominates latent transcription is duplicated, spans a segment
boundary, and is split across the two ends of the reference sequence.

## What this provides

| Component | What it is |
|---|---|
| Gene dictionary | Every gene with its systematic name, protein name and historical aliases, each alias carrying a citation |
| Annotation | 1,779 feature rows, each carrying its class, its source record and that record's provenance |
| Strain crosswalk | Coordinate lift-over across 183 complete genomes, declining rather than guessing |
| Repeat map | Which positions are not well defined, and which other positions they could equally be |
| Convention resolver | Whether a position is usable under each published handling of the repeats, and what it is called there |
| Version crosswalk | `NC_001806.1` against `.2`, whose offset is not constant |
| Documentation gate | Every load-bearing figure in these documents derived from its artefact, so a stale number fails the suite |

Static data validated in continuous integration. No server, no database, no hosted service.

## Getting started

```
make inputs     # fetch and checksum the deposited records
make genomes    # fetch the pinned set of 183 complete genomes
make all        # build every derived file and run the gates
```

Inputs are fetched rather than committed, so `make test` on a fresh clone skips the two
test classes that need them and tells you which command to run.

## How this repository treats its own mistakes

`FINDINGS.md` records what building this established, including corrections to claims
this repository previously made, with the old wording kept rather than quietly replaced.

Not one of those corrections was found by the test suite passing. Several of the tests
were themselves the defect: a gate that is a mathematical identity returns 100 per cent
on anything. So the gates that carry the most weight are mutation-tested against
deliberately broken inputs: the strain crosswalk, the convention resolver and the
documentation gate, which reports how many figure-in-document mutations it catches on
every run. The M1, M2 and M4 gates are not yet mutation-tested, and a check measured as
having no power is recorded as a sanity check rather than counted as a gate.

## Credit

The annotation work this builds on is Whisnant et al. 2020, whose re-annotation
(Zenodo `10.5281/zenodo.3724023`) is the best existing description of the HSV-1
transcriptome and is an input here, not a thing being replaced. Provenance records on
derived files name it per feature.

## Licence

MIT, see `LICENSE`. Input data retains the terms of its source database.
