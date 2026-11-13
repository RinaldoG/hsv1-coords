# hsv1-coords

A coordinate and annotation layer for herpes simplex virus 1, so that a position in the
HSV-1 genome means one thing.

**Status: early, and usable in part.** Six components are built and gated: the gene
dictionary, the annotation, the strain crosswalk, the repeat map, the convention resolver
and the version crosswalk. The annotation is still a classed feature table rather than
GFF3, which is the one acceptance criterion outstanding. See `ROADMAP.md` for the criteria
each component must clear, and `FINDINGS.md` for what building it established, including
corrections to this repository's own earlier claims.

```python
from hsv1coords import ConventionSet, Dictionary, RepeatMap, lift

Dictionary.load().resolve("VP16")              # 'UL48'
Dictionary.load().resolve("alpha-TIF")         # 'UL48'

lift(118777, "JQ673480.1", "JQ780693.1")       # 118295, or None where it cannot say

RepeatMap.load().classify(118777)
# not well defined (inside IRL, which is present in two copies); could equally be 7597

ConventionSet.load().report(118777)            # the LAT TATA box
# usable under full_length, whisnant2020_masked, szpara2014_trimmed
# NOT usable under lassalle2020_excluded: excluded from the analysis
# NOT usable under morse2017_nr: not in the NR-reference
```

So one position in the gene that dominates latent transcription is printed as 109,564 in
one published analysis, excluded from a second, absent from a third, and retained by a
fourth, while also naming no unique locus in the reference.

The lift-over declines rather than guessing. Between the two deposited genomes labelled
strain KOS the offset is not constant: measured across the alignment it ranges from -79 to
-887 bp depending on position, so a fixed correction would be wrong almost everywhere.

## The problem this addresses

The HSV-1 genome is well sequenced and poorly described.

- The RefSeq record `NC_001806.2` is still flagged `PROVISIONAL REFSEQ`, was last modified
  on 13 August 2018, and names `Maq v. 0.7.1` as its assembler.
- `BK012101.1` annotates the identical 152,222 bases with 300 CDS features against the
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
| Annotation | 1,779 features, each carrying its class, its source record and that record's provenance |
| Strain crosswalk | Coordinate lift-over across 183 complete genomes, declining rather than guessing |
| Repeat map | Which positions are not well defined, and which other positions they could equally be |
| Convention resolver | Whether a position is usable under each published handling of the repeats, and what it is called there |
| Version crosswalk | `NC_001806.1` against `.2`, whose offset is not constant |
| Un-callable mask | The regions where a single linear coordinate is not well defined, enforced in code rather than noted in prose |

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
on anything. So the gates here are mutation-tested, each against deliberately broken
inputs, and any check measured as having no power is recorded as a sanity check rather
than counted as a gate.

## Credit

The annotation work this builds on is Whisnant et al. 2020, whose re-annotation
(Zenodo `10.5281/zenodo.3724023`) is the best existing description of the HSV-1
transcriptome and is an input here, not a thing being replaced. Provenance records on
derived files name it per feature.

## Licence

MIT, see `LICENSE`. Input data retains the terms of its source database.
