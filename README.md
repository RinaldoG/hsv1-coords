# hsv1-coords

A coordinate and annotation layer for herpes simplex virus 1, so that a position in the
HSV-1 genome means one thing.

**Status: early. Nothing here is usable yet.** Milestone 0 of the plan: inputs and
skeleton only. See `PLAN.md` for what is being built and the acceptance threshold each
component must clear before it is released.

## The problem this addresses

The HSV-1 genome is well sequenced and poorly described.

- The RefSeq record `NC_001806.2` is still flagged `PROVISIONAL REFSEQ`, was last modified
  on 13 August 2018, and names `Maq v. 0.7.1` as its assembler.
- `BK012101.1` annotates the identical 152,222 bases with 300 CDS features against the
  reference's 77.
- Three deposited "strain KOS" complete genomes measure 152,011, 151,024 and 151,974 bp,
  from three different groups using three different assemblers, and none links to its raw
  reads. A coordinate cited in a KOS paper is ambiguous across roughly 987 bp.
- `UL48`, `VP16` and `alpha-TIF` are one gene, and no machine-readable dictionary says so.
- Both published HSV-1 assembly pipelines are unmaintained (VirAmp last pushed 8 April
  2015, VirGA 25 July 2018), and VirGA is still cited as the assembly method in later
  GenBank deposits.

This matters most exactly where it is least convenient. LAT, which dominates latent
transcription, lies wholly inside the inverted repeats: two indistinguishable 9,213 bp
copies, GC-rich, straddling a segment boundary that differs across the genome's four
isomers, with one annotated copy running off the start of the reference sequence (the
record itself flags this as `complement(<1..7569)`). The CTRL2 insulator sits inside the
repeat too. So the least reliable region of the sequence is the region latency happens in.

## What this provides

| Component | What it is |
|---|---|
| Gene dictionary | Every gene with its systematic name, protein name and historical aliases, each alias carrying a citation |
| Canonical annotation | One semantically versioned annotation with a changelog and per-feature provenance |
| Strain crosswalk | Coordinate lift-over between complete genomes, so a position in one strain resolves in another |
| Un-callable mask | The regions where a single linear coordinate is not well defined, enforced in code rather than noted in prose |

Static data validated in continuous integration. No server, no database, no hosted service.

## Credit

The annotation work this builds on is Whisnant et al. 2020, whose re-annotation
(Zenodo `10.5281/zenodo.3724023`) is the best existing description of the HSV-1
transcriptome and is an input here, not a thing being replaced. Provenance records on
derived files name it per feature.

## Licence

MIT, see `LICENSE`. Input data retains the terms of its source database.
