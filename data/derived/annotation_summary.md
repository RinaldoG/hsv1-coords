# Annotation v0.1.0: features by class and source

Produced by `src/build_annotation.py`. Every feature carries the class of thing it is and the record it came from.

**Why this table exists.** The two deposited annotations of this genome are often described as disagreeing, with 77 coding sequences against 300 on the identical 152,222 bases. They do not disagree. They describe different classes of feature, and a count is only comparable within a class. See `FINDINGS.md`.

**How to read it.** The first table counts **distinct names**, which is the only figure comparable across these sources. A row count is not: `NC_001806.2` states a gene and a CDS and sometimes an mRNA for one locus, `BK012101.1` states a gene and an mRNA, and the GTF states a transcript and each of its exons, so the same locus contributes a different number of rows to each. The second table gives the rows by feature type, where a count is comparable because the type is fixed. Until round 10 this file blended the two, and reported 30 against 479 in one class for two formats of the same annotation by the same authors.

## Distinct names by class

| Class | NC_001806.2 | BK012101.1 | whisnant2020-gtf | What it is |
|---|---|---|---|---|
| `dORF` | 0 | 1 | 0 | A downstream open reading frame, 3-prime of a main coding sequence |
| `iORF` | 0 | 1 | 0 | An internal open reading frame, inside another coding sequence |
| `lat_region` | 1 | 4 | 5 | A named transcript or open reading frame in the LAT region, which is named by function rather than by position |
| `misc_feature` | 2 | 0 | 0 | A feature the source record did not type further |
| `named_gene` | 74 | 6 | 77 | A classically named gene in the systematic positional scheme. The trailing letter is part of the scheme: UL49A (gN) and US8A are genes, and were classed `other` |
| `novel_orf_2020` | 0 | 18 | 85 | An open reading frame resolved by the 2020 re-annotation and numbered by position, with no independent name in the literature. A decimal name that the 2018 reference record already carries is a named gene instead, which is what UL26.5 (ICP35, VP22a, described in 1974) is |
| `origin_locus` | 0 | 1 | 1 | A locus named after an origin of replication. Note that `Ori-S` here is the 2020 re-annotation's locus spanning the TRS copy and TRS1, not the oriS palindrome itself, which is a `stem_loop` feature |
| `primary_mirna` | 0 | 0 | 2 | A primary microRNA transcript |
| `regulatory` | 57 | 0 | 0 | A regulatory element. In this reference record 53 are polyadenylation signals and 19 are TATA boxes, which the `subtype` column now distinguishes; the earlier description said "such as a TATA box" of a set that is three quarters polyA signals |
| `repeat_region` | 21 | 0 | 0 | A repeat, including the inverted repeats and the reiteration sets nested inside them |
| `sORF` | 0 | 11 | 0 | A small open reading frame named by the 2020 re-annotation |
| `stem_loop` | 2 | 0 | 0 | A stem-loop structure, such as an origin palindrome |
| `transcript_span` | 0 | 37 | 0 | A transcript unit spanning more than one named gene, so not a gene itself. Both halves must be systematic names: the earlier rule matched any hyphenated name, so AL-RNA, ORF-O, ORF-P, Ori-S and gG-1 were all declared multi-gene spans |
| `uORF` | 0 | 3 | 0 | An upstream open reading frame, 5-prime of a main coding sequence |
| `uoORF` | 0 | 1 | 0 | An upstream open reading frame overlapping the main coding sequence |

## Rows by feature type

| Feature type | NC_001806.2 | BK012101.1 | whisnant2020-gtf |
|---|---|---|---|
| `CDS` | 77 | 300 | 309 |
| `exon` | 0 | 0 | 230 |
| `gene` | 79 | 83 | 170 |
| `mRNA` | 17 | 203 | 0 |
| `misc_feature` | 2 | 0 | 0 |
| `ncRNA` | 2 | 0 | 0 |
| `regulatory` | 72 | 0 | 0 |
| `repeat_region` | 31 | 0 | 0 |
| `stem_loop` | 3 | 0 | 0 |
| `transcript` | 0 | 0 | 201 |

## Rows by class, and the feature types in them

Stated so that nobody has to infer the mix from the counts above.

| Class | Source | Rows | Feature types |
|---|---|---|---|
| `dORF` | `BK012101.1` | 2 | CDS 1, gene 1 |
| `iORF` | `BK012101.1` | 2 | CDS 1, gene 1 |
| `lat_region` | `NC_001806.2` | 4 | gene 2, ncRNA 2 |
| `lat_region` | `BK012101.1` | 14 | CDS 4, gene 4, mRNA 6 |
| `lat_region` | `whisnant2020-gtf` | 18 | CDS 4, exon 5, gene 5, transcript 4 |
| `misc_feature` | `NC_001806.2` | 2 | misc_feature 2 |
| `named_gene` | `NC_001806.2` | 171 | CDS 77, gene 77, mRNA 17 |
| `named_gene` | `BK012101.1` | 30 | CDS 13, gene 6, mRNA 11 |
| `named_gene` | `whisnant2020-gtf` | 483 | CDS 166, exon 129, gene 77, transcript 111 |
| `novel_orf_2020` | `BK012101.1` | 53 | CDS 27, gene 18, mRNA 8 |
| `novel_orf_2020` | `whisnant2020-gtf` | 399 | CDS 138, exon 93, gene 85, transcript 83 |
| `origin_locus` | `BK012101.1` | 6 | CDS 2, gene 1, mRNA 3 |
| `origin_locus` | `whisnant2020-gtf` | 4 | CDS 1, exon 1, gene 1, transcript 1 |
| `primary_mirna` | `whisnant2020-gtf` | 6 | exon 2, gene 2, transcript 2 |
| `regulatory` | `NC_001806.2` | 72 | regulatory 72 |
| `repeat_region` | `NC_001806.2` | 31 | repeat_region 31 |
| `sORF` | `BK012101.1` | 25 | CDS 12, gene 11, mRNA 2 |
| `stem_loop` | `NC_001806.2` | 3 | stem_loop 3 |
| `transcript_span` | `BK012101.1` | 446 | CDS 236, gene 37, mRNA 173 |
| `uORF` | `BK012101.1` | 6 | CDS 3, gene 3 |
| `uoORF` | `BK012101.1` | 2 | CDS 1, gene 1 |

## Provenance

- **`NC_001806.2`**: RefSeq reference annotation, provisional, last modified 13 August 2018
- **`BK012101.1`**: Whisnant et al. 2020 re-annotation, third-party annotation on JN555585.1, Nat Commun 11:2038, PMID 32341360
- **`whisnant2020-gtf`**: Whisnant et al. 2020 annotation as distributed on Zenodo, 10.5281/zenodo.3724023
