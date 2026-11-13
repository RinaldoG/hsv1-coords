# Annotation v0.1.0: features by class and source

Produced by `src/build_annotation.py`. Every feature carries the class of thing it is and the record it came from.

**Why this table exists.** The two deposited annotations of this genome are often described as disagreeing, with 77 coding sequences against 300 on the identical 152,222 bases. They do not disagree. They describe different classes of feature, and a count is only comparable within a class. See `FINDINGS.md`.

| Class | NC_001806.2 | BK012101.1 | whisnant2020-gtf | What it is |
|---|---|---|---|---|
| `dORF` | 0 | 2 | 0 | A downstream open reading frame, 3-prime of a main coding sequence |
| `iORF` | 0 | 2 | 0 | An internal open reading frame, inside another coding sequence |
| `named_gene` | 215 | 30 | 479 | A classically named gene in the systematic positional scheme |
| `novel_orf_2020` | 2 | 53 | 403 | An open reading frame resolved by the 2020 re-annotation and numbered by position, with no independent name in the literature |
| `other` | 10 | 14 | 11 | Not matched by any class rule; see the notes column |
| `primary_mirna` | 0 | 0 | 6 | A primary microRNA transcript |
| `sORF` | 0 | 25 | 0 | A small open reading frame named by the 2020 re-annotation |
| `transcript_span` | 0 | 452 | 11 | A transcript unit spanning more than one named gene, so not a gene itself |
| `uORF` | 0 | 6 | 0 | An upstream open reading frame, 5-prime of a main coding sequence |
| `uoORF` | 0 | 2 | 0 | An upstream open reading frame overlapping the main coding sequence |

## Provenance

- **`NC_001806.2`**: RefSeq reference annotation, provisional, last modified 13 August 2018
- **`BK012101.1`**: Whisnant et al. 2020 re-annotation, third-party annotation on JN555585.1, Nat Commun 11:2038, PMID 32341360
- **`whisnant2020-gtf`**: Whisnant et al. 2020 annotation as distributed on Zenodo, 10.5281/zenodo.3724023
