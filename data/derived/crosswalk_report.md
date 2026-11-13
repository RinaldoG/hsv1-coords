# Coordinate crosswalk: coverage

Every genome aligned to `NC_001806.2` with minimap2 (`asm10`, `--cs`, `--secondary=no`), reduced to ungapped blocks by `src/build_crosswalk.py`.

- Genomes mapped: **183**
- Median mappable fraction: **98.84%**
- Range: 69.87% to 100.00%
- Genomes below 95% mappable: **8**

- Genomes with regions this alignment placed in more than one way: **22**, median 205 bp removed where present
- Genomes with blocks dropped by chaining: **4**, 78,667 bp in total

**What this does not do.** minimap2 runs with `--secondary=no`, so a second placement of the same query region is discarded by the aligner before the ambiguity filter sees it. The counts above are therefore about regions a single alignment pass placed twice, and they are not a measure of repeat ambiguity: round 7 measured roughly 92 per cent of positions inside the inverted repeats still receiving one confident block. Repeat ambiguity is handled by `liftover.lift` consulting the repeat map, which declines where a position has more than one valid copy, and not by this build.

A position that falls outside every block has no reference coordinate, and the lift-over returns nothing for it rather than the nearest match.

## Ten least mappable genomes

Excluding any genome that has alignments but maps at 0 per cent, listed above where such genomes exist.

| Accession | Length | Mapped | Blocks |
|---|---|---|---|
| `PQ569920.1` | 151,817 | 69.87% | 112 |
| `ON783216.1` | 151,352 | 79.40% | 745 |
| `ON783217.1` | 154,773 | 79.67% | 739 |
| `ON783215.1` | 153,542 | 80.62% | 840 |
| `OR405519.1` | 157,769 | 92.09% | 154 |
| `AB618031.1` | 149,643 | 92.70% | 105 |
| `ON960060.1` | 152,428 | 94.15% | 156 |
| `ON783214.1` | 152,108 | 94.31% | 1088 |
| `OP297884.1` | 156,387 | 95.07% | 155 |
| `OP871843.1` | 159,092 | 95.18% | 158 |
