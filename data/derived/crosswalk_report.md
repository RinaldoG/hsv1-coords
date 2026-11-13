# Coordinate crosswalk: coverage

Every genome aligned to `NC_001806.2` with minimap2 (`asm5`, `--cs`, `--secondary=no`), reduced to ungapped blocks by `src/build_crosswalk.py`.

- Genomes mapped: **183**
- Median mappable fraction: **98.08%**, excluding the isomer-orientation cases below
- Range: 67.86% to 100.00%
- Genomes below 95% mappable: **13**

- Genomes with regions that align in more than one place: **36**, median 312 bp dropped where present

Those regions are dropped rather than resolved. They are the inverted repeats, where a position genuinely corresponds to two places in the genome, and picking one silently is the error this layer exists to prevent.

A position that falls outside every block has no reference coordinate, and the lift-over returns nothing for it rather than the nearest match. The unmappable fraction is concentrated in the inverted repeats, which is the region milestone 4 marks explicitly.

## Genomes in a different isomer orientation

**3 genomes map at 0 per cent while having alignments.** These are not poor assemblies. The HSV-1 genome exists as four isomers, because the long and short segments invert independently, and a genome deposited in a different arrangement from the reference has its alignments on the reverse strand, which this build discards.

Handling them needs per-segment orientation, which belongs with the isomer-aware mask in milestone 4. Until then the crosswalk declines for them entirely, which is correct but incomplete.

| Accession | Length | Alignments found | Blocks kept |
|---|---|---|---|
| `ON783216.1` | 151,352 | 43 | 0 |
| `ON783217.1` | 154,773 | 51 | 0 |
| `PQ569920.1` | 151,817 | 5 | 0 |

## Ten least mappable genomes

Excluding the isomer-orientation cases above.

| Accession | Length | Mapped | Blocks |
|---|---|---|---|
| `ON783215.1` | 153,542 | 67.86% | 269 |
| `ON783214.1` | 152,108 | 79.19% | 293 |
| `OR405519.1` | 157,769 | 91.62% | 132 |
| `AB618031.1` | 149,643 | 91.84% | 88 |
| `MH999849.1` | 155,011 | 93.50% | 144 |
| `ON960060.1` | 152,428 | 93.52% | 138 |
| `OP297861.1` | 149,598 | 93.56% | 107 |
| `OP871843.1` | 159,092 | 94.38% | 146 |
| `OP871859.1` | 156,328 | 94.91% | 136 |
| `ON960056.1` | 152,443 | 94.98% | 82 |
