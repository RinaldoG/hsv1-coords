# Coordinate crosswalk: coverage

Every genome aligned to `NC_001806.2` with minimap2 (`asm10`, `--cs`, `--secondary=no`), reduced to ungapped blocks by `src/build_crosswalk.py`.

- Genomes mapped: **183**
- Median mappable fraction: **98.84%**, excluding the isomer-orientation cases below
- Range: 69.87% to 100.00%
- Genomes below 95% mappable: **8**

- Genomes with regions that align in more than one place: **24**, median 430 bp dropped where present

Those regions are dropped rather than resolved. They are the inverted repeats, where a position genuinely corresponds to two places in the genome, and picking one silently is the error this layer exists to prevent.

A position that falls outside every block has no reference coordinate, and the lift-over returns nothing for it rather than the nearest match. The unmappable fraction is concentrated in the inverted repeats, which is the region milestone 4 marks explicitly.

## Ten least mappable genomes

Excluding the isomer-orientation cases above.

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
