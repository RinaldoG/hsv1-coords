# Complete HSV-1 genomes: strain inventory

Fetched 2026-11-14 by `src/fetch_genomes.py`, NCBI nuccore query:

```
"Human alphaherpesvirus 1"[Organism] AND complete genome[Title] AND 140000:160000[SLEN]
```

- Records: **183**
- Distinct strain labels, as deposited: **180**
- Distinct laboratory strains, after dropping qualifiers such as `isolate` and `clone`: **177**
- Laboratory strains deposited more than once: **3**
- Of those, strains whose records **disagree on length**: **2**

The two counts differ because a depositor may write `KOS`, `isolate KOS` and `isolate KOS 1.1` for records of the same laboratory strain. Grouping on the label as deposited splits those three ways, which is what this file did until round 11, and it reported one disputed strain while missing the KOS family entirely. The `strain` column in the TSV keeps the label as deposited; `strain_group` is the unit counted here.

The grouping drops leading qualifiers only. A hyphenated construct such as `McKrae-dsncRNA2` keeps its own identity rather than joining `McKrae`, because merging constructs with their parents would also merge isolate designations such as `HSV1/USA/WA-UW-2L9/2020`, which are not laboratory strains. So these counts are a floor on how much the same material has been deposited at different lengths, not a complete reckoning.

- Length range: 147,898 to 159,092 bp, a spread of 11,194 bp

## Same strain label, different length

These are the records a coordinate crosswalk has to reconcile. A position cited against one of them does not mean the same position in another.

### `KOS` - 4 records, spread 987 bp

| Accession | Length | Ambiguous bases |
|---|---|---|
| `JQ673480.1` | 152,011 | 0 |
| `KT899744.1` | 151,974 | 0 |
| `KT887224.1` | 151,794 | 0 |
| `JQ780693.1` | 151,024 | 0 |

### `MACINTYRE` - 2 records, spread 659 bp

| Accession | Length | Ambiguous bases |
|---|---|---|
| `MN136523.1` | 152,476 | 0 |
| `PQ569920.1` | 151,817 | 0 |

