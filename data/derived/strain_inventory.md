# Complete HSV-1 genomes: strain inventory

Fetched 2026-11-13 by `src/fetch_genomes.py`, NCBI nuccore query:

```
"Human alphaherpesvirus 1"[Organism] AND complete genome[Title] AND 140000:160000[SLEN]
```

- Records: **183**
- Distinct strain labels: **180**
- Strain labels deposited more than once: **2**
- Of those, labels whose records **disagree on length**: **1**

- Length range: 147,898 to 159,092 bp, a spread of 11,194 bp

## Same strain label, different length

These are the records a coordinate crosswalk has to reconcile. A position cited against one of them does not mean the same position in another.

### `KOS` - 2 records, spread 987 bp

| Accession | Length | Ambiguous bases |
|---|---|---|
| `JQ673480.1` | 152,011 | 0 |
| `JQ780693.1` | 151,024 | 0 |

