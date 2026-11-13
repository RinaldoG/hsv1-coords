# The repeat and diploid-locus map

Built by `src/build_repeat_map.py` from the reference record's own feature table. Nothing here is hand-listed.

- Repeated or structurally ambiguous regions: **31**
- Bases inside at least one of them: **32,019** of 152,222, which is 21.0 per cent of the genome
- Diploid locus pairs derived: **14**

## The four inverted repeats

| Repeat | Start | End | Length |
|---|---|---|---|
| `TRL` | 1 | 9,213 | 9,213 |
| `IRL` | 117,161 | 126,373 | 9,213 |
| `IRS` | 125,975 | 132,607 | 6,633 |
| `TRS` | 145,590 | 152,222 | 6,633 |

`TRL` and `IRL` are the same 9,213 bp of sequence in two copies, and `IRS` and `TRS` likewise at 6,633 bp. A position inside any of them names a place in the linear record but not which copy of the locus is meant.

## Diploid loci

Features present in both copies of a repeat, matched by feature type, equal length and offset within the repeat.

| Locus | Repeat pair | Type | Copy 1 | Copy 2 |
|---|---|---|---|---|
| `RL1` | TRL/IRL | CDS | 513..1,259 | 125,115..125,861 |
| `RL1` | TRL/IRL | gene | 513..1,540 | 124,834..125,861 |
| `RL1` | TRL/IRL | regulatory | 1,535..1,540 | 124,834..124,839 |
| `RL2` | TRL/IRL | gene | 2,087..5,699 | 120,675..124,287 |
| `RL2` | TRL/IRL | regulatory | 2,087..2,092 | 124,282..124,287 |
| `RL2` | TRL/IRL | CDS | 2,262..5,490 | 120,884..124,112 |
| `RL2` | TRL/IRL | regulatory | 5,677..5,682 | 120,692..120,697 |
| `regulatory@complement(7592..7597) / LAT` | TRL/IRL | regulatory | 7,592..7,597 | 118,777..118,782 |
| `LAT` | IRS/TRS | regulatory | 127,146..127,151 | 151,046..151,051 |
| `RS1` | IRS/TRS | regulatory | 127,187..127,192 | 151,005..151,010 |
| `RS1` | IRS/TRS | CDS | 127,235..131,131 | 147,066..150,962 |
| `RS1` | IRS/TRS | regulatory | 131,455..131,460 | 146,737..146,742 |
| `oriS palindrome` | IRS/TRS | stem_loop | 131,980..132,024 | 146,173..146,217 |
| `US1 / US12` | IRS/TRS | regulatory | 132,101..132,106 | 146,091..146,096 |

## Regions by kind

| Kind | Count |
|---|---|
| `repeat_region` | 28 |
| `stem_loop` | 3 |
