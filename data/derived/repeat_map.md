# The repeat and diploid-locus map

Built by `src/build_repeat_map.py` from the reference record's own feature table. Nothing here is hand-listed.

- Repeated or structurally ambiguous regions: **34**
- Bases inside at least one of them: **32,019** of 152,222, which is 21.0 per cent of the genome
- Diploid locus pairs derived: **32**

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
| `'a' sequence` | TRL/IRL | repeat_region | 1..399 | 125,975..126,373 |
| `LAT` | TRL/IRL | gene | 1..7,569 | 118,777..127,151 |
| `LAT` | TRL/IRL | ncRNA | 1..7,569 | 118,805..127,151 |
| `'a' sequence reiteration set` | TRL/IRL | repeat_region | 98..320 | 126,054..126,276 |
| `RL1` | TRL/IRL | CDS | 513..1,259 | 125,115..125,861 |
| `RL1` | TRL/IRL | gene | 513..1,540 | 124,834..125,861 |
| `RL1 / RL reiteration set 6` | TRL/IRL | repeat_region | 988..1,040 | 125,334..125,386 |
| `RL1` | TRL/IRL | regulatory | 1,535..1,540 | 124,834..124,839 |
| `RL2` | TRL/IRL | gene | 2,087..5,699 | 120,675..124,287 |
| `RL2` | TRL/IRL | regulatory | 2,087..2,092 | 124,282..124,287 |
| `RL2` | TRL/IRL | mRNA | 2,114..5,699 | 120,675..124,287 |
| `RL2` | TRL/IRL | CDS | 2,262..5,490 | 120,884..124,112 |
| `RL2 / RL reiteration set 5` | TRL/IRL | repeat_region | 2,466..2,624 | 123,750..123,908 |
| `RL2` | TRL/IRL | regulatory | 5,677..5,682 | 120,692..120,697 |
| `RL reiteration set 4` | TRL/IRL | repeat_region | 5,732..5,878 | 120,496..120,642 |
| `regulatory@complement(7592..7597) / LAT` | TRL/IRL | regulatory | 7,592..7,597 | 118,777..118,782 |
| `RL reiteration set 3` | TRL/IRL | repeat_region | 8,554..8,585 | 117,789..117,820 |
| `RL reiteration set 2` | TRL/IRL | repeat_region | 8,806..8,826 | 117,548..117,568 |
| `RL reiteration set 1` | TRL/IRL | repeat_region | 9,033..9,213 | 117,160..117,341 |
| `'a' sequence` | IRS/TRS | repeat_region | 125,975..126,373 | 151,824..152,222 |
| `'a' sequence reiteration set` | IRS/TRS | repeat_region | 126,054..126,276 | 151,921..152,143 |
| `RS reiteration set 1` | IRS/TRS | repeat_region | 126,573..126,712 | 151,486..151,624 |
| `RS reiteration set 2` | IRS/TRS | repeat_region | 126,813..127,145 | 151,052..151,384 |
| `LAT` | IRS/TRS | regulatory | 127,146..127,151 | 151,046..151,051 |
| `RS1` | IRS/TRS | gene | 127,173..131,460 | 146,737..151,023 |
| `RS1` | IRS/TRS | mRNA | 127,173..131,431 | 146,766..151,023 |
| `RS1` | IRS/TRS | regulatory | 127,187..127,192 | 151,005..151,010 |
| `RS1` | IRS/TRS | CDS | 127,235..131,131 | 147,066..150,962 |
| `RS1` | IRS/TRS | regulatory | 131,455..131,460 | 146,737..146,742 |
| `oriS palindrome` | IRS/TRS | stem_loop | 131,980..132,024 | 146,173..146,217 |
| `US1 / US12` | IRS/TRS | regulatory | 132,101..132,106 | 146,091..146,096 |
| `US1 / RS reiteration set 3` | IRS/TRS | repeat_region | 132,391..132,516 | 145,681..145,806 |

## Regions by kind

| Kind | Count |
|---|---|
| `repeat_region` | 31 |
| `stem_loop` | 3 |
