# Findings from the inputs themselves

Facts established by running `make check` against the deposited records, not taken from a
paper or a search result. Each is asserted in `src/check_inputs.py`, so if one stops
holding, continuous integration fails.

Checked 13 November 2026.

## `NC_001806.2` and `JN555585.1` hold the identical sequence

Both are 152,222 bp and the bases are byte-for-byte identical. The Whisnant et al. 2020
re-annotation uses `JN555585` as its sequence name while the RefSeq reference is
`NC_001806.2`, which reads like a coordinate-base mismatch and is not one. **It is an
accession naming difference, not a coordinate difference.**

This materially simplifies the plan: milestone 2 does not have to reconcile two coordinate
systems before it can use the 2020 annotation, and milestone 3's lift-over does not need an
entry for this pair. It also means anyone who has avoided the 2020 annotation because its
sequence name differs from the reference has avoided it for no reason, which is worth
stating plainly in the release.

## The 77 against 300 CDS discrepancy is real and is not a length difference

`NC_001806.2` carries 77 CDS features and 79 gene features. `BK012101.1` carries 300 CDS
features and 83 gene features. Both describe 152,222 bp. So the two records disagree about
what is coding on identical sequence, which is a feature-level disagreement, not an
assembly one. Milestone 2 documents both rather than silently choosing.

## LAT is annotated twice, and one copy is truncated by the sequence boundary

In `NC_001806.2`:

- First copy: `gene complement(<1..7569)`, with
  `ncRNA complement(join(<1..4953,6910..7569))`. The `<` is the record's own flag that the
  feature runs off the start of the sequence.
- Second copy: `gene 118777..127151`, 8,375 bp, with a TATA box at `118777..118782` and
  `ncRNA join(118805..119464,121421..127151)`.

So the full-length annotation of the gene that dominates latent transcription exists only
in the second copy. The first is clipped by where the linear record happens to begin. Any
analysis using the first copy's coordinates is working with a truncated feature, and
nothing in the file stops it.

## `TRL` is annotated at `1..9213`

Confirmed directly. The reference also carries 31 `repeat_region` features in total,
including nested reiteration sets inside `TRL` (the 'a' sequence at `1..399`, a reiteration
set at `98..320`, and further sets at `988..1040`, `2466..2624`, `5732..5878` and
`8554..8585`). The mask in milestone 4 has more structure available to it than the four
large repeats alone, and should use it.

## Correction: `BK012101.1` is the 2020 re-annotation, not a competing record

Established 13 November 2026 by reading the record header, which the research strands had
not done.

`BK012101.1` is a **Third Party Annotation** deposit by Whisnant, Juerges, Erhard and
Doelken (submitted 7 November 2019, released 7 November 2020), citing Whisnant et al.,
"Integrative functional genomics decodes herpes simplex virus 1", *Nature Communications*
11(1):2038, 2020, PMID 32341360. Its `COMMENT` block states outright that it uses data
from `JN555585.1`, which this repository has separately confirmed is byte-identical to
`NC_001806.2`.

**So the "77 CDS against 300 CDS on identical bases" is not two laboratories disagreeing.**
It is the 2018 RefSeq annotation against the 2020 re-annotation, correctly deposited as TPA
on the same sequence and properly linked to its paper. The earlier framing, carried from
the research strands, overstated this. The 2020 work is in GenBank, discoverable and
attributed.

The 300 CDS features are also not 300 genes. Across the record's 83 gene names:

| Kind of name | Count |
|---|---|
| Hyphenated span, e.g. `UL11-UL14`, `TRL1-TRL2` | 42 |
| Plain gene name | 24 |
| `sORF` (small open reading frame) | 11 |
| `uORF` (upstream) | 3 |
| `iORF`, `uoORF`, `dORF` (internal, upstream-overlapping, downstream) | 1 each |

So the record annotates nested and alternative reading frames and multi-gene transcript
units that the RefSeq annotation does not attempt. The two records describe different
**classes** of feature, and the counts are not in contradiction.

## What the naming problem actually is

The premise needed correcting, but the naming problem is real and is sharper than stated.

Of **232 distinct gene names** across the four sources, only **5 appear in all four**:
`UL22`, `UL36`, `UL41`, `UL43`, `UL48`. 135 names appear in exactly one source. The
`product` text differs between sources for **92** genes.

The clearest case, and the one that makes the project's case on its own: **the same group's
own two artefacts use different naming for the same features.** The Whisnant 2020 GenBank
TPA deposit carries 83 gene names, the Whisnant 2020 Zenodo GTF carries 170, and they share
**27**. The GTF names novel open reading frames individually (`UL5.5`, `UL11.5`,
`IRL1.5`), while the TPA deposit names transcript units spanning them (`UL11-UL14`,
`TRL1-TRL2`, `IRL2-IRL1`). `LAT_gene` in the GTF is `LAT-ORF-O` in the TPA record.

Neither is wrong. They are describing the same biology at two levels, for two different
submission formats, and nothing anywhere states the mapping. A reader with both files in
hand, from one paper, cannot join them.

**Consequence for the plan.** Milestone 2 shrinks: there is no annotation to arbitrate,
because the good annotation already exists and is deposited twice. What is missing is the
**crosswalk between naming systems**, including between one paper's own two deposits. That
is milestone 1's job, and it has just become the centre of the project rather than its
warm-up.
