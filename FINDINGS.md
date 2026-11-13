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
  feature is incomplete at that end. **Corrected 13 November 2026:** an earlier version of
  this file said the feature was clipped by "where the linear record happens to begin".
  That is the wrong reason. Position 1 is the L terminus, and the missing portion is
  present in the same record inside `TRS` at the other end, because the genome is
  circular. The truncation is a consequence of linearising a circle, not of an arbitrary
  start point.
- Second copy: `gene 118777..127151`, 8,375 bp, with a TATA box at `118777..118782` and
  `ncRNA join(118805..119464,121421..127151)`.

So the full-length annotation of the gene that dominates latent transcription exists only
in the second copy. Any analysis using the first copy's coordinates is working with a
truncated feature, and nothing in the file stops it.

**Also corrected:** LAT does not lie wholly inside the long inverted repeat. Copy 2 ends
at 127,151 while `IRL` ends at 126,373, so it extends 778 bp into `IRS`. Earlier versions
of this file and the README asserted otherwise. The overhang is the reason the `TRL` copy
is flagged incomplete, so the two facts are the same fact.

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

## A name in an inverted repeat is not an ambiguous name

Found 13 November 2026 by the milestone 1 gate, which is the point of writing the gate
before the code.

Three alias strings in the dictionary map to two canonical genes each. The gate initially
treated all three as ambiguous. Two of them are not:

| Alias | Maps to | What it is |
|---|---|---|
| `LAT 0.7kb` | `TRL2.5`, `IRL2.5` | One locus, present in both copies of the long inverted repeat |
| `pri-miRNA-H1` | `TRL2.5`, `IRL2.5` | The same |
| `IE63` | `UL54`, `US1` | A genuine collision across two different genes |

`TRL2.5` and `IRL2.5` are the same locus. It appears twice because it lies inside the
inverted repeat, which is also why LAT itself appears twice. So a name mapping to both is
not under-determined in the way an ambiguous name is: it is fully determined apart from
which repeat copy is meant, and that is a coordinate question, not a naming one.

The resolver now separates them. A repeat-paired name raises `RepeatPairedName`, which
tells the caller it must choose a copy or accept both. A genuinely ambiguous name raises
`AmbiguousName`. Collapsing the two would have meant either silently picking one copy of
LAT, which is the error the un-callable mask exists to prevent, or declaring a
well-determined locus unusable.

`IE63` is the real collision, and the dictionary documents it rather than resolving it:
`IE63` on `UL54` (ICP27) is well attested in HSV-1 papers, while `IE63` on `US1` (ICP22)
traces to a varicella-zoster virus paper on that virus's ORF63, which is the likely source
of the confusion. Both usages are recorded with citations and neither is deleted.

**One correction to the gate itself.** The first version of the test demanded that every
row for an ambiguous alias be marked `ambiguous`. That was wrong, and the data was right:
`confidence` describes how well sourced one gene-and-alias pair is, which is a different
axis from whether the alias collides across genes. Flattening every `IE63` row would have
thrown away the fact that one of the two assignments is far better supported than the
other. The test now requires the collision to be documented in at least one row, with an
explanation, which is what the data already did.


## Corrections from the review of 13 November 2026

Three reviewers were run over this repository before the next component was built: one on
whether the lift-over is actually correct, one on the engineering, one on whether a
specialist would accept the claims. Between them they found the following, all now fixed,
and all recorded here rather than quietly corrected.

| What was wrong | What is true |
|---|---|
| The annotation dropped every feature without a `/gene` qualifier, omitting all four inverted repeats, both `'a'` sequences and the `oriL` palindrome, while the changelog said nothing was dropped | They are carried, with their own feature classes. 1,779 rows, up from 1,723 |
| The GenBank parser dropped wrapped location lines, shipping 8 truncated coordinates including `RL2` (ICP0) | Fixed, and the gate now checks brackets balance and coordinates parse |
| "LAT lies wholly inside the inverted repeats" | It extends 778 bp past `IRL` into `IRS`, which is why the other copy is truncated |
| The truncated copy is clipped "by where the linear record happens to begin" | It is clipped by linearising a circular genome; the remainder is in `TRS` at the other end |
| "No machine-readable dictionary says `UL48`, `VP16` and `alpha-TIF` are one gene" | UniProt `P06492` does, and this repository's own alias table cites it |
| Four records of strain "B" spanning 3,421 bp, described as legitimate isomer variation | Inversion preserves length, so that spread is assembly disagreement. The single label was this repository's regex truncating `B^3x1.1A` |
| Three genomes "in a different isomer orientation", unmappable | They are opposite-strand deposits with a rotated start. Reverse complemented, they map at 85.8, 94.7 and 99.4 per cent. The one real segmental isomer is `ON783215.1` |
| A coordinate cited against KOS is "ambiguous across 987 bp" | That was a length difference, not a coordinate offset. The measured offset between the two plain-KOS records runs from -79 to -887 bp depending on position. Two of the four KOS deposits are a deliberate KOS and KOS1.1 pair with published provenance markers (Colgrove et al., *Virology* 2016;487:215-221) |
| `asm5` as the alignment preset | Too strict for 0.7 to 1.3 per cent divergence. `asm10` raises the median mappable fraction from 98.08 to 98.84 per cent |
| A position 900 bp into `OP297884.1` resolved to 151,879 | A cross-repeat placement, and a genuinely wrong coordinate. Blocks are now chained to be monotonic, and that position declines |
| The round-trip gate could not fail, being a mathematical identity over the same blocks | The gate now checks the base at the lifted position against the reference, with a deliberate one-base offset as a control proving the check has power |
| Continuous integration appeared to verify the build and did not: the manifest check was circular, the diff gate saw 6 of 196 derived files, and `minimap2` was never installed | All three closed |

The pattern worth noting: **every one of these was found by someone checking against the
primary source rather than against this repository's own tests.** The tests passed
throughout.
