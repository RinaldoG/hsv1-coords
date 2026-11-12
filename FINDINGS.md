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
