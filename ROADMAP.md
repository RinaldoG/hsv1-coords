# Roadmap

Each component ships only when it clears the acceptance criteria written below, which were
fixed before the code was written. A component that cannot clear its criteria is reported
as such rather than relaxed.

## M0: inputs and skeleton - done

Eight inputs fetched and checksummed into `data/raw/MANIFEST.tsv`, with a gate
(`make check`) that asserts the facts the rest of the work depends on rather than assuming
them. See `FINDINGS.md` for what the gate established, including one correction to the
premise this repository started from.

**Criteria, met:** every input fetches and verifies from a clean machine with one command.

## M1: gene name dictionary - met, with a recorded gap list

One row per gene: canonical identifier, systematic name (`UL48`), protein name (`VP16`),
and every historical alias (`alpha-TIF`, `Vmw65`), each alias carrying a citation.

The mechanical half is done: `src/extract_names.py` extracts every gene name from the
deposited records and reconciles them across sources, without inventing synonymy. Results
in `data/derived/`. The finding that shapes this milestone is in `FINDINGS.md`: of 232
distinct gene names across the four sources, five appear in all four.

**Criteria, all met** and enforced by `tests/test_m1_dictionary.py`, which runs in CI:

- Name to canonical to name is lossless for every entry.
- Reproduces the synonymy used in at least three published papers: 13 papers, 154
  assertions, generated from `data/derived/alias_testcases.md` by `make testcases` so the
  prose and the tests cannot drift.
- No entry carries an uncited alias: 683 rows, 166 genes, 103 citations, zero uncited.
- Where a name has been used for more than one locus, the collision is recorded rather
  than resolved, and `Dictionary.resolve` refuses it instead of guessing.

**Gaps are listed, not hidden.** `tests/known_gaps.tsv` holds 27 published name strings
the dictionary does not resolve, each with a category and a reason. A new gap fails the
suite, and so does a listed gap that has quietly started resolving, so the list cannot
rot. The largest category is the viral microRNAs, deferred because their reference
sequences need a sourcing decision of their own (see `FINDINGS.md`).

## M2: canonical versioned annotation - v0.1.0 released, GFF3 generation still to come

One semantically versioned annotation with a changelog and per-feature provenance, stating
which annotation class each feature belongs to, so that counts from differently scoped
annotations are not compared as though they measured the same thing.

**Criteria, two of three met** and enforced by `tests/test_m2_annotation.py`:

- **Met:** every feature carries a provenance record naming its source, plus the class of
  feature it is, plus its location. 1,723 rows across three sources.
- **Met, and reworded:** divergences are listed and justified in `CHANGELOG.md`. The
  criterion anticipated arbitrating between two conflicting annotations, and there is no
  conflict to arbitrate: one record is the other's successor on the same sequence with a
  wider scope. So it is met by recording that the difference is one of scope, and by making
  the class column stop the two counts being compared directly again. 452 of the 2020
  record's features are transcript units spanning more than one gene, and a further 37 are
  nested reading frames; the reference annotation describes none of either.
- **Outstanding:** parses with two independent GFF3 parsers, zero warnings. The current
  release is a classed feature table, not yet GFF3. GFF3 generation waits on the coordinate
  crosswalk (M3), because emitting one merged GFF3 across sources before the crosswalk
  exists would assert equivalences that cannot yet be checked.

## M3: strain coordinate crosswalk and lift-over - built, criteria met except the hand-checked published coordinates

Pairwise alignments across the complete genomes, chain files, and a lift-over tool, so that
a coordinate cited against one strain resolves against another.

**Inventory built** (`make genomes`): 183 complete genomes, 147,898 to 159,092 bp, a spread
of 11,194 bp. See `data/derived/strain_inventory.md`.

Two things the inventory establishes, both of which the crosswalk has to handle:

- **The laboratory strain KOS is deposited four times at four lengths** (152,011, 151,974,
  151,794 and 151,024 bp) under **three different label strings** (`KOS`, `isolate KOS`,
  `isolate KOS 1.1`), from different groups. A coordinate cited against "KOS" is ambiguous
  across up to 987 bp. Strain label normalisation is therefore part of this milestone, not
  a separate concern.
- **Corrected after review:** an earlier version of this file said the `ON783214` to
  `ON783217` set spans 3,421 bp "legitimately, because they are the four isomers". That
  was wrong twice over. Segment inversion preserves length, so a length spread is
  assembly disagreement and not isomerism. And the four records appeared under one strain
  label only because this repository's own regex truncated `B^3x1.1A` to `B`, which is
  now fixed. Three of those genomes turned out to be deposited on the opposite strand
  with a rotated start, which the build now detects and reverse complements; the one
  genuine segmental isomer in the set is `ON783215.1`, with its long segment forward and
  its short segment reversed.

**Dependency:** `minimap2` (tested against 2.31-r1302), preset `asm10`. Declared rather
than hand-rolled: a bespoke aligner is the wrong risk to take in a repository whose purpose
is that others can check the coordinates. `asm5` was used first and was too strict for this
divergence, costing real alignments.

**Built.** 183 genomes aligned to the reference and reduced to ungapped coordinate blocks
(`make crosswalk`), with `hsv1coords.lift` resolving a position between any two of them
through the reference. Median mappable fraction 98.08 per cent.

**Criteria:**
- **Met:** round-trip returns the original coordinate for at least 99 per cent of positions,
  enforced in `tests/test_m3_crosswalk.py` over 10,000 sampled positions with a fixed seed,
  failures enumerated rather than averaged.
- **Met:** returns no confident coordinate where it has none. Positions before the first
  block, past the last, and inside insertions relative to the reference all return nothing.
- **Outstanding:** resolving a published coordinate from at least five papers across at
  least three strains, each checked by hand. This needs the papers read, not the code run.

**Two things the build established, both recorded in `data/derived/crosswalk_report.md`:**

- **Regions that align in more than one place are dropped**, so the lift-over declines
  there. These are the inverted repeats, where a position genuinely corresponds to two
  places. An earlier version kept both and the milestone gate caught it.
- **Blocks are chained to be monotonic in both sequences.** Review found a position 900 bp
  into one genome resolving to 151,879, the far copy of the terminal repeat: a genuinely
  wrong coordinate, and exactly the error this layer exists to prevent. The build now
  keeps the largest chain that cannot jump backwards, and that position declines.
- **Genomes deposited on the opposite strand are reverse complemented** rather than
  reported as unmappable. Three were previously at 0 per cent and are now at 85.8, 94.7
  and 99.4 per cent. No genome in the set is unmapped.

**Worked example.** Between the two deposited genomes of strain KOS, the offset is not
constant: a position cited against `JQ673480.1` resolves in `JQ780693.1` at -122 bp near
the start, -482 bp at the LAT TATA box, and -626 bp by position 150,000. A fixed correction
would be wrong everywhere except by accident.

## M4: the repeat and diploid-locus map - built, with its limits stated

**The premise this was built on was wrong, and is corrected here.** An earlier review
told me the field masks the repeats with Ns and that no machine-readable description
exists. A later review checked both and neither holds.

- **Masking is one of at least three practices, and the named example does the opposite
  of what was claimed.** Whisnant et al. 2020 mask the terminal repeats and then
  explicitly restore the read mappings into the masked regions. Szpara and colleagues
  (*J Virol* 2014;88:1209-1227) trim to one copy of each inverted repeat and deposit both
  forms. Lassalle et al. (*Virus Evolution* 2020;6:veaa012) exclude by hard-coded
  coordinates, which disagree with the reference record at all three blocks.
- **Machine-readable descriptions do exist**, in at least five places: the RefSeq GFF3 for
  `GCF_000859985.2`, the BV-BRC API, the UCSC GenArk hub, several GitHub repositories, and
  the Whisnant GTF sitting in this repository's own `data/raw/`, which names every diploid
  locus in both copies. The README already credits that file as the best existing
  annotation, so the original framing contradicted this repository on the same commit.

**What does hold, and is the honest reason to ship this:** Dfam has no viral families,
RepeatMasker finds only simple repeats here, no purpose-built repeat interval file exists,
and **nobody publishes the counterpart arithmetic**. Knowing a region is repeated is
available everywhere. Being told which other position a coordinate could equally be is
not available anywhere.

What it is: every annotated repeated region of the genome as data, with loci present in
two copies paired to each other, and an interface answering **"is this position well
defined, and if not, which copies could it be?"**

```python
from hsv1coords import RepeatMap, lift_detailed

RepeatMap.load().classify(900)
# 900: not well defined (inside TRL, which is present in two copies); could equally be 125474

lift_detailed(118777, "JQ673480.1", "JQ780693.1")
# ... 118295, but NOT well defined: inside IRL, which is present in two copies.
#     The reference position could equally be 7487
```

**Built** (`make repeatmap`): 34 regions covering 32,019 bp, 21.0 per cent of the genome,
and 32 diploid locus pairs derived from the annotation rather than hand-listed, 4 of them
marked partial where the two copies differ in length.

**Criteria, met:**

- Reproduces all four inverted repeat boundaries from the reference's own feature table,
  exactly. Independently verified base by base.
- Every diploid pair is derived, not asserted, and both members have equal length.
- Machine-readable and enforced in code: `classify` returns a structured status.
- The lift-over consults it. `lift_detailed` reports whether a lifted coordinate
  identifies a unique locus and names the alternatives when it does not. **This was
  marked met before it was true**: `liftover.py` did not import the map at all, and the
  review caught it.
- Direct tests against the three defects review found in M3, including that 151,879, the
  cross-repeat coordinate that was being returned confidently, is a position the map knows
  is not well defined.
- A worked example shows a tool declining and naming the alternatives.

**Independently verified.** The claim that `TRL` and `IRL` are the same sequence inverted,
and likewise `IRS` and `TRS`, was checked base by base at every position in both pairs:
**15,846 of 15,846 reverse-complementary, no mismatches.** The alternatives are decisively
excluded, since a one-base offset or a direct-orientation reading both score 29 to 38 per
cent. `TRL` and `IRL` are byte-identical across all 9,213 bases, as are `IRS` and `TRS`
across 6,633.

**What the review corrected:**

- Three `repeat_region` features were being dropped because they carry a `/gene` qualifier
  and so are classed `named_gene` upstream. The region count of 31 matched the reference's
  31 repeat features by coincidence, which hid it.
- Excluding every `repeat_region` from pairing lost nine genuine pairs.
- The `'a'` sequence is present **three** times, not twice: `1..399` and
  `151824..152222` identical and in the same orientation, and `125975..126373` inverted
  relative to both. Verified base by base. Two successive versions of this map gave those
  positions an incomplete list of alternatives, first by modelling only inverted repeats
  and then by modelling only one direct pair.
- **LAT was not paired at all**, though this roadmap claimed it was. Only its 6 bp TATA
  boxes were. Containment-based matching could never find the gene, because LAT's
  long-repeat copy runs 778 bp past the end of `IRL`. Pairing now allows overlap and
  marks the result partial.
- The docstring example gave the wrong counterpart for its own worked position.
- `classify` asserted that a position outside the large repeats is "still unique in the
  genome". That is false for every tandem reiteration set, and the map cannot certify
  uniqueness in any case: it describes annotated repeat structure and now says so.
- Out-of-range coordinates were reported as well defined.
- `IRL` and `IRS` genuinely overlap across `125,975..126,373`, so a position there belongs
  to both and has a counterpart in each. Found by the gate during the build; verified
  correct, and no other pair overlaps.

**Known limits, stated rather than fixed:**

- The map covers the reference only. A coordinate in another strain is related to it
  through the crosswalk, not described directly.
- Two further diploid pairs are missed because the reference annotates them with a 1 bp
  asymmetry.
- Output is TSV only. A specialist would want BED and GFF3, and that is cheap to add.
- **The highest-value thing this does not yet do**, from the usefulness review: resolve a
  coordinate between the different conventions the field actually uses for these regions,
  full-length against Szpara-trimmed against Lassalle-excluded. That, not the repeat
  intervals, is the real reason an HSV-1 coordinate cannot be trusted, and M3's lift-over
  is a genuine advantage there rather than a reimplementation of something that exists.
- `classify` returns no strand, has no interval API, and does not label which repeat each
  alternative belongs to.
- Tandem reiteration sets are flagged as repeat features but no alternative positions are
  computed for them, because a tandem array has no single counterpart.

## M5: the coordinate-convention resolver - built, criteria met

Different groups handle the repeated regions differently, so the same printed coordinate
means different things in different papers, and nothing reconciles them. This holds each
handling as data, with its definition quoted from its own paper, and answers whether a
position is present and what its equivalent coordinate is.

```python
ConventionSet.load().report(118777)      # the LAT TATA box
# usable under full_length
# usable under whisnant2020_masked
# usable under szpara2014_trimmed
# NOT usable under lassalle2020_excluded: inside 116,926-132,647, excluded from the analysis
# NOT usable under morse2017_nr: inside 117,161-126,373
# caveat: ... not well defined in the reference: inside IRL ... could equally be 7597
```

**Five handlings, each read from its own source**, with the date it was read recorded:

| Name | What the source says it does | Coordinates |
|---|---|---|
| `full_length` | The deposited record, 152,222 bp | unchanged |
| `lassalle2020_excluded` | Excludes `1-9,338`, `116,926-132,647`, `145,582-150,962` from analysis | unchanged |
| `whisnant2020_masked` | Masks the two terminal repeats, `TRL` and `TRS`, then copies the read mappings back in | unchanged |
| `szpara2014_trimmed` | Keeps one copy of each inverted repeat, `IRL` and `IRS`, removing `TRL` and `TRS` | **shifted** |
| `morse2017_nr` | Removes the `IRL` and `TRS` blocks | **shifted**, and defined on `NC_001806.1`, so refused here |

**The worked example.** The LAT TATA box at reference 118,777 is printed as **109,564**
in a Szpara-trimmed analysis, **excluded** by Lassalle, and **retained** by Whisnant.
Three published handlings, three different answers, for one position in the gene that
dominates latent transcription.

**The clearest illustration is a pair with identical intervals.** Szpara and Whisnant
both act on `TRL` and `TRS`. Szpara removes them, so every downstream coordinate shifts.
Whisnant masks them and copies the read mappings back, so no coordinate moves. Same two
intervals, opposite consequences for a printed number.

**And a second problem, orthogonal to the repeats, which review surfaced.**
`NC_001806.1` is 152,261 bp and `NC_001806.2` is 152,222 bp, and **the offset between
them is not constant**: roughly 0 near the start, -1 through `UL` and `IRL`, -3 to -5
near `IRS`, and +39 through most of `TRS`. Morse's NR-reference is built on `.1`. So a
coordinate can be wrong simply because two papers used two versions of the same
accession, with no repeat involved at all, and nothing in this repository could detect
it until now. `morse2017_nr` is therefore recorded with its definition verified and
**refused for conversion**, because applying it in `.2` space would give confident wrong
answers.

**Criteria, written before the code, all met** and enforced in
`tests/test_m5_conventions.py`:

- Every definition is quoted from its own primary source, with the source and the date
  it was read, or the convention is marked unverified and refused.
- A convention that only drops data leaves coordinates unchanged, and "unusable here" is
  a distinct answer from "moves to a different number".
- A convention that removes sequence shifts coordinates, and the shift is **computed from
  the retained blocks, not assumed constant**. Morse removes an internal block, so a
  single offset cannot describe it, and a test asserts that.
- Round trip: every retained position converts out and back to itself, in both collapsed
  references.
- The resolver states which handlings a position is usable under.
- A worked example reconciles one locus across the handlings, above.

**What review corrected, and it is the important part:**

- **Two of the five intervals were wrong.** `whisnant2020_masked` was encoded as `TRL`
  plus `IRL`, which is both copies of the long repeat and neither of the short. That is
  impossible, because the paper restores read mappings from the unmasked copy and this
  encoding leaves none. The sentence immediately before the quoted one settles it: "The
  HSV-1 genome consists of two components (L and S) that are both flanked by long
  repeats." It masks `TRL` and `TRS`. `morse2017_nr` over-removed by 399 bases, because
  their `IRL` block stops before the a-sequence that belongs to the `IRS` block they
  kept.
- **The reference-version problem above was invisible to every check here.**

- **`szpara2014_trimmed` was marked unverified and refused, on the strength of one
  failed API call.** The paper is free to read on the publisher's own site, and another
  review in this same project had quoted its methods verbatim the same day. A failed
  request is not a paywall. It is now read, verified and implemented.
- **A second collapsed reference existed and was missed**: Morse et al. 2017. It removes
  a different pair of blocks from Szpara, which is what makes the disagreement real
  rather than theoretical.
- **The criteria said "all met" while the coordinate-shifting criterion was unimplemented
  and listed as a known limit in the same section.** That is the same error M4 was
  corrected for, repeated one milestone later. It is now actually implemented.
- **A gap in Lassalle's exclusion list was being reported as comparable.** Their list
  stops at 150,962 while `TRS` runs to 152,222, so positions in that tail were called
  usable without comment. `report` now consults the repeat map and adds a caveat.
- `comparable_everywhere` counted refused conventions, making it False at every position
  in the genome.

**Known limits:**

- Five handlings, one of them refused. There are certainly more: review documented five
  further distinct ones, of which reference-fill and synthetic duplication matter most,
  because both produce records that look like complete 152 kb data and are not.
- **No crosswalk between `NC_001806.1` and `.2` exists here**, which is what blocks
  Morse. Building one is the obvious next piece and would also let the repo detect a
  class of error it currently cannot see.
- `lassalle2020_excluded` is transcribed verbatim, but review found it is one of at
  least three frames in that paper and not the one behind its variant calls. Its numbers
  are gene boundaries. The row is correct as quoted and narrower than its name suggests.
- "Convention" is this repository's word, not the field's. Papers say trimmed,
  NR-reference, masked, excluded.
- The resolver reasons about the reference only. Other strains reach it through M3.

## M5b: release

Versioned release with a dataset DOI, a container image, and a README from which someone who
has not seen the project can reproduce a release.

**Criteria:** every figure, table and number in the repository regenerates from one command,
tested by someone who has not seen it before.

## Under consideration, not committed

A repeat-aware pangenome graph over the complete genomes and the isomer-resolved set. This
becomes tractable once M3 and M4 exist and not before. The method is established for other
large DNA viruses; no published graph representation of a herpesvirus genome was found as of
November 2026.

## Contributing

Issues and corrections are welcome, particularly: a gene alias with a citation that is
missing from the dictionary, a published coordinate that the crosswalk resolves incorrectly,
and any case where the mask is wrong in either direction.
