# Changelog

Semantic versioning. The annotation is versioned because the thing being fixed is partly
that the reference has no versioning discipline.

## v0.2.0 - 17 November 2026

The first tagged release. Every acceptance criterion written in `research/PICK.md`
before any code was written is met, the last of them being the one a test run cannot
settle on its own: published coordinates resolved by hand out of the papers.

**Why 0.2.0 and not 0.1.0.** `v0.1.0` appears above as the annotation artefact's data
version, which is what `annotation_v0.1.0.tsv` is named for, and it was also the package
version. It was never a tag: until today this repository had no tag at all. There is
real added behaviour since that version string was set, so the first tag takes the next
minor number rather than reusing one that already means something else in this project.
Still `0.x`, so nothing here is a stability promise, and the annotation table keeps its
own `v0.1.0` data version because that artefact has not changed.

### The by-hand criterion closed, and the two calls a reader asked for

**Added**

- `data/published_coordinates.tsv` and `data/published_coordinate_checks.tsv`: the
  record of M3's by-hand criterion, one row per paper opened and one row per mechanical
  check. Rendered as `data/derived/published_coordinates.md` by
  `make published-coordinates` and gated by `tests/test_m3_published_coordinates.py`.
  Six published coordinates resolved by hand across four strains, against a threshold
  of five and three. The published coordinate, the quote and the fetch route in each
  row are hand transcription, flagged in a `transcribed` column; everything downstream
  of the transcription is checked against the records and the artefacts, and no
  expected value was produced by calling the function it checks.
- `hsv1coords.features_at`, `introns_at` and `describe`: what is at a reference
  position. An intron is answered separately from an exon, and `describe` leads with
  the repeat ambiguity when a position names more than one place.
- `hsv1coords.lift_range` and `RangeLift`: an interval lifted as an interval, with
  `length_change` as a field and in the printed form whenever it is not zero. Lifting
  two endpoints and subtracting them gives the wrong length when an indel lies between
  them, and two independent reviewers hit that on two different papers.
- `hsv1coords.AnnotatedFeature` and `hsv1coords.sources()`, the types behind the above.

**Changed, and a caller may see different text**

- The `paper_transcribed` caveat printed by `ConventionSet.report` no longer says the
  numbers were transcribed "from a paper this repository does not hold". That
  contradicted the `source` column of the one row carrying the basis, which records the
  full text as read. Two facts had been conflated and only one is what the basis
  records.

**Not changed, and listed because it was asked for**

- `lift()` still returns a bare `None` where it declines, with the explanation only in
  `lift_detailed`. All three reviewers raised it. Changing what `lift` returns is an
  API break and belongs in a release rather than a patch; `ROADMAP.md` records it with
  the rest of what the by-hand work asked for and did not get.

**A note on versions, since this file claims semantic versioning.** There is no git tag
in this repository at all. `v0.1.0` above is the annotation artefact's data version,
which is what `annotation_v0.1.0.tsv` is named for, and the package version in
`pyproject.toml`. Nothing has been released, so the sections above are what a first
release would contain.

### The annotation as GFF3, and the repeat map as BED

**Added**

- `data/derived/annotation_v0.1.0_<source>.gff3`, one file per annotation source
  (`NC_001806.2`, `BK012101.1`, `whisnant2020-gtf`), each under that source's own
  sequence name, rendered from the v0.1.0 table by `make gff3` (`src/build_gff3.py`).
  The table is unchanged: these are its rows in a format genome browsers and GFF3 tools
  read. How joins, partial ends, phase and parents are represented is stated in the
  builder's docstring.
- `data/derived/repeat_map.bed` (BED6+2) and `data/derived/diploid_pairs.bedpe`
  (BEDPE): the repeat map and the diploid pairs in BED's 0-based, half-open frame.
- `tests/test_m2_gff3.py`: each table row in its GFF3 exactly once, compared as whole
  tuples; every CDS phase checked by translating against the record's `/translation`;
  the BED frame checked by reading bases.
- `make validate-gff3` (`src/validate_gff3.py`), also run in CI: GenomeTools
  `gt gff3validator`, bcbio-gff and the `gff3` validator. They are validation tools
  only; the library still needs nothing beyond the standard library.

**Changed, with no change to any artefact**

- `build_annotation.annotation_rows()` returns each table row with the source feature
  it came from, and `extract_names.parse_gtf` also keeps the GTF's sequence name, phase
  and `orf_id`. The table and its summary rebuild byte-identical.

### Library behaviour changed in rounds 16 and 17

Round 17: the library's own API changed in four ways across rounds 16 and 17 and none
of them was recorded here, in a file whose second line claims semantic versioning. The
package version stays `0.1.0` because nothing is released yet; these are the changes a
caller would notice if it were.

**Changed, and a caller may see different answers**

- `Convention.from_record(record, offset)` returns a coordinate in NC_001806.2, this
  library's base version, where it used to return the convention's own NC_001806.1
  position under a docstring promising "back to the reference". `native=True` gives the
  untranslated one. 109,639 positions disagreed with `ConventionSet.coordinate_in`
  before this, out of 152,222 bp of sequence.
- `Convention.locate`, `coordinate` and `from_record` refuse a base-version position on
  a convention defined against another reference version, naming the `ConventionSet`
  method that translates, unless passed `native=True`.
- `Convention.to_reference` and `to_convention` return None for a position outside the
  convention's own range, where they used to compute one from it.
- `Convention.from_record` raises `ValueError` for an offset past the end of a record,
  where it returned None, which it still returns for the seven positions that exist in
  the convention's own version and have no NC_001806.2 equivalent.
- `Crosswalk.load` raises `WrongReferenceVersion` for `NC_001806` or `NC_001806.1`,
  naming the translation to run, where it raised a bare `FileNotFoundError` pointing at
  `Crosswalk.available()`.
- `hsv1coords.codon(gene, number)` is new: the reference positions of a codon, in
  reading order, with the strand, for the resistance literature that cites amino acid
  positions.
- `Region.within_gene` and `PositionReport.checks_not_run` are new fields; the first
  was in the artefact and dropped by the loader, the second records a check that could
  not run so that silence is not read as a clean answer.
- `Crosswalk.unmapped_stretch` raises `ValueError` for a position outside the deposit,
  where it returned None, which that method's docstring defines as "that position lifts
  fine".
- `ConventionSet.convert` checks both ends, so a position either convention drops has no
  equivalent under the pair; previously only the target was checked and the two
  directions disagreed.
- `RepeatMap.load` refuses an empty or incomplete map with `RepeatMapIncomplete`, and
  `data/derived/diploid_pairs.tsv` is required rather than optional.
- Every coordinate-taking entry point accepts any integral type that is not a bool,
  `numpy.int64` included, and refuses everything else with one message.

## v0.1.0 - 13 November 2026

First annotation release. Not yet a GFF3: this is the classed and sourced feature table
that a GFF3 will be generated from.

**Added**

- `data/derived/annotation_v0.1.0.tsv`: 1,779 feature rows across three sources (1,723 as first released, before review found the inverted repeats were being dropped), each row
  carrying the class of feature it is, its source record, and that record's provenance.
- `data/derived/annotation_summary.md`: feature counts per class per source.

**Divergences from the source annotations, and why**

- **Nothing is merged across sources.** Each source's features are carried separately.
  Merging requires the coordinate crosswalk (M3), which does not exist yet. A merged
  annotation released before the crosswalk would be asserting equivalences it cannot check.
- **Features without a gene name are kept, from v0.1.0 as corrected.** The first build of
  this release kept only features carrying a `/gene` qualifier, which silently dropped
  56 features without a `/gene` qualifier: 28 repeat regions, among them all four
  inverted repeats (`TRL`, `IRL`, `IRS`, `TRS`) and all three copies of the `'a'`
  sequence, 23 regulatory elements, 3 stem-loops, among them the `oriL` palindrome, and
  2 untyped features, while this changelog claimed nothing was dropped. (Round 15: the
  enumeration here read "both `'a'` sequences ... and 50 other features", which adds to
  57, against "52 other unnamed features" eleven lines below and a true 56; and "both"
  was the two-copy framing PITFALLS root cause 1 retracted. The 56 is derived from the
  record now.) That
  was the worst defect in the release: the inverted repeats are the reason the coordinate
  problem exists. Repeats, stem-loops, origins, regulatory elements and untyped features
  are now carried with their own classes, named from their notes where they have one.
  Only the whole-genome `source` record is excluded, deliberately. 1,779 rows, up from
  1,723.
- Features that no class rule matched are labelled `other` and keep their original notes,
  so they are visible rather than quietly lost.
- **Feature classes are derived from names and qualifiers by stated rules**, in
  `src/build_annotation.py`, not assigned by hand. The rules are in the source and in the
  summary table, so a disputed classification can be traced to the rule that made it.

**The documented discrepancy**

The coding sequence count difference between `NC_001806.2` and `BK012101.1` is
**documented, not resolved**, as the M2 criteria require. `NC_001806.2` carries 77 CDS
features and `BK012101.1` carries 300 over the identical bases. It is not a disagreement:
`BK012101.1` is the 2020 re-annotation deposited as third-party annotation on the same
sequence, and it names 37 named transcript units spanning more than one gene, stated
across 446 rows, counting each gene, CDS and mRNA feature of a unit, and a further 17
nested reading frames (small, upstream, internal and downstream).
The RefSeq annotation does not attempt those classes. A count is only comparable within a
class, which is what the class column exists to enforce. See `FINDINGS.md`.

**Corrections applied to this release after review, 13 November 2026**

- The GenBank parser dropped wrapped location lines, so 8 features shipped with truncated
  coordinates, including `RL2` (ICP0) and the `IRL2-IRL1` span. Fixed, and the M2 gate now
  checks that every location has balanced brackets and parses to coordinates inside the
  genome, which checking only for non-empty did not.
- The 56 features without a `/gene` qualifier, the inverted repeats among them, were
  missing entirely. Fixed above.

## Unreleased - 13 November 2026

Milestones 3 to 6, not yet cut as a version because the acceptance criteria for the
annotation's GFF3 output are still outstanding.

**Added**

- **Strain crosswalk (M3)**: 183 complete genomes aligned to the reference and reduced to
  ungapped coordinate blocks, with `lift` resolving a position between any two of them.
  Median mappable fraction 98.84 per cent.
- **Repeat and diploid-locus map (M4)**: 34 regions covering 21.0 per cent of the genome,
  32 diploid locus pairs derived from the annotation, and `RepeatMap.classify`, which
  answers whether a position is well defined and names the copies it could otherwise be.
- **Convention resolver (M5)**: five published handlings of the repeat regions, each
  quoted from its own paper, with `ConventionSet.report`.
- **Version crosswalk (M6)**: `NC_001806.1` against `NC_001806.2`, nine distinct offsets
  across 22 aligned blocks.

**Corrected during construction**, with the old wording kept in `FINDINGS.md` rather than
quietly replaced: the inverted repeats were missing from the annotation entirely; the
GenBank parser truncated wrapped locations; a lift-over returned a confident coordinate on
the wrong copy of a terminal repeat; LAT was described as lying wholly inside the long
repeat when it extends past it; the `'a'` sequence was modelled with two copies when the
reference carries three; and the Morse convention was encoded in the wrong reference
version, then with boundaries wrong at three of four edges, before being derived from the
paper's stated rule.

## v0.0.0 - 12 November 2026

Inputs and skeleton (M0), gene name dictionary (M1). See `ROADMAP.md`.
