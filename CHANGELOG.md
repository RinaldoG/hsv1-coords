# Changelog

Semantic versioning. The annotation is versioned because the thing being fixed is partly
that the reference has no versioning discipline.

## v0.1.0 - 13 November 2026

First annotation release. Not yet a GFF3: this is the classed and sourced feature table
that a GFF3 will be generated from.

**Added**

- `data/derived/annotation_v0.1.0.tsv`: 1,723 feature rows across three sources, each row
  carrying the class of feature it is, its source record, and that record's provenance.
- `data/derived/annotation_summary.md`: feature counts per class per source.

**Divergences from the source annotations, and why**

- **Nothing is merged across sources.** Each source's features are carried separately.
  Merging requires the coordinate crosswalk (M3), which does not exist yet. A merged
  annotation released before the crosswalk would be asserting equivalences it cannot check.
- **No feature is dropped**, including the 10 to 14 per source that no class rule matched.
  Those are labelled `other` and keep their original notes, so they are visible rather
  than quietly lost.
- **Feature classes are derived from names and qualifiers by stated rules**, in
  `src/build_annotation.py`, not assigned by hand. The rules are in the source and in the
  summary table, so a disputed classification can be traced to the rule that made it.

**The documented discrepancy**

The 77-against-300 coding sequence difference between `NC_001806.2` and `BK012101.1` is
**documented, not resolved**, as the M2 criteria require. It is not a disagreement:
`BK012101.1` is the 2020 re-annotation deposited as third-party annotation on the same
sequence, and 452 of its named features are transcript units spanning more than one gene,
with a further 37 being nested reading frames (small, upstream, internal and downstream).
The RefSeq annotation does not attempt those classes. A count is only comparable within a
class, which is what the class column exists to enforce. See `FINDINGS.md`.

## v0.0.0 - 12 November 2026

Inputs and skeleton (M0), gene name dictionary (M1). See `ROADMAP.md`.
