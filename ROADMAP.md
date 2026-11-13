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

## M4: repeat and un-callable region mask

The regions where a single linear coordinate is not well defined: the inverted repeats
(`TRL` and `IRL` at 9,213 bp each, `IRS` and `TRS` at 6,633 bp each), the nested reiteration
sets inside them, the GC-extreme stretches, and the segment boundaries that differ between
the genome's four isomers.

This matters most where it is least convenient. LAT sits in the long inverted repeat and
extends 778 bp past it into `IRS`, exists in two copies, and in the reference record one of
those copies is truncated because the genome is circular and its remainder lies in `TRS` at
the other end of the linear sequence. So the mask has to express three different things
about one gene: duplicated, spanning a segment boundary, and split across the two ends of
the record.

**Criteria:**
- Reproduces the repeat boundaries from the reference record's own feature table.
- Machine-readable and enforced in code, not described in prose.
- A worked example demonstrates a downstream tool declining to make a claim inside the mask.

## M5: release

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
