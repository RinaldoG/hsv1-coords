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

## M1: gene name dictionary - in progress

One row per gene: canonical identifier, systematic name (`UL48`), protein name (`VP16`),
and every historical alias (`alpha-TIF`, `Vmw65`), each alias carrying a citation.

The mechanical half is done: `src/extract_names.py` extracts every gene name from the
deposited records and reconciles them across sources, without inventing synonymy. Results
in `data/derived/`. The finding that shapes this milestone is in `FINDINGS.md`: of 232
distinct gene names across the four sources, five appear in all four.

**Criteria:**
- Name to canonical to name is lossless for every entry.
- Reproduces the synonymy used in at least three published papers, each named in the test
  suite.
- No entry carries an uncited alias. Names believed to exist but lacking a citation are
  listed separately as unverified.
- Where a name has been used for more than one gene, the ambiguity is recorded rather than
  resolved.

## M2: canonical versioned annotation

One semantically versioned annotation with a changelog and per-feature provenance, stating
which annotation class each feature belongs to, so that counts from differently scoped
annotations are not compared as though they measured the same thing.

**Criteria:**
- Parses with two independent GFF3 parsers, zero warnings.
- Reproduces its source annotations exactly where it claims to, with every divergence listed
  and justified in the changelog.
- Every feature carries a provenance record naming its source.

## M3: strain coordinate crosswalk and lift-over

Pairwise alignments across the complete genomes, chain files, and a lift-over tool, so that
a coordinate cited against one strain resolves against another. Three deposited genomes of
one commonly used laboratory strain differ in length by up to 987 bp with no published
mapping between them.

**Criteria:**
- Round-trip lift-over returns the original coordinate for at least 99 per cent of positions
  outside masked regions, with failures enumerated individually rather than averaged.
- Resolves a published coordinate from at least five papers across at least three strains,
  each checked by hand against the paper.
- Returns no confident coordinate inside a masked region.

## M4: repeat and un-callable region mask

The regions where a single linear coordinate is not well defined: the inverted repeats
(`TRL` and `IRL` at 9,213 bp each, `IRS` and `TRS` at 6,633 bp each), the nested reiteration
sets inside them, the GC-extreme stretches, and the segment boundaries that differ between
the genome's four isomers.

This matters most where it is least convenient. LAT lies wholly inside the inverted repeats,
in two copies, and in the reference record one of those copies is truncated by the start of
the sequence.

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
