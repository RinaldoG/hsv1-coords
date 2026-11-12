# Plan A: the HSV-1 coordinate and annotation layer

Written 12 November 2026. Build plan, progress only. Acceptance thresholds are written
before the code, per `CLAUDE.md`, and each milestone's gate is pass or fail.

**Repo:** `HSV1_Latency/hsv1-coords/`, its own git repo (conventions: code repos inside a
folder are their own repos). Public from the first commit, for the reason in M0.

**Working name:** `hsv1-coords`.

## What it is, in one paragraph

A machine-readable layer that makes a position in the HSV-1 genome mean one thing. Four
components: a gene name and alias dictionary; one canonical versioned annotation; a
coordinate translator between strains; and an enforceable mask over the regions where a
single linear coordinate is not well defined. Static data plus continuous integration. No
server, no database, no hosting bill, so upkeep is a scheduled job and a few hours when
something upstream moves.

## Design rules, fixed now

1. **Static data, CI-validated.** Everything is a file in git. Every number regenerates
   from one command in a container. Nothing requires a running service.
2. **Provenance per row.** Every gene, alias, feature and alignment carries where it came
   from. A row without a source does not ship.
3. **The mask is enforceable, not advisory.** The translator refuses to return a confident
   coordinate inside a masked region. Refusing is a feature, and a worked example
   demonstrates it refusing.
4. **Build on the 2020 re-annotation, cited, not re-derived.** Whisnant et al. 2020
   (Zenodo `10.5281/zenodo.3724023`, 113 kB, 170 genes / 201 transcripts / 309 CDS, parsed
   and confirmed in strand `03`) is the best existing annotation. It is the input, credited
   in the README, the paper and every derived file's provenance. Re-deriving it would be
   worse and slower.
5. **Semantic versioning with a changelog.** The thing being fixed is partly that the
   reference has no versioning discipline. This one does from v0.1.0.

## Milestones

### M0: skeleton and inputs (one evening)

Public repo, container, CI, licence (MIT), and the inputs pulled and checksummed:

| Input | Where | Size |
|---|---|---|
| Whisnant 2020 annotation | Zenodo `10.5281/zenodo.3724023`, note its seqname is `JN555585`, not `NC_001806.2` | 113 kB |
| Reference record | `NC_001806.2`, 77 CDS, provisional, last modified 13 Aug 2018, assembler `Maq v. 0.7.1` | small |
| The 300-CDS record on identical bases | `BK012101.1`, same 152,222 bp | small |
| `JN555585` | the Whisnant coordinate base | small |
| Complete genomes | 185 of the 597 GenBank records labelled complete | ~92 MB for all 597 |
| Isomer-resolved genomes | `ON783214` to `ON783217` from `PRJNA846339` | small |

**Why public tonight, and it is a scheduling fact not a preference:** JOSS treats commits
concentrated in the weeks before submission, or a repo made public just before submission,
as not acceptable, and it requires roughly six months of history plus documented external
use for a single author (strand `04`). The commit clock starts when the repo is public. A
repo started tonight is submittable in May; one started in April is not.

**Gate:** `docker run ... make inputs` pulls and verifies every input from a clean machine.

### M1: the gene name dictionary (1 to 2 weeks)

One row per gene: canonical identifier, systematic name (`UL48`), protein name (`VP16`),
every historical alias (`alpha-TIF`, `Vmw65`), and a citation per alias. Sources: the
Whisnant annotation, the GenBank feature tables of all three records above, and the
synonym usage in the literature.

Note from strand `03`: the canonical protein count is contested. "About 80 HSV-1 proteins"
is the 1988 list; Whisnant resolved 284 open reading frames, and the 46 novel large ones
have no interaction data and appeared in no screen. The dictionary states which inventory
each row belongs to.

**Gate:**
- Name to canonical to name is lossless for every entry.
- Reproduces the synonymy used in at least three published papers, each named in the test
  suite.
- Zero rows with an uncited alias.

### M2: the canonical annotation (2 to 3 weeks)

One versioned GFF3, released with a changelog, reconciling the 77-CDS and 300-CDS views of
the identical 152,222 bases, and stating which coordinate base is canonical and why.

**Gate:**
- Parses with two independent GFF3 parsers, zero warnings.
- Reproduces the Whisnant 2020 feature set exactly where it claims to, with every
  divergence listed and justified in the changelog.
- Every feature carries a provenance record.
- The 77-against-300 discrepancy is documented, not silently resolved. Silent resolution is
  a fail even if the answer is right.

### M3: the strain crosswalk and lift-over (2 to 3 weeks)

Pairwise whole-genome alignments across the complete genomes, chain files, and a lift-over
command line tool. This is where the three "strain KOS" records at 152,011, 151,024 and
151,974 bp, none linked to raw reads, stop being a silent error across roughly 987 bp.

**Gate:**
- Round-trip A to B to A returns the original coordinate for at least 99 per cent of
  positions outside masked regions, with the failures enumerated individually, not averaged.
- Resolves a published coordinate from at least five papers across at least three strains,
  each checked by hand against the paper.
- Returns no confident coordinate anywhere inside a masked region.

### M4: the repeat and un-callable mask (1 to 2 weeks)

The regions where one linear coordinate is not well defined: `TRL` 1..9,213 and `IRL`
117,161..126,373 (9,213 bp each, identical), `IRS` 125,975..132,607 and `TRS`
145,590..152,222 (6,633 bp each), the GC-extreme stretches (LAT sits at 71.6 per cent GC
against 66.4 per cent in mid-UL), and the segment boundaries that differ across the four
isomers (roughly 27/21/23/29 per cent).

This is the component that makes every downstream conservation and guide-design claim
honest, because LAT lies wholly inside `TRL`/`IRL` in two indistinguishable copies, with
one annotated copy running off the start of the reference as `complement(<1..7569)`, and
CTRL2 inside the repeat as well.

**Gate:**
- Reproduces the four repeat boundaries from the reference record's own feature table.
- A worked example shows a downstream tool refusing to make a claim inside the mask.
- The mask is machine-readable and enforced in code, not documented in prose.

### M5: release (1 week)

Zenodo DOI, container image, README that a stranger can follow to reproduce a release, and
the v1.0.0 tag. Then offered to `DepledgeLab/HerpesvirusGenomicsWiki` as a pull request
filling its empty "Software", "Databases" and "Workflows" headings, and announced where the
field reads: a preprint, and an abstract to the International Herpesvirus Workshop, which
strand `04` found accepts abstracts on the work rather than the credential.

**Gate:** every figure, table and number in the repo and the preprint regenerates from one
command. Tested on someone who has not seen it.

### M6: the pangenome graph (separate decision, not part of v1)

Becomes possible once M3 and M4 exist, and not before. `vg` or `minigraph-cactus` over the
complete genomes plus the isomer-resolved set. Searches on 12 November 2026 found no `vg`,
`minigraph` or `minigraph-cactus` work on any herpesvirus; the method is proven on other
large DNA viruses (a 2025 bioRxiv preprint on lumpy skin disease virus, DOI
10.1101/2025.11.26.690900, and an HBV variation graph, both search-summary-only in strand
`03`). Decided in `DECISIONS.md` after M5, not assumed here.

## Timeline

Seven to eleven weeks to v1.0.0 at a part-time pace, M1 to M5. M0 is tonight.

## What tells you it worked

In order of how much it would mean:

1. Someone else's paper or pipeline resolves coordinates with it.
2. A gene-editing programme uses the mask or the crosswalk in target design. There is a
   live buyer: `NCT06474442`, BD111, Shanghai BDgene, Phase IIa, recruiting since 28 April
   2025.
3. The wiki pull request is merged.
4. A citation in a methods section.

Not stars, not downloads.

## Known risks, and what is done about each

| Risk | What is done |
|---|---|
| The 77-against-300 CDS reconciliation turns out to be a genuine scientific disagreement rather than a bookkeeping one | Document both, pick one as canonical, make the other derivable. The layer does not need to settle it, only to stop it being invisible |
| Published coordinates cannot be resolved because papers do not say which record they used | That is a finding, logged. M3's hand-check set is the evidence, and it becomes part of the write-up |
| Dead upstream tooling cannot be revived | Not required. Assembly revival is deliberately out of v1: the layer works on deposited sequence |
| The work goes unused | M5 pushes it through three channels, not one. The adoption test is explicit above and failing it is recorded, not hidden |
