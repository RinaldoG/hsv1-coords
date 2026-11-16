#!/usr/bin/env python3
"""Render the M2 annotation as GFF3, and the repeat map as BED (roadmap M2 and M4).

The annotation table carries GenBank location syntax (`complement(join(...))`, with
fuzzy `<` and `>` ends), which no genome browser reads. This renders the same rows in
the formats tools do read. It adds no feature and drops none: every row of
`annotation_v0.1.0.tsv` appears in its source's GFF3 exactly once, which
`tests/test_m2_gff3.py` checks tuple by tuple.

**One GFF3 per source, not one merged file.** The three sources describe the identical
152,222 bp under three sequence names (`NC_001806.2`, `BK012101.1`, and `JN555585` in
the GTF) and with different scope. One file per source keeps each a faithful rendering
of one record under that record's own sequence name, and does not present the same
locus stated by two annotations as one annotation.

How the GenBank constructs are represented, and why:

- **Strand** comes from `complement(...)`.
- **A spliced CDS** (`join(...)`) is one feature over several lines sharing one ID, each
  line carrying its own phase. That is the GFF3 1.26 convention for a CDS.
- **A spliced mRNA or ncRNA** is one line spanning the whole feature plus one `exon` line
  per segment, each with `Parent` set to it. The join itself states those segments, so
  the exons are a rendering of the source, not an invented relationship; a transcript
  spread over several lines sharing an ID is legal GFF3 but is drawn badly or refused by
  most browsers. An unspliced mRNA gets no exon line. No other feature type is spliced in
  these records, and one that is fails the build rather than being guessed at.
- **Phase** is computed from the segments and `/codon_start` (the reference states 1 on
  all of its CDS today, and it is read rather than assumed). A CDS with a partial 5' end
  takes its phase from `/codon_start` like any other. The GTF states its own phase, which
  is written as stated after being checked against the same arithmetic. The suite checks
  every phase against the sequence by translation, not against this arithmetic.
- **Fuzzy ends** follow the NCBI GFF3 convention: `partial=true` plus
  `start_range=.,N` for a `<` end and `end_range=N,.` for a `>` end. Like columns 4 and
  5, these name the left and right ends on the sequence, whatever the strand.
- **Parent** is written only where a source states one. The GenBank records state none
  (a shared `/gene` value is a name, not a stated relationship), so their genes, mRNAs
  and CDS are flat features with IDs. The GTF states `gene_id` and `transcript_id`, so
  its transcripts point at genes, exons and CDS at transcripts, and a CDS with no
  `transcript_id` (47 lines do this) at its gene.

The repeat map goes to BED6+2 (`repeat_map.bed`) and the diploid pairs to BEDPE
(`diploid_pairs.bedpe`). BED is 0-based and half-open where every other file here is
1-based and inclusive, so a start is written as `start - 1` and an end as is; the suite
proves that frame by reading bases with the BED coordinates.

Run: python3 -I src/build_gff3.py
"""
from __future__ import annotations

import collections
import csv
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"

sys.path.insert(0, str(ROOT / "src"))
from build_annotation import SOURCES, VERSION, annotation_rows  # noqa: E402
from build_crosswalk import genbank_sequence  # noqa: E402

GFF3_FILES = {source: f"annotation_v{VERSION}_{source}.gff3"
              for source, _, _, _ in SOURCES}
REPEAT_BED = "repeat_map.bed"
PAIRS_BEDPE = "diploid_pairs.bedpe"
OUTPUTS = tuple(sorted(GFF3_FILES.values())) + (REPEAT_BED, PAIRS_BEDPE)
"""Every file this build writes, for the gates and the checkout check to read."""

SEQUENCE_RECORD = {
    # The record each source's coordinates are on, read for its length. The GTF's
    # sequence name is `JN555585`, and that record's sequence is byte-identical to
    # NC_001806.2's (`make check` asserts it).
    "NC_001806.2": "ref_NC_001806.2.gb",
    "BK012101.1": "alt_BK012101.1.gb",
    "whisnant2020-gtf": "base_JN555585.1.gb",
}

SO_TYPES: dict[tuple[str, str], str] = {
    # (source format, source feature type) -> Sequence Ontology term. Declared one by
    # one. An unmapped type stops the build rather than defaulting to a generic term,
    # because a default would publish a guess with nothing to say so.
    ("genbank", "gene"): "gene",                          # SO:0000704
    ("genbank", "mRNA"): "mRNA",                          # SO:0000234
    ("genbank", "CDS"): "CDS",                            # SO:0000316
    ("genbank", "ncRNA"): "ncRNA",                        # SO:0000655
    ("genbank", "repeat_region"): "repeat_region",        # SO:0000657
    ("genbank", "stem_loop"): "stem_loop",                # SO:0000313
    ("genbank", "rep_origin"): "origin_of_replication",   # SO:0000296
    # `misc_feature` says only that something is there. NCBI's own GFF3 renders it as
    # the root term, and so does this.
    ("genbank", "misc_feature"): "sequence_feature",      # SO:0000110
    ("gtf", "gene"): "gene",                              # SO:0000704
    # Not `mRNA`: the GTF's `gene_type` is "protein_coding" on every line, LAT and the
    # pri-miRNAs included, so it cannot be read as saying which transcripts code.
    ("gtf", "transcript"): "transcript",                  # SO:0000673
    ("gtf", "exon"): "exon",                              # SO:0000147
    ("gtf", "CDS"): "CDS",                                # SO:0000316
}

REGULATORY_CLASSES: dict[str, str] = {
    # A GenBank `regulatory` feature is typed by its `/regulatory_class`, which SO has
    # a term for in each case present. An unlisted class stops the build.
    "TATA_box": "TATA_box",                               # SO:0000174
    "polyA_signal_sequence": "polyA_signal_sequence",     # SO:0000551
}

SO_ACCESSIONS = {
    "gene": "SO:0000704", "mRNA": "SO:0000234", "CDS": "SO:0000316",
    "ncRNA": "SO:0000655", "repeat_region": "SO:0000657", "stem_loop": "SO:0000313",
    "origin_of_replication": "SO:0000296", "sequence_feature": "SO:0000110",
    "transcript": "SO:0000673", "exon": "SO:0000147", "TATA_box": "SO:0000174",
    "polyA_signal_sequence": "SO:0000551",
}

SPLICED_AS_EXONS = {"mRNA", "ncRNA", "transcript"}
"""Source types whose joins are written as a spanning line plus `exon` children."""

BED_STRAND = {"forward": "+", "reverse": "-", "single": ".", "unequal_length": "."}
"""`repeat_map.tsv`'s orientation of each copy, as a BED strand.

`forward` and `reverse` say which way a copy of a multi-copy repeat reads relative to
its first copy, which is what a BED strand expresses. `single` (a region with one copy
per orientation, such as TRL) and `unequal_length` (copies whose orientation was not
measured because their lengths differ) state no strand, so they get `.`."""


class UnsupportedLocation(ValueError):
    """A location this build does not know how to render, named rather than guessed."""


SEGMENT = re.compile(r"^(<?)(\d+)\.\.(>?)(\d+)$")


def parse_location(location: str) -> tuple[str, list[tuple[int, int]], bool, bool]:
    """(strand, segments in sequence order, fuzzy left end, fuzzy right end).

    Understands `a..b`, `<a..b`, `a..>b`, `join(...)` and `complement(...)` around
    either. Anything else (`order`, a single base, a `join` of complements, a fuzzy
    end inside a join) raises, because each would need a decision this build has not
    made.
    """
    text = location.strip()
    strand = "+"
    if text.startswith("complement(") and text.endswith(")"):
        strand, text = "-", text[len("complement("):-1]
    if text.startswith("join(") and text.endswith(")"):
        parts = text[len("join("):-1].split(",")
    else:
        parts = [text]
    segments: list[tuple[int, int]] = []
    fuzzy: list[tuple[bool, bool]] = []
    for part in parts:
        m = SEGMENT.match(part.strip())
        if not m:
            raise UnsupportedLocation(f"cannot render location {location!r}")
        start, end = int(m.group(2)), int(m.group(4))
        if start > end:
            raise UnsupportedLocation(f"segment runs backwards in {location!r}")
        segments.append((start, end))
        fuzzy.append((bool(m.group(1)), bool(m.group(3))))
    for (_, end), (start, _) in zip(segments, segments[1:]):
        if start <= end:
            raise UnsupportedLocation(f"segments overlap or are out of order in "
                                      f"{location!r}")
    if any(f[0] for f in fuzzy[1:]) or any(f[1] for f in fuzzy[:-1]):
        raise UnsupportedLocation(f"a fuzzy end inside a join in {location!r}")
    return strand, segments, fuzzy[0][0], fuzzy[-1][1]


def phases(segments: list[tuple[int, int]], strand: str, codon_start: int) -> list[int]:
    """GFF3 phase of each segment, in the order `segments` is given (sequence order).

    Phase is the number of bases to skip at a segment's 5' end, on its own strand,
    before the first complete codon. The first segment in reading order skips
    `codon_start - 1`; each later one skips whatever completes the codon the previous
    segments left open.
    """
    if codon_start not in (1, 2, 3):
        raise ValueError(f"codon_start {codon_start} is not 1, 2 or 3")
    reading = segments if strand == "+" else list(reversed(segments))
    out: dict[tuple[int, int], int] = {}
    before = 0
    for segment in reading:
        out[segment] = (codon_start - 1 - before) % 3
        before += segment[1] - segment[0] + 1
    return [out[s] for s in segments]


RESERVED = re.compile(r"[\x00-\x1f\x7f%;=&,]")


def escape(value: str) -> str:
    """Percent-encode what GFF3 1.26 reserves in column 9.

    Tab, newline, carriage return, control characters and `%` must be encoded
    anywhere; `;`, `=`, `&` and `,` have meanings in column 9 and must be encoded
    inside a value. A comma left in a `Note` would split it into several values.
    """
    return RESERVED.sub(lambda m: f"%{ord(m.group(0)):02X}", value)


def unescape(value: str) -> str:
    return re.sub(r"%([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), value)


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", name).strip("_") or "unnamed"


def attributes(pairs: list[tuple[str, str]]) -> str:
    return ";".join(f"{key}={escape(value)}" for key, value in pairs if value != "")


class Feature:
    """One GFF3 feature: one line, or several sharing an ID, plus its children."""

    def __init__(self, seqid: str, source: str, so_type: str, strand: str,
                 segments: list[tuple[int, int]], attrs: list[tuple[str, str]],
                 phase: list[int] | None = None):
        self.seqid, self.source, self.so_type = seqid, source, so_type
        self.strand, self.segments, self.attrs = strand, segments, attrs
        self.phase = phase
        self.children: list[Feature] = []

    @property
    def id(self) -> str:
        return dict(self.attrs)["ID"]

    @property
    def sort_key(self) -> tuple:
        rank = ["gene", "transcript", "mRNA", "ncRNA", "exon", "CDS"]
        position = rank.index(self.so_type) if self.so_type in rank else len(rank)
        return (self.segments[0][0], -self.segments[-1][1], position, self.id)

    def lines(self) -> list[str]:
        out = []
        for i, (start, end) in enumerate(self.segments):
            phase = str(self.phase[i]) if self.phase is not None else "."
            out.append("\t".join([self.seqid, self.source, self.so_type, str(start),
                                  str(end), ".", self.strand, phase,
                                  attributes(self.attrs)]))
        for child in sorted(self.children, key=lambda c: c.sort_key):
            out.extend(child.lines())
        return out


def fuzzy_attributes(left: bool, right: bool, start: int, end: int) -> list[tuple]:
    out = []
    if left or right:
        out.append(("partial", "true"))
    if left:
        out.append(("start_range", f".,{start}"))
    if right:
        out.append(("end_range", f"{end},."))
    return out


def common_attributes(row: dict[str, str]) -> list[tuple[str, str]]:
    return [("Name", row["name"]), ("Note", row["notes"]),
            ("feature_class", row["feature_class"]), ("product", row["product"]),
            ("subtype", row["subtype"]), ("transcript_id", row["transcript_id"])]


def so_type_of(kind: str, row: dict[str, str]) -> str:
    if kind == "genbank" and row["feature_type"] == "regulatory":
        if row["subtype"] not in REGULATORY_CLASSES:
            raise KeyError(f"regulatory class {row['subtype']!r} of {row['name']} has "
                           f"no Sequence Ontology mapping; add one to "
                           f"REGULATORY_CLASSES")
        return REGULATORY_CLASSES[row["subtype"]]
    key = (kind, row["feature_type"])
    if key not in SO_TYPES:
        raise KeyError(f"{kind} feature type {row['feature_type']!r} of {row['name']} "
                       f"has no Sequence Ontology mapping; add one to SO_TYPES")
    return SO_TYPES[key]


def unique_id(wanted: str, taken: set[str]) -> str:
    candidate, n = wanted, 1
    while candidate in taken:
        n += 1
        candidate = f"{wanted}-{n}"
    taken.add(candidate)
    return candidate


def genbank_features(source: str, seqid: str,
                     pairs: list[tuple[dict, dict]]) -> list[Feature]:
    """Flat features, in record order, from one GenBank record's rows."""
    taken: set[str] = set()
    features = []
    for row, raw in pairs:
        so_type = so_type_of("genbank", row)
        strand, segments, left, right = parse_location(row["location"])
        start, end = segments[0][0], segments[-1][1]
        protein_id = (raw.get("protein_id") or "").strip()
        if row["feature_type"] == "CDS" and protein_id:
            wanted = f"cds-{protein_id}"
        else:
            wanted = f"{so_type}-{slug(row['name'])}-{start}-{end}"
        if wanted in taken and row["feature_type"] == "CDS" and protein_id:
            raise ValueError(f"protein_id {protein_id} is stated twice in {source}")
        fid = unique_id(wanted, taken)
        attrs = ([("ID", fid)] + common_attributes(row)
                 + [("protein_id", protein_id)]
                 + fuzzy_attributes(left, right, start, end)
                 + [("provenance", row["provenance"])])
        if row["feature_type"] == "CDS":
            codon_start = int((raw.get("codon_start") or "1").strip())
            features.append(Feature(seqid, source, so_type, strand, segments, attrs,
                                    phases(segments, strand, codon_start)))
            continue
        if len(segments) > 1 and row["feature_type"] not in SPLICED_AS_EXONS:
            raise UnsupportedLocation(f"{row['feature_type']} {row['name']} is spliced "
                                      f"({row['location']}), and only a CDS or a "
                                      f"transcript has a GFF3 form decided here")
        parent = Feature(seqid, source, so_type, strand, [(start, end)], attrs)
        if len(segments) > 1:
            reading = segments if strand == "+" else list(reversed(segments))
            for number, (a, b) in enumerate(reading, 1):
                exon_attrs = ([("ID", unique_id(f"{fid}-exon{number}", taken)),
                               ("Parent", fid), ("Name", row["name"])]
                              + fuzzy_attributes(left and a == start,
                                                 right and b == end, a, b)
                              + [("provenance", row["provenance"])])
                parent.children.append(
                    Feature(seqid, source, "exon", strand, [(a, b)], exon_attrs))
        features.append(parent)
    return features


def gtf_features(source: str, seqid: str,
                 pairs: list[tuple[dict, dict]]) -> list[Feature]:
    """The GTF's own hierarchy: gene, transcript, then exon and CDS."""
    genes: dict[str, Feature] = {}
    transcripts: dict[str, Feature] = {}
    # Keyed by (orf_id, gene_id, transcript_id, strand), not orf_id alone: the GTF
    # reuses `ORF-O_CDS` and `ORF-P_CDS` for the TRL copy (gene `ORF-O`, minus strand)
    # and the IRL copy (gene `ORF-O/P`, plus strand), which are different features.
    cds_lines: dict[tuple, list[tuple[dict, dict]]] = collections.defaultdict(list)
    exons: list[tuple[dict, dict]] = []
    for row, raw in pairs:
        if raw["_seqid"] != seqid:
            raise ValueError(f"GTF line on {raw['_seqid']}, expected {seqid}")
        kind = row["feature_type"]
        so_type = so_type_of("gtf", row)
        strand, segments, left, right = parse_location(row["location"])
        if len(segments) != 1 or left or right:
            raise UnsupportedLocation(f"GTF line with location {row['location']!r}")
        attrs = common_attributes(row) + [("gene_id", raw["gene"])]
        if kind == "gene":
            fid = f"gene-{raw['gene']}"
            if raw["gene"] in genes:
                raise ValueError(f"gene_id {raw['gene']} has two gene lines")
            genes[raw["gene"]] = Feature(
                seqid, source, so_type, strand, segments,
                [("ID", fid)] + attrs + [("provenance", row["provenance"])])
        elif kind == "transcript":
            tid = raw["transcript"]
            if not tid or tid in transcripts:
                raise ValueError(f"transcript line with id {tid!r} missing or repeated")
            transcripts[tid] = Feature(
                seqid, source, so_type, strand, segments,
                [("ID", f"rna-{tid}"), ("Parent", f"gene-{raw['gene']}")] + attrs
                + [("provenance", row["provenance"])])
            transcripts[tid].gene_id = raw["gene"]
        elif kind == "exon":
            exons.append((row, raw))
        elif kind == "CDS":
            if not raw["orf_id"]:
                raise ValueError(f"CDS line of {raw['gene']} has no orf_id")
            cds_lines[(raw["orf_id"], raw["gene"], raw["transcript"],
                       raw["strand"])].append((row, raw))
        else:                                        # pragma: no cover
            raise KeyError(kind)

    for tid, transcript in transcripts.items():
        if transcript.gene_id not in genes:
            raise ValueError(f"transcript {tid} names gene_id {transcript.gene_id}, "
                             f"which has no gene line")
        genes[transcript.gene_id].children.append(transcript)

    by_transcript: dict[str, list[tuple[dict, dict]]] = collections.defaultdict(list)
    for row, raw in exons:
        if raw["transcript"] not in transcripts:
            raise ValueError(f"exon of {raw['gene']} names transcript "
                             f"{raw['transcript']!r}, which has no transcript line")
        by_transcript[raw["transcript"]].append((row, raw))
    for tid, members in sorted(by_transcript.items()):
        transcript = transcripts[tid]
        ordered = sorted(members, key=lambda p: parse_location(p[0]["location"])[1])
        if transcript.strand == "-":
            ordered.reverse()
        for number, (row, raw) in enumerate(ordered, 1):
            strand, segments, _, _ = parse_location(row["location"])
            transcript.children.append(Feature(
                seqid, source, "exon", strand, segments,
                [("ID", f"exon-{tid}-{number}"), ("Parent", f"rna-{tid}")]
                + common_attributes(row) + [("gene_id", raw["gene"]),
                                            ("provenance", row["provenance"])]))

    taken: set[str] = set()
    groups = sorted(cds_lines.items(),
                    key=lambda kv: (kv[0][0], min(parse_location(r["location"])[1][0]
                                                  for r, _ in kv[1])))
    for (orf_id, gene_id, tid, _), members in groups:
        located = sorted(((parse_location(r["location"]), r, raw)
                          for r, raw in members), key=lambda x: x[0][1])
        strand = located[0][0][0]
        segments = [loc[1][0] for loc, _, _ in located]
        if any(b >= a2 for (_, b), (a2, _) in zip(segments, segments[1:])):
            raise ValueError(f"orf_id {orf_id} of {gene_id}: segments overlap")
        stated = [int(raw["_phase"]) for _, _, raw in located]
        reading_first = located[0] if strand == "+" else located[-1]
        computed = phases(segments, strand, int(reading_first[2]["_phase"]) + 1)
        if stated != computed:
            raise ValueError(f"orf_id {orf_id}: the GTF states phases {stated} and the "
                             f"segments give {computed}")
        # One ID over several lines, and the attributes of a multi-line feature must
        # agree, so a line-specific value (the GTF's `comment` is on some lines of a
        # spliced ORF and not others) is not carried as an attribute of the whole.
        notes = {r["notes"] for _, r, _ in located}
        products = {r["product"] for _, r, _ in located}
        row = dict(located[0][1])
        if len(notes) != 1 or len(products) != 1:
            raise ValueError(f"orf_id {orf_id}: lines disagree on notes or product")
        parent = f"rna-{tid}" if tid else f"gene-{gene_id}"
        feature = Feature(seqid, source, "CDS", strand, segments,
                          [("ID", unique_id(f"cds-{orf_id}", taken)),
                           ("Parent", parent)]
                          + common_attributes(row)
                          + [("gene_id", gene_id), ("orf_id", orf_id),
                             ("provenance", row["provenance"])],
                          stated)
        if tid:
            if tid not in transcripts:
                raise ValueError(f"CDS {orf_id} names transcript {tid}, which has no "
                                 f"transcript line")
            transcripts[tid].children.append(feature)
        else:
            if gene_id not in genes:
                raise ValueError(f"CDS {orf_id} names gene {gene_id}, which has no "
                                 f"gene line")
            genes[gene_id].children.append(feature)
    return list(genes.values())


def sequence_length(source: str, raw: pathlib.Path = RAW) -> int:
    return len(genbank_sequence(raw / SEQUENCE_RECORD[source]))


def render_gff3(raw: pathlib.Path = RAW) -> dict[str, str]:
    by_source: dict[str, list[tuple[dict, dict]]] = collections.defaultdict(list)
    for row, feature in annotation_rows(raw):
        by_source[row["source"]].append((row, feature))
    out = {}
    for source, _, kind, provenance in SOURCES:
        pairs = by_source[source]
        if kind == "genbank":
            seqid = source
            features = genbank_features(source, seqid, pairs)
        else:
            seqids = {feature["_seqid"] for _, feature in pairs}
            if len(seqids) != 1:
                raise ValueError(f"{source} names {len(seqids)} sequences: {seqids}")
            seqid = seqids.pop()
            features = gtf_features(source, seqid, pairs)
        length = sequence_length(source, raw)
        ids = collections.Counter()
        lines = [
            "##gff-version 3",
            f"##sequence-region {seqid} 1 {length}",
            f"# hsv1-coords annotation v{VERSION}, source {source}: {provenance}.",
            f"# Rendered by src/build_gff3.py from data/derived/annotation_v{VERSION}"
            f".tsv; one line per row of that table for this source, except that a "
            f"spliced CDS takes one line per segment under one ID and a spliced mRNA "
            f"or ncRNA takes one exon line per segment.",
            "# Coordinates are 1-based and inclusive. Fuzzy ends: partial=true, with "
            "start_range=.,N or end_range=N,. (NCBI convention).",
        ]
        # An ID belongs to one feature. Several lines may share it only because they
        # are the segments of that one feature, so the check is per Feature object.
        def walk(feature):
            yield feature
            for child in feature.children:
                yield from walk(child)
        for feature in features:
            for one in walk(feature):
                ids[one.id] += 1
        repeated = sorted(i for i, n in ids.items() if n > 1)
        if repeated:
            raise ValueError(f"{source}: IDs given to more than one feature: "
                             f"{repeated[:5]}")
        for feature in sorted(features, key=lambda f: f.sort_key):
            lines.extend(feature.lines())
        out[GFF3_FILES[source]] = "\n".join(lines) + "\n"
    return out


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def render_bed(derived: pathlib.Path = DERIVED) -> dict[str, str]:
    """The repeat map as BED6+2 and the diploid pairs as BEDPE, on NC_001806.2."""
    chrom = "NC_001806.2"
    regions = read_tsv(derived / "repeat_map.tsv")
    bed = ["# hsv1-coords repeat map on NC_001806.2, rendered by src/build_gff3.py "
           "from data/derived/repeat_map.tsv.",
           "# BED6+2: 0-based half-open (chromStart = the map's 1-based start - 1, "
           "chromEnd = its end). Extra columns: kind, orientation.",
           "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tkind\torientation"]
    rows = []
    for r in regions:
        if r["orientation"] not in BED_STRAND:
            raise KeyError(f"orientation {r['orientation']!r} of {r['name']} has no "
                           f"BED strand; add it to BED_STRAND")
        start, end = int(r["start"]), int(r["end"])
        if not 1 <= start <= end:
            raise ValueError(f"{r['name']}: {start}..{end} is not a 1-based interval")
        rows.append((start - 1, end, r["name"], r["kind"], r["orientation"]))
    for start0, end, name, kind, orientation in sorted(rows):
        bed.append("\t".join([chrom, str(start0), str(end), name, "0",
                              BED_STRAND[orientation], kind, orientation]))

    pairs = read_tsv(derived / "diploid_pairs.tsv")
    bedpe = ["# hsv1-coords diploid pairs on NC_001806.2, rendered by "
             "src/build_gff3.py from data/derived/diploid_pairs.tsv.",
             "# BEDPE: each copy 0-based half-open. strand1/strand2 are +/- where the "
             "two copies were measured as inverted (each the reverse complement of "
             "the other) and ./. otherwise. Extra columns: repeat_pair, feature_type, "
             "strand_relationship, offset_delta, copy_1_name, copy_2_name.",
             "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\t"
             "strand2\trepeat_pair\tfeature_type\tstrand_relationship\toffset_delta\t"
             "copy_1_name\tcopy_2_name"]
    lines = []
    for p in pairs:
        a1, b1 = int(p["copy_1_start"]), int(p["copy_1_end"])
        a2, b2 = int(p["copy_2_start"]), int(p["copy_2_end"])
        inverted = p["strand_relationship"] == "inverted"
        lines.append(((a1 - 1, b1, a2 - 1, b2, p["locus"], p["feature_type"]), [
            chrom, str(a1 - 1), str(b1), chrom, str(a2 - 1), str(b2), p["locus"], "0",
            "+" if inverted else ".", "-" if inverted else ".", p["repeat_pair"],
            p["feature_type"], p["strand_relationship"], p["offset_delta"],
            p["copy_1_name"], p["copy_2_name"]]))
    for _, fields in sorted(lines, key=lambda x: x[0]):
        bedpe.append("\t".join(fields))
    return {REPEAT_BED: "\n".join(bed) + "\n", PAIRS_BEDPE: "\n".join(bedpe) + "\n"}


def render(raw: pathlib.Path = RAW, derived: pathlib.Path = DERIVED) -> dict[str, str]:
    """Every output, as text, keyed by file name. The suite compares it with disk."""
    out = render_gff3(raw)
    out.update(render_bed(derived))
    return out


def main() -> int:
    try:
        outputs = render()
    except FileNotFoundError as missing:
        print(f"{missing}; run `make inputs` and `make annotation repeatmap` first",
              file=sys.stderr)
        return 1
    for name in OUTPUTS:
        (DERIVED / name).write_text(outputs[name])
        print(f"  wrote {name}: {outputs[name].count(chr(10))} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
