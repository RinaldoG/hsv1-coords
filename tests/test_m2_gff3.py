"""The M2 GFF3 gate, and the BED rendering of the repeat map.

ROADMAP M2's last criterion is that the annotation parses with two independent GFF3
parsers with zero warnings. Those parsers are external tools and run under
`make validate-gff3`, not here, because the library and this suite are standard library
only. What this file gates is what a parser cannot see:

- every row of the annotation table is in its source's GFF3 exactly once, and nothing
  else is, compared as whole tuples (type, name, strand, every segment, fuzzy ends,
  class, product, notes, provenance) rather than as counts;
- every CDS phase is right AGAINST THE SEQUENCE: each CDS is translated from the
  GFF3's own segments and phases and must reproduce the record's `/translation`, and
  each segment's own phase must land on the record's reading frame;
- the BED files convert back to the repeat map and the diploid pairs exactly, and their
  coordinates, used as 0-based half-open slices, read bases that prove the frame;
- the committed files are what the builder produces now.

The location parser, the GFF3 reader, the BED reader and the translation here are
written independently of `src/build_gff3.py`, so an error in the builder's arithmetic is
not reproduced by the check. The type mapping is imported, because a second typed copy
of a declaration goes stale (see `test_classes_are_from_the_declared_set` in the M2
annotation gate).
"""
from __future__ import annotations

import collections
import csv
import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
TABLE = DERIVED / "annotation_v0.1.0.tsv"

import build_gff3  # noqa: E402

RECORDS = {"NC_001806.2": RAW / "ref_NC_001806.2.gb",
           "BK012101.1": RAW / "alt_BK012101.1.gb",
           "whisnant2020-gtf": RAW / "base_JN555585.1.gb"}
GTF = RAW / "HSV1_Annotation.gtf"
GENBANK_SOURCES = ("NC_001806.2", "BK012101.1")
GTF_SOURCE = "whisnant2020-gtf"


def _require_inputs() -> None:
    """Absent inputs fail in a checkout and skip outside one, as the other gates do."""
    needed = ([TABLE, GTF, DERIVED / "repeat_map.tsv", DERIVED / "diploid_pairs.tsv"]
              + list(RECORDS.values())
              + [DERIVED / name for name in build_gff3.OUTPUTS])
    absent = [str(p.relative_to(ROOT)) for p in needed if not p.exists()]
    if not absent:
        return
    if (ROOT / ".git").exists():
        raise AssertionError(f"{absent} are absent in a git checkout. Run `make all`.")
    raise unittest.SkipTest(f"not a checkout, and these are absent: {absent}")


# ---------------------------------------------------------------- independent readers

def sequence_of(record: pathlib.Path) -> str:
    body = record.read_text().split("\nORIGIN", 1)[1].split("\n//", 1)[0]
    return re.sub(r"[^acgtnACGTN]", "", body).upper()


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE_1 = {a + b + c: AMINO[16 * i + 4 * j + k]
           for i, a in enumerate(BASES) for j, b in enumerate(BASES)
           for k, c in enumerate(BASES)}
STARTS_TABLE_1 = {"ATG", "TTG", "CTG"}


def parse_location(text: str) -> tuple[str, tuple[tuple[int, int], ...], bool, bool]:
    """GenBank location to (strand, segments ascending, fuzzy left, fuzzy right)."""
    strand = "-" if text.startswith("complement(") else "+"
    spans = re.findall(r"(<?)(\d+)\.\.(>?)(\d+)", text)
    segments = tuple(sorted((int(a), int(b)) for _, a, _, b in spans))
    lefts = {int(a) for lt, a, _, _ in spans if lt}
    rights = {int(b) for _, _, gt, b in spans if gt}
    return (strand, segments, segments[0][0] in lefts, segments[-1][1] in rights)


def unescape(value: str) -> str:
    return re.sub(r"%([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), value)


class Line:
    def __init__(self, number: int, text: str):
        fields = text.split("\t")
        self.number, self.fields = number, fields
        (self.seqid, self.source, self.type, start, end, self.score, self.strand,
         self.phase, self.raw_attributes) = fields
        self.start, self.end = int(start), int(end)
        self.attrs = {}
        for pair in self.raw_attributes.split(";"):
            key, _, value = pair.partition("=")
            self.attrs[key] = unescape(value)

    @property
    def id(self) -> str:
        return self.attrs.get("ID", "")


def read_gff3(path: pathlib.Path) -> tuple[list[str], list[Line]]:
    directives, lines = [], []
    for number, text in enumerate(path.read_text().splitlines(), 1):
        if text.startswith("##"):
            directives.append(text)
        elif text.startswith("#") or not text.strip():
            continue
        else:
            lines.append(Line(number, text))
    return directives, lines


def table_rows(source: str) -> list[dict[str, str]]:
    with TABLE.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["source"] == source]


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_bed(path: pathlib.Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()
            if line and not line.startswith(("#", "track", "browser"))]


def gff3_path(source: str) -> pathlib.Path:
    return DERIVED / build_gff3.GFF3_FILES[source]


def expected_type(source: str, row: dict[str, str]) -> str:
    kind = "gtf" if source == GTF_SOURCE else "genbank"
    return build_gff3.so_type_of(kind, row)


def fuzzy_from(line: Line) -> tuple[bool, bool]:
    """The fuzzy ends a GFF3 line declares, each checked against its own coordinate."""
    left = line.attrs.get("start_range")
    right = line.attrs.get("end_range")
    if left is not None:
        assert left == f".,{line.start}", f"line {line.number}: start_range {left}"
    if right is not None:
        assert right == f"{line.end},.", f"line {line.number}: end_range {right}"
    if left is not None or right is not None:
        assert line.attrs.get("partial") == "true", f"line {line.number}: no partial"
    return left is not None, right is not None


def genbank_features(lines: list[Line]) -> dict[str, list[Line]]:
    by_id: dict[str, list[Line]] = collections.defaultdict(list)
    for line in lines:
        by_id[line.id].append(line)
    return by_id


def described(row_or_line, kind: str) -> tuple:
    """The attributes a row and a GFF3 feature must agree on, besides the location."""
    if kind == "row":
        r = row_or_line
        return (r["name"], r["feature_class"], r["product"], r["subtype"],
                r["transcript_id"], r["notes"], r["provenance"])
    a = row_or_line.attrs
    return (a.get("Name", ""), a.get("feature_class", ""), a.get("product", ""),
            a.get("subtype", ""), a.get("transcript_id", ""), a.get("Note", ""),
            a.get("provenance", ""))


# ---------------------------------------------------------------- the gates

class TheFilesAreWellFormedGff3(unittest.TestCase):
    """What a GFF3 parser checks, checked here too so the suite fails without one."""

    @classmethod
    def setUpClass(cls):
        _require_inputs()
        cls.files = {s: read_gff3(gff3_path(s)) for s in RECORDS}

    def test_each_file_names_its_version_and_its_own_sequence(self):
        gtf_seqids = {line.split("\t")[0] for line in GTF.read_text().splitlines()
                      if line and not line.startswith("#")}
        self.assertEqual(len(gtf_seqids), 1)
        seqids = {"NC_001806.2": "NC_001806.2", "BK012101.1": "BK012101.1",
                  GTF_SOURCE: gtf_seqids.pop()}
        for source, (directives, lines) in self.files.items():
            with self.subTest(source=source):
                length = len(sequence_of(RECORDS[source]))
                self.assertEqual(directives[0], "##gff-version 3")
                self.assertIn(f"##sequence-region {seqids[source]} 1 {length}",
                              directives)
                self.assertEqual({line.seqid for line in lines}, {seqids[source]})
                bad = [line.number for line in lines
                       if not 1 <= line.start <= line.end <= length]
                self.assertEqual(bad, [], f"lines outside the sequence: {bad[:5]}")

    def test_every_line_has_nine_columns_a_strand_and_a_mapped_type(self):
        allowed = (set(build_gff3.SO_TYPES.values())
                   | set(build_gff3.REGULATORY_CLASSES.values()))
        for source, (_, lines) in self.files.items():
            with self.subTest(source=source):
                self.assertEqual([l.number for l in lines if len(l.fields) != 9], [])
                self.assertEqual({l.type for l in lines} - allowed, set())
                self.assertEqual({l.strand for l in lines} - {"+", "-"}, set())
                self.assertEqual({l.score for l in lines}, {"."})

    def test_phase_is_on_every_cds_line_and_on_nothing_else(self):
        for source, (_, lines) in self.files.items():
            with self.subTest(source=source):
                bad = [l.number for l in lines
                       if (l.type == "CDS") != (l.phase in ("0", "1", "2"))
                       or (l.type != "CDS" and l.phase != ".")]
                self.assertEqual(bad, [], f"phase wrong on lines {bad[:5]}")

    def test_column_9_is_percent_escaped(self):
        """No value holds a raw `;`, `=`, `&`, `,` or control character.

        None of these attributes is multi-valued, so a raw comma anywhere would be
        read as a list separator. A `%` must start a two-digit escape.
        """
        for source, (_, lines) in self.files.items():
            with self.subTest(source=source):
                bad = []
                for line in lines:
                    for pair in line.raw_attributes.split(";"):
                        key, sep, value = pair.partition("=")
                        if (not sep or not key or re.search(r"[,=&\x00-\x1f]", value)
                                or re.search(r"%(?![0-9A-F]{2})", value)):
                            bad.append((line.number, pair[:40]))
                self.assertEqual(bad, [], f"unescaped attribute values: {bad[:5]}")

    def test_an_id_is_one_feature_and_only_a_cds_spans_lines(self):
        for source, (_, lines) in self.files.items():
            with self.subTest(source=source):
                by_id = genbank_features(lines)
                self.assertNotIn("", by_id, "a line has no ID")
                for fid, members in by_id.items():
                    if len(members) == 1:
                        continue
                    self.assertEqual({m.type for m in members}, {"CDS"}, fid)
                    self.assertEqual(len({m.raw_attributes for m in members}), 1,
                                     f"{fid}: the lines of one feature disagree")
                    self.assertEqual(len({m.strand for m in members}), 1, fid)
                    spans = sorted((m.start, m.end) for m in members)
                    self.assertTrue(all(b < c for (_, b), (c, _) in
                                        zip(spans, spans[1:])),
                                    f"{fid}: segments overlap")

    def test_every_parent_exists_and_contains_its_child(self):
        for source, (_, lines) in self.files.items():
            with self.subTest(source=source):
                by_id = genbank_features(lines)
                problems = []
                for line in lines:
                    parent = line.attrs.get("Parent")
                    if parent is None:
                        continue
                    if parent not in by_id:
                        problems.append((line.number, f"no feature {parent}"))
                        continue
                    span = (min(p.start for p in by_id[parent]),
                            max(p.end for p in by_id[parent]))
                    if not span[0] <= line.start <= line.end <= span[1]:
                        problems.append((line.number, f"outside {parent} {span}"))
                    if {p.strand for p in by_id[parent]} != {line.strand}:
                        problems.append((line.number, f"strand differs from {parent}"))
                self.assertEqual(problems, [])


class EveryRowAppearsOnceAndNothingElseDoes(unittest.TestCase):
    """Whole tuples, as multisets, per source. A count would not see a swap."""

    @classmethod
    def setUpClass(cls):
        _require_inputs()

    def _genbank_tuples(self, source: str) -> tuple[collections.Counter, collections.Counter]:
        rows = collections.Counter()
        for r in table_rows(source):
            strand, segments, left, right = parse_location(r["location"])
            rows[(expected_type(source, r), strand, segments, left, right)
                 + described(r, "row")] += 1
        _, lines = read_gff3(gff3_path(source))
        by_id = genbank_features(lines)
        children: dict[str, list[Line]] = collections.defaultdict(list)
        for line in lines:
            if "Parent" in line.attrs:
                self.assertEqual(line.type, "exon",
                                 f"{source} line {line.number}: a GenBank record states "
                                 f"no parent, so only an exon of a join may have one")
                children[line.attrs["Parent"]].append(line)
        rendered = collections.Counter()
        for fid, members in by_id.items():
            first = members[0]
            if first.type == "exon":
                continue
            if fid in children:
                exons = sorted(children[fid], key=lambda l: l.start)
                self.assertGreater(len(exons), 1, f"{fid} has a single exon")
                self.assertEqual((first.start, first.end),
                                 (exons[0].start, exons[-1].end),
                                 f"{fid} does not span its exons")
                segments = tuple((e.start, e.end) for e in exons)
                for exon in exons:
                    self.assertEqual(exon.strand, first.strand)
                left, right = fuzzy_from(first)
                self.assertEqual((left, right), (fuzzy_from(exons[0])[0],
                                                 fuzzy_from(exons[-1])[1]),
                                 f"{fid}: the fuzzy ends of the feature and its "
                                 f"outer exons disagree")
            else:
                segments = tuple(sorted((m.start, m.end) for m in members))
                left, right = fuzzy_from(first)
            rendered[(first.type, first.strand, segments, left, right)
                     + described(first, "line")] += 1
        return rows, rendered

    def test_each_genbank_row_is_one_gff3_feature(self):
        for source in GENBANK_SOURCES:
            with self.subTest(source=source):
                rows, rendered = self._genbank_tuples(source)
                self.assertGreater(sum(rows.values()), 200)
                self.assertEqual(rows - rendered, collections.Counter(),
                                 "rows of the table missing from the GFF3")
                self.assertEqual(rendered - rows, collections.Counter(),
                                 "features in the GFF3 that are not rows of the table")

    def test_each_gtf_row_is_one_gff3_line(self):
        rows = collections.Counter()
        for r in table_rows(GTF_SOURCE):
            strand, segments, left, right = parse_location(r["location"])
            self.assertEqual(len(segments), 1)
            rows[(expected_type(GTF_SOURCE, r), strand, segments, left, right)
                 + described(r, "row")] += 1
        _, lines = read_gff3(gff3_path(GTF_SOURCE))
        rendered = collections.Counter(
            (l.type, l.strand, ((l.start, l.end),)) + fuzzy_from(l)
            + described(l, "line") for l in lines)
        self.assertEqual(sum(rows.values()), len(GTF.read_text().splitlines()))
        self.assertEqual(rows - rendered, collections.Counter(),
                         "GTF rows of the table missing from the GFF3")
        self.assertEqual(rendered - rows, collections.Counter(),
                         "lines in the GTF's GFF3 that are not rows of the table")

    def test_the_gtf_hierarchy_is_the_one_the_gtf_states(self):
        """Every Parent in the GTF's GFF3 is a `gene_id` or `transcript_id` it states.

        Read from the GTF directly: each line's (type, span, strand, gene_id,
        transcript_id) must be reached by walking the GFF3's Parent links.
        """
        stated = collections.Counter()
        for text in GTF.read_text().splitlines():
            if not text or text.startswith("#"):
                continue
            f = text.split("\t")
            gene = re.search(r'gene_id "([^"]*)"', f[8]).group(1)
            tid = re.search(r'transcript_id "([^"]*)"', f[8])
            stated[(f[2], int(f[3]), int(f[4]), f[6], gene,
                    tid.group(1) if tid and f[2] != "gene" else "")] += 1
        _, lines = read_gff3(gff3_path(GTF_SOURCE))
        by_id = genbank_features(lines)
        walked = collections.Counter()
        for line in lines:
            gene, tid, node = "", "", line
            while True:
                if node.type == "gene":
                    gene = node.attrs["gene_id"]
                elif node.type == "transcript":
                    tid = node.attrs["transcript_id"]
                parent = node.attrs.get("Parent")
                if parent is None:
                    break
                node = by_id[parent][0]
            walked[(line.type, line.start, line.end, line.strand, gene, tid)] += 1
        self.assertEqual(stated - walked, collections.Counter())
        self.assertEqual(walked - stated, collections.Counter())


def transl_excepts(feature: dict) -> dict[tuple[int, ...], str]:
    out = {}
    three = {"Met": "M"}
    for m in re.finditer(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)",
                         feature.get("transl_except", "")):
        a, b = int(m.group(2)), int(m.group(3))
        positions = tuple(range(a, b + 1))
        out[positions if not m.group(1) else positions[::-1]] = three[m.group(4)]
    return out


def translate(sequence: str, positions: list[int], exceptions: dict,
              five_prime_complete: bool, strand: str) -> str:
    out = []
    for i in range(0, len(positions) - 2, 3):
        codon_at = tuple(positions[i:i + 3])
        bases = "".join(sequence[p - 1] for p in codon_at)
        if strand == "-":
            bases = bases.translate(COMPLEMENT)
        if codon_at in exceptions:
            out.append(exceptions[codon_at])
        elif i == 0 and five_prime_complete and bases in STARTS_TABLE_1:
            out.append("M")
        else:
            out.append(TABLE_1.get(bases, "X"))
    return "".join(out)


def cds_translation_problems(lines: list[Line], sequence: str, expected: str,
                             exceptions: dict) -> list[str]:
    """Translate one CDS from its GFF3 lines alone, and compare with `expected`.

    Two checks. The whole CDS, read from the first segment's phase, must give
    `expected`. And each later segment, read from ITS OWN phase, must give the stretch of
    `expected` its position says it codes, so a phase that is wrong on an inner segment
    is caught even though the whole-CDS reading never uses it.
    """
    strand = lines[0].strand
    reading = sorted(lines, key=lambda l: l.start, reverse=(strand == "-"))

    def bases_of(line):
        return (list(range(line.start, line.end + 1)) if strand == "+"
                else list(range(line.end, line.start - 1, -1)))

    first = reading[0]
    five_prime = "start_range" if strand == "+" else "end_range"
    complete = five_prime not in first.attrs
    positions = [p for line in reading for p in bases_of(line)]
    protein = translate(sequence, positions[int(first.phase):], exceptions, complete,
                        strand)
    problems = []
    if protein.endswith("*") and len(protein) - 1 == len(expected):
        protein = protein[:-1]
    if protein != expected:
        problems.append(f"translation differs: {protein[:30]}... against "
                        f"{expected[:30]}...")
    before = 0
    for line in reading:
        offset = before + int(line.phase) - int(first.phase)
        before += line.end - line.start + 1
        if line is first:
            continue
        if offset % 3:
            problems.append(f"line {line.number}: phase {line.phase} puts a codon "
                            f"boundary off the reading frame")
            continue
        own = translate(sequence, bases_of(line)[int(line.phase):], exceptions, False,
                        strand)
        start = offset // 3
        if start + len(own) == len(expected) + 1 and own.endswith("*"):
            own = own[:-1]
        if own != expected[start:start + len(own)]:
            problems.append(f"line {line.number}: read from its phase {line.phase} "
                            f"it gives {own[:20]}..., the record has "
                            f"{expected[start:start + 20]}...")
    return problems


class EveryPhaseIsRightAgainstTheSequence(unittest.TestCase):
    """CDS phases checked by translation, never by the builder's arithmetic."""

    @classmethod
    def setUpClass(cls):
        _require_inputs()
        from extract_names import parse_genbank
        cls.record_cds = {s: [f for f in parse_genbank(RECORDS[s]) if f["_type"] == "CDS"]
                          for s in GENBANK_SOURCES}
        cls.sequences = {s: sequence_of(p) for s, p in RECORDS.items()}

    def test_every_genbank_cds_reproduces_its_records_translation(self):
        for source in GENBANK_SOURCES:
            with self.subTest(source=source):
                by_protein = {f["protein_id"]: f for f in self.record_cds[source]}
                _, lines = read_gff3(gff3_path(source))
                cds = collections.defaultdict(list)
                for line in lines:
                    if line.type == "CDS":
                        cds[line.attrs["protein_id"]].append(line)
                translated = [f for f in self.record_cds[source]
                              if f.get("translation")]
                self.assertEqual(sorted(cds), sorted(by_protein),
                                 "the GFF3's CDS are not the record's CDS")
                problems = {}
                for protein_id, members in cds.items():
                    feature = by_protein[protein_id]
                    expected = re.sub(r"\s", "", feature.get("translation", ""))
                    if not expected:
                        continue
                    _, segments, _, _ = parse_location(feature["_location"])
                    self.assertEqual(tuple(sorted((l.start, l.end) for l in members)),
                                     segments, protein_id)
                    found = cds_translation_problems(
                        members, self.sequences[source], expected,
                        transl_excepts(feature))
                    if found:
                        problems[protein_id] = found
                self.assertEqual(problems, {},
                                 f"{len(problems)} of {len(translated)} CDS do not "
                                 f"translate to their record's protein")
                self.assertEqual(len(translated), len(cds))
                print(f"\n    {source}: {len(translated)} CDS translated from their "
                      f"GFF3 segments and phases, all equal to /translation")

    def test_every_gtf_cds_matching_a_record_cds_reproduces_its_translation(self):
        """The GTF gives no translation, so the BK012101.1 record stands in for it.

        Both are the 2020 re-annotation of the same sequence. Each GTF CDS whose strand
        and segments are exactly those of a record CDS must translate to that CDS's
        `/translation` from the GTF's own phases. Every GTF CDS must also read through
        without a stop codon before its last codon from those phases.
        """
        sequence = self.sequences[GTF_SOURCE]
        self.assertEqual(sequence, self.sequences["BK012101.1"],
                         "JN555585.1 and BK012101.1 no longer hold one sequence, so "
                         "the record's translations cannot stand in for the GTF")
        record = collections.defaultdict(list)
        for f in self.record_cds["BK012101.1"]:
            strand, segments, _, _ = parse_location(f["_location"])
            record[(strand, segments)].append(f)
        _, lines = read_gff3(gff3_path(GTF_SOURCE))
        cds = collections.defaultdict(list)
        for line in lines:
            if line.type == "CDS":
                cds[line.id].append(line)
        matched, problems = 0, {}
        for fid, members in sorted(cds.items()):
            key = (members[0].strand, tuple(sorted((l.start, l.end) for l in members)))
            candidates = record.get(key, [])
            if candidates:
                matched += 1
                found = []
                for feature in candidates:
                    found = cds_translation_problems(
                        members, sequence, re.sub(r"\s", "", feature["translation"]),
                        transl_excepts(feature))
                    if not found:
                        break
                if found:
                    problems[fid] = found
            strand = members[0].strand
            reading = sorted(members, key=lambda l: l.start, reverse=(strand == "-"))
            positions = [p for l in reading for p in
                         (range(l.start, l.end + 1) if strand == "+"
                          else range(l.end, l.start - 1, -1))]
            protein = translate(sequence, positions[int(reading[0].phase):], {}, True,
                                strand)
            if "*" in protein[:-1]:
                problems.setdefault(fid, []).append(
                    f"a stop codon at residue {protein.index('*') + 1} of "
                    f"{len(protein)}")
        self.assertEqual(problems, {}, f"{len(problems)} GTF CDS fail: "
                                       f"{dict(list(problems.items())[:3])}")
        # Measured on 17 November 2026: every one of the GTF's CDS has a record CDS
        # with the identical strand and segments, so every GTF phase is checked by
        # translation and none only by the weaker stop-codon read.
        self.assertEqual(matched, len(cds),
                         f"only {matched} of {len(cds)} GTF CDS have a record CDS to "
                         f"translate against")
        print(f"\n    {GTF_SOURCE}: {matched} of {len(cds)} CDS match a BK012101.1 CDS "
              f"exactly and translate to its /translation; all {len(cds)} read "
              f"without an internal stop")


class TheBedFilesKeepTheFrame(unittest.TestCase):
    """BED is 0-based half-open; everything else here is 1-based inclusive."""

    STRAND = {"forward": "+", "reverse": "-", "single": ".", "unequal_length": "."}

    @classmethod
    def setUpClass(cls):
        _require_inputs()
        cls.reference = sequence_of(RECORDS["NC_001806.2"])
        cls.bed = read_bed(DERIVED / build_gff3.REPEAT_BED)
        cls.bedpe = read_bed(DERIVED / build_gff3.PAIRS_BEDPE)

    def test_the_bed_converts_back_to_the_repeat_map(self):
        rows = collections.Counter(
            (r["name"], int(r["start"]), int(r["end"]), r["kind"], r["orientation"],
             self.STRAND[r["orientation"]])
            for r in read_tsv(DERIVED / "repeat_map.tsv"))
        back = collections.Counter(
            (f[3], int(f[1]) + 1, int(f[2]), f[6], f[7], f[5]) for f in self.bed)
        self.assertEqual({f[0] for f in self.bed}, {"NC_001806.2"})
        self.assertEqual(rows - back, collections.Counter(), "map rows missing")
        self.assertEqual(back - rows, collections.Counter(), "BED rows not in the map")

    def test_the_bedpe_converts_back_to_the_diploid_pairs(self):
        rows = collections.Counter(
            (p["locus"], int(p["copy_1_start"]), int(p["copy_1_end"]),
             int(p["copy_2_start"]), int(p["copy_2_end"]), p["repeat_pair"],
             p["feature_type"], p["strand_relationship"], p["offset_delta"],
             p["copy_1_name"], p["copy_2_name"])
            for p in read_tsv(DERIVED / "diploid_pairs.tsv"))
        back = collections.Counter(
            (f[6], int(f[1]) + 1, int(f[2]), int(f[4]) + 1, int(f[5]), f[10], f[11],
             f[12], f[13], f[14], f[15]) for f in self.bedpe)
        self.assertEqual(rows - back, collections.Counter(), "pairs missing")
        self.assertEqual(back - rows, collections.Counter(), "BEDPE rows not pairs")
        strands = {(f[12] == "inverted", f[8], f[9]) for f in self.bedpe}
        self.assertEqual(strands - {(True, "+", "-"), (False, ".", ".")}, set())

    def test_the_bedpe_slices_of_each_inverted_pair_are_reverse_complements(self):
        """The frame proof. A start not shifted by one reads the wrong bases.

        Only pairs of equal length (`offset_delta` 0): `RL reiteration set 1` and `RS
        reiteration set 1` are inverted over their shared length and differ by 1 bp,
        so their full slices cannot be reverse complements whatever the frame.
        """
        checked = []
        for f in self.bedpe:
            if f[8] != "+" or f[9] != "-" or f[13] != "0":
                continue
            first = self.reference[int(f[1]):int(f[2])]
            second = self.reference[int(f[4]):int(f[5])]
            self.assertEqual(first, revcomp(second),
                             f"{f[6]} {f[1]}-{f[2]} / {f[4]}-{f[5]} read through the "
                             f"BED coordinates are not reverse complements")
            checked.append(f[6])
        measurable = [p for p in read_tsv(DERIVED / "diploid_pairs.tsv")
                      if p["strand_relationship"] == "inverted"
                      and p["offset_delta"] == "0"]
        self.assertEqual(len(checked), len(measurable),
                         "an inverted pair of equal length was not read through the "
                         "BEDPE, so its strands were not written as + and -")
        self.assertGreater(len(checked), 10)

    def test_the_bed_slices_of_every_measured_copy_agree_in_orientation(self):
        """Every same-named pair of + and - copies in the BED reads as reverse
        complements, and every pair of + copies reads identically."""
        copies = collections.defaultdict(list)
        for f in self.bed:
            if f[5] in "+-":
                copies[f[3]].append((f[5], self.reference[int(f[1]):int(f[2])]))
        compared = 0
        for name, found in copies.items():
            anchor = [seq for strand, seq in found if strand == "+"]
            self.assertTrue(anchor, f"{name} has no forward copy")
            for strand, seq in found:
                expected = anchor[0] if strand == "+" else revcomp(anchor[0])
                self.assertEqual(seq, expected, f"{name}: a {strand} copy read through "
                                                f"the BED does not match")
                compared += 1
        self.assertEqual(compared, sum(1 for f in self.bed if f[5] in "+-"))
        self.assertGreater(compared, 15)


class TheCommittedFilesAreWhatTheBuilderWrites(unittest.TestCase):
    """Stale or hand-edited output fails here rather than only in CI's rebuild diff."""

    @classmethod
    def setUpClass(cls):
        _require_inputs()

    def test_every_output_is_current_and_the_build_is_deterministic(self):
        first = build_gff3.render()
        second = build_gff3.render()
        self.assertEqual(first, second, "two renders in one process differ")
        self.assertEqual(sorted(first), sorted(build_gff3.OUTPUTS))
        for name in build_gff3.OUTPUTS:
            with self.subTest(file=name):
                self.assertEqual((DERIVED / name).read_text(), first[name],
                                 f"{name} is not what `make gff3` writes now")

    def test_an_unmapped_feature_type_stops_the_build(self):
        row = {"name": "x", "feature_type": "mat_peptide", "subtype": ""}
        with self.assertRaises(KeyError):
            build_gff3.so_type_of("genbank", row)
        with self.assertRaises(KeyError):
            build_gff3.so_type_of("genbank", dict(row, feature_type="regulatory",
                                                  subtype="enhancer"))
        with self.assertRaises(build_gff3.UnsupportedLocation):
            build_gff3.parse_location("order(1..10,20..30)")


if __name__ == "__main__":
    unittest.main(verbosity=2)
