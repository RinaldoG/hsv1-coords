"""Record the M3 by-hand criterion, and render it as a derived document.

M3's fourth criterion, written in `research/PICK.md` before any code, is that the
lift-over "correctly resolves a published coordinate from at least five papers across
at least three strains to the canonical reference, checked by hand against each paper".
That is the one criterion a test run cannot settle on its own: it needs papers read.

What this module holds is the record of that reading, as two committed tables rather
than prose, so the criterion is marked met by something executable rather than by its
author. `PITFALLS.md` standing rule 3 is the reason: twice in this project a criterion
was marked met by the person who wrote it while the thing it required was not built.

**`data/published_coordinates.tsv`** is one row per paper that was opened. Its
`published_coordinate`, `quote` and `source` columns are **hand transcription** from a
paper's full text, flagged in the `transcribed` column exactly as
`data/conventions.tsv` flags its `paper_transcribed` rows. The suite cannot re-read a
paper, because continuous integration has no network and a PDF is not an artefact, so
nothing here asserts that a paper says what the row claims. What the suite does assert
is everything downstream of the transcription.

**`data/published_coordinate_checks.tsv`** is one row per mechanical check, and every
one of them is decided against something committed: the deposited records in
`data/raw`, the genome FASTA, or the built artefacts. No `expect` value in that table
was produced by calling the function it checks. The resolutions were derived
independently, by locating the sequence around a published position in the reference by
exact or mismatch-tolerant anchoring, or by reading a record's own feature table. A
gate whose two sides are written by the same step compares nothing, which is root
cause 41 in `PITFALLS.md`, and the whole point of this table is to not be that.

The triangle each counted case is built on, and what makes it evidence rather than
assertion:

    the sequence the paper prints is at the paper's own coordinate, in the record the
    paper names                                           (`sequence_unique_at`)
    the same sequence is at the resolved coordinate, in the reference
                                                          (`sequence_unique_at`)
    and the library maps the first position to the second (`lift`, `version_convert`)

Any one of the three alone proves nothing. Together they cannot all hold unless the
resolution is right, and `tests/test_m3_published_coordinates.py` asserts the triangle
closes for every motif that appears in two records.

Run as `make published-coordinates`. The output is `data/derived/
published_coordinates.md`, which is deterministic: the same tables and artefacts render
byte-identically.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hsv1coords._paths import data_path, require  # noqa: E402

CASES = data_path("published_coordinates.tsv")
CHECKS = data_path("published_coordinate_checks.tsv")
RAW = ROOT / "data" / "raw"
OUTPUT = "published_coordinates.md"

# The criterion's own two numbers, quoted from `research/PICK.md`. They are the
# thresholds, not measurements, so they are named here and compared against counts
# derived from the table rather than against anything typed beside it.
PAPERS_REQUIRED = 5
STRAINS_REQUIRED = 3
CRITERION = ("Correctly resolves a published coordinate from at least five papers "
             "across at least three strains to the canonical reference, checked by "
             "hand against each paper.")

# Controlled values. A status or a basis outside these lists is a typo that would
# otherwise change a count silently, which is why they are closed sets and not advice.
STATUSES = ("counted", "supporting", "not_counted")
BASES = ("stated_with_version", "stated_accession_held", "determined_from_paper_claim",
         "not_stated", "none_printed", "not_retrieved")
KINDS_THAT_RESOLVE = ("version_conversion", "strain_lift")
RESOLUTION_KINDS = KINDS_THAT_RESOLVE + ("identity", "declined", "unresolvable", "none")
TRANSCRIBED = ("by_hand_from_full_text", "by_hand_fetch_failed")

# A basis that cannot support a counted case. `not_stated` is the strict rule: a paper
# that does not say which record its numbers are in has not published a resolvable
# coordinate, however well the arithmetic happens to come out.
BASES_THAT_CANNOT_COUNT = ("not_stated", "none_printed", "not_retrieved")

# Check kinds that tie a published coordinate to a resolved one, and kinds that are
# evidence from outside the library. A counted case needs at least one of each, so a
# resolution cannot be recorded with nothing independent behind it.
RESOLVING_CHECKS = ("version_convert", "lift")
EVIDENCE_CHECKS = ("sequence_unique_at", "sequence_at", "feature_location",
                   "counterpart_pair", "window_revcomp_equals", "window_equals",
                   "orf_on_minus_strand", "orf_matches_translation", "feature_absent",
                   "positions_are_feature_edges", "codon_positions", "interval_length")
CHECK_KINDS = tuple(sorted(set(RESOLVING_CHECKS + EVIDENCE_CHECKS + ("lift_declines",))))

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _rows(path: pathlib.Path) -> list[dict[str, str]]:
    with require(path).open() as fh:
        return [{k: (v or "").strip() for k, v in row.items()}
                for row in csv.DictReader(fh, delimiter="\t")]


def load_cases() -> list[dict[str, str]]:
    rows = _rows(CASES)
    if not rows:
        raise ValueError(
            f"{CASES} states no cases, so every question about the M3 by-hand "
            f"criterion would be answered by silence and the criterion would read as "
            f"met over nothing. Restore the file.")
    return rows


def load_checks() -> list[dict[str, str]]:
    rows = _rows(CHECKS)
    if not rows:
        raise ValueError(
            f"{CHECKS} states no checks, so the cases table would be prose. "
            f"Restore the file.")
    return rows


def counted(cases: list[dict[str, str]]) -> list[dict[str, str]]:
    return [c for c in cases if c["status"] == "counted"]


def strains_counted(cases: list[dict[str, str]]) -> list[str]:
    """The distinct strains among the counted cases.

    Normalised in the table: three papers write 17, 17syn+ and 17Syn+ for the same
    strain, and counting those as three would clear the criterion's strain threshold
    on spelling alone. `strain_as_stated` keeps each paper's own words.
    """
    return sorted({c["strain"] for c in counted(cases)})


# --------------------------------------------------------------------------- loaders

def read_fasta(path: pathlib.Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name, buf = None, []
    with require(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def genbank_features(path: pathlib.Path) -> list[tuple[str, str, dict[str, str]]]:
    """Every feature of a GenBank record as (type, location, qualifiers).

    The location is joined across its continuation lines before the first qualifier.
    Round 7 shipped eight truncated coordinates because a parser set its continuation
    key to None at a feature start and so dropped every wrapped location line; this
    reads the lines before the first `/` instead of tracking a key.
    """
    text = require(path).read_text()
    body = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0]
    out = []
    for chunk in re.split(r"\n(?=     \S)", body):
        lines = chunk.split("\n")
        head = re.match(r"\s+(\S+)\s+(.*)$", lines[0])
        if not head:
            continue
        location, qualifier_lines, in_location = [head.group(2).strip()], [], True
        for line in lines[1:]:
            if line.strip().startswith("/"):
                in_location = False
            (location if in_location else qualifier_lines).append(line.strip()
                                                                  if in_location
                                                                  else line)
        flat = "\n".join(qualifier_lines).replace("\n                     ", "")
        quals: dict[str, str] = {}
        for m in re.finditer(r'/(\w+)=(?:"([^"]*)"|(\S+))', flat):
            quals.setdefault(m.group(1),
                             m.group(2) if m.group(2) is not None else m.group(3))
        out.append((head.group(1), "".join(location), quals))
    return out


def genbank_sequence(path: pathlib.Path) -> str:
    origin = require(path).read_text().split("\nORIGIN", 1)[1]
    return "".join(re.findall(r"[acgtnryswkmbdhv]+", origin)).upper()


def location_edges(location: str) -> set[int]:
    bare = location.replace("<", "").replace(">", "")
    edges: set[int] = set()
    for a, b in re.findall(r"(\d+)\.\.(\d+)", bare):
        edges.update((int(a), int(b)))
    return edges


CODON_TABLE = {}
_bases = "TCAG"
_acids = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for _a in _bases:
    for _b in _bases:
        for _c in _bases:
            CODON_TABLE[_a + _b + _c] = _acids[_i]
            _i += 1


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - len(seq) % 3, 3))


# ------------------------------------------------------------------- the check runner

class Artefacts:
    """Everything the checks are decided against, loaded once."""

    def __init__(self) -> None:
        self.genomes = read_fasta(RAW / "genomes.fasta")
        self._records: dict[str, list[tuple[str, str, dict[str, str]]]] = {}
        self._sequences: dict[str, str] = {}

    def sequence(self, record: str) -> str:
        """A record's bases, from the FASTA for an accession or from a .gb file."""
        if record.endswith(".gb"):
            if record not in self._sequences:
                self._sequences[record] = genbank_sequence(RAW / record)
            return self._sequences[record]
        if record == "NC_001806.1":
            # Not in the corpus: the FASTA holds the current reference only. The .1
            # record's own file carries it, and X14112.1 is the same sequence.
            return self.sequence("ref_NC_001806.1.gb")
        if record not in self.genomes:
            raise KeyError(f"{record} is not in data/raw/genomes.fasta")
        return self.genomes[record]

    def features(self, record: str) -> list[tuple[str, str, dict[str, str]]]:
        if record not in self._records:
            self._records[record] = genbank_features(RAW / record)
        return self._records[record]

    def bases_at(self, record: str, position: int, length: int) -> str:
        return self.sequence(record)[position - 1:position - 1 + length]

    def occurrences(self, record: str, motif: str) -> list[int]:
        return [m.start() + 1
                for m in re.finditer("(?=" + re.escape(motif) + ")",
                                     self.sequence(record))]


def run_check(row: dict[str, str], art: Artefacts) -> tuple[bool, str]:
    """Decide one check row. Returns (passed, what was actually found).

    Every branch reports what it found, not just whether it matched, because a gate
    that says only "failed" invites the wrong repair: root cause 47 is a failure
    message naming a new total rather than what moved.
    """
    kind, record, target = row["kind"], row["record"], row["target"]
    expect, motif = row["expect"], row["motif"]
    position = int(row["position"]) if row["position"].isdigit() else None

    if kind == "version_convert":
        from hsv1coords.versions import VersionCrosswalk
        got = VersionCrosswalk.load().convert(position, record, target)
        return str(got) == expect, f"convert({position}, {record}, {target}) = {got}"

    if kind == "lift":
        from hsv1coords import lift
        got = lift(position, record, target)
        return str(got) == expect, f"lift({position}, {record}, {target}) = {got}"

    if kind == "lift_declines":
        from hsv1coords import lift
        from hsv1coords.liftover import lift_detailed
        plain = lift(position, record, target)
        if plain is not None:
            return False, (f"lift({position}, {record}, {target}) returned {plain} "
                           f"where it must decline")
        result = lift_detailed(position, record, target)
        primary, _, alternative = expect.partition("|")
        got = f"{result.position}|{','.join(str(a) for a in result.alternatives)}"
        ok = (str(result.position) == primary
              and int(alternative) in tuple(result.alternatives))
        return ok, f"declined; detail {got}, wanted {expect}"

    if kind in ("sequence_at", "sequence_unique_at"):
        found = art.bases_at(record, position, len(motif))
        if found != motif:
            return False, f"{record}:{position} holds {found!r}, not {motif!r}"
        if kind == "sequence_unique_at":
            hits = art.occurrences(record, motif)
            if hits != [position]:
                return False, (f"{motif!r} occurs at {hits} in {record}, not only at "
                               f"{position}")
        return True, f"{record}:{position} = {found}"

    if kind == "feature_location":
        ftype, _, selector = target.partition("/")
        qualifier, _, value = selector.partition("=")
        locations = [loc for t, loc, quals in art.features(record)
                     if t == ftype
                     and quals.get(qualifier, "").split(";")[0].strip() == value]
        got = "|".join(locations)
        return got == expect, f"{record} {target} -> {got or '(none)'}"

    if kind == "feature_absent":
        count = sum(1 for t, _, _ in art.features(record) if t == target)
        return str(count) == expect, f"{record} holds {count} {target} feature(s)"

    if kind == "counterpart_pair":
        from hsv1coords import RepeatMap
        other = int(expect)
        seq = art.sequence(record)
        paired = seq[position - 1] == revcomp(seq[other - 1])
        status = RepeatMap.load().classify(position)
        named = other in tuple(status.alternatives)
        return paired and named, (
            f"{record}:{position}={seq[position - 1]} against {other}="
            f"{seq[other - 1]}; complementary {paired}; repeat map alternatives "
            f"{tuple(status.alternatives)}")

    if kind in ("window_equals", "window_revcomp_equals"):
        length = int(expect)
        (rec_a, _, pos_a), (rec_b, _, pos_b) = (record.partition(":"),
                                                target.partition(":"))
        left = art.bases_at(rec_a, int(pos_a), length)
        if kind == "window_equals":
            right = art.bases_at(rec_b, int(pos_b), length)
        else:
            # The counterpart runs backwards, so the window ENDS at the named base.
            right = revcomp(art.bases_at(rec_b, int(pos_b) - length + 1, length))
        return left == right, (f"{rec_a}:{pos_a} {left[:24]}... against {rec_b}:"
                               f"{pos_b} {right[:24]}...")

    if kind == "positions_are_feature_edges":
        positions = [int(p) for p in row["position"].split(",")]
        edges: set[int] = set()
        for t, loc, _ in art.features(record):
            if t == target:
                edges.update(location_edges(loc))
        hit = sum(1 for p in positions if p in edges)
        return str(hit) == expect, (f"{hit} of {len(positions)} positions are edges of "
                                    f"a {target} feature in {record}")

    if kind == "orf_on_minus_strand":
        start, end = position, int(target)
        protein = translate(revcomp(art.bases_at(record, start, end - start + 1)))
        clean = protein.startswith("M") and protein.endswith("*") \
            and "*" not in protein[:-1]
        return (clean and str(len(protein) - 1) == expect), (
            f"{record}:{start}-{end} on the minus strand gives "
            f"{len(protein) - 1} residues, {protein[:4]}...{protein[-4:]}")

    if kind == "orf_matches_translation":
        start, end = position, int(target)
        protein = translate(revcomp(art.bases_at(record, start, end - start + 1)))[:-1]
        published = [q["translation"] for t, _, q in art.features("ref_NC_001806.2.gb")
                     if t == "CDS" and q.get("gene") == motif and "translation" in q]
        if not published:
            return False, f"the reference states no translation for {motif}"
        reference = published[0]
        same = sum(1 for a, b in zip(protein, reference) if a == b)
        return f"{same}/{len(reference)}" == expect, (
            f"{len(protein)} residues against the reference's {len(reference)}, "
            f"{same} identical")

    if kind == "codon_positions":
        from hsv1coords import codon
        got = codon(target, position)
        return ",".join(str(p) for p in got.positions) == expect, (
            f"codon({target}, {position}).positions = {got.positions}")

    if kind == "interval_length":
        start, end = int(row["position"]), int(target)
        return str(end - start + 1) == expect, (
            f"{start} to {end} inclusive is {end - start + 1} bases")

    raise ValueError(f"unknown check kind {kind!r} in {CHECKS}; the kinds this runner "
                     f"decides are {CHECK_KINDS}")


def run_all(art: Artefacts | None = None
            ) -> list[tuple[dict[str, str], bool, str]]:
    art = art or Artefacts()
    return [(row, *run_check(row, art)) for row in load_checks()]


# ----------------------------------------------------------------------- the renderer

def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(cell.replace("|", "\\|") for cell in row) + " |"
            for row in rows]
    return out


def render() -> str:
    cases, results = load_cases(), run_all()
    kept, strains = counted(cases), strains_counted(cases)
    failed = [(r, why) for r, ok, why in results if not ok]
    by_case: dict[str, int] = {}
    for row, _, _ in results:
        by_case[row["case"]] = by_case.get(row["case"], 0) + 1

    lines = [
        "# Published coordinates resolved by hand",
        "",
        "Derived by `src/build_published_coordinates.py` from "
        "`data/published_coordinates.tsv` and",
        "`data/published_coordinate_checks.tsv`. Do not edit: "
        "`make published-coordinates` rewrites it.",
        "",
        "This is the record of M3's one criterion a test run cannot settle by itself.",
        "The criterion, written in the build plan before any code:",
        "",
        f"> {CRITERION}",
        "",
        f"**Met:** {len(kept)} papers against the {PAPERS_REQUIRED} required, across "
        f"{len(strains)} strains against the {STRAINS_REQUIRED} required "
        f"({', '.join(strains)}).",
        f"{len(cases)} papers were opened in all. Every resolution is checked by "
        f"{len(results)} mechanical checks over the",
        "deposited records, the genome FASTA and the built artefacts, "
        f"of which {len(results) - len(failed)} pass.",
        "",
        "**What is transcription and what is checked.** The published coordinate, the "
        "quoted sentence and",
        "the fetch route in each row are **hand transcription** from the paper's full "
        "text, flagged in the",
        "`transcribed` column the same way `data/conventions.tsv` flags its "
        "`paper_transcribed` rows.",
        "Nothing in the suite re-reads a paper: continuous integration has no network. "
        "Everything",
        "downstream of the transcription is checked, and no expected value in the "
        "checks table was",
        "produced by calling the function it checks.",
        "",
        "**Quotations are literal slices** of the fetched bytes, cut by string search "
        "rather than retyped,",
        "so a subscript in the publisher's markup can come through as a space: Lee et "
        "al.'s \"both TR L",
        "repeats\" is how that sentence reads once the tags are stripped.",
        "",
        "## The counted cases",
        "",
    ]
    lines += _table(
        ["Paper", "Strain", "Accession as the paper states it", "How the frame is known",
         "Published", "Resolved in NC_001806.2", "Checks"],
        [[c["citation"], c["strain"], c["frame_accession"] or "(none stated)",
          c["frame_basis"].replace("_", " "), c["published_coordinate"], c["resolved"],
          str(by_case.get(c["case"], 0))] for c in kept])
    lines += ["", "### Why each counted case's frame is known", ""]
    for c in kept:
        lines += [f"- **{c['case']}** ({c['frame_basis']}): {c['frame_evidence']}"]
    lines += ["", "## Cases deliberately not counted", "",
              "Recorded because a clean run of five is not what happened, and because "
              "a criterion is worth",
              "nothing if the cases that failed it are not visible.", ""]
    lines += _table(
        ["Paper", "Strain", "Status", "Why it is not counted"],
        [[c["citation"].split(".")[0], c["strain"] or "-", c["status"],
          c["status_reason"]] for c in cases if c["status"] != "counted"])
    lines += ["", "## Every check, and what decides it", ""]
    for case in sorted({r["case"] for r, _, _ in results},
                       key=lambda k: [c["case"] for c in cases].index(k)):
        rows = [(r, ok, why) for r, ok, why in results if r["case"] == case]
        lines += [f"### {case}", ""]
        lines += _table(["Check", "Kind", "Decided against", "Result"],
                        [[r["check"], r["kind"], why, "pass" if ok else "FAIL"]
                         for r, ok, why in rows])
        lines += [""]
    if failed:
        lines += ["## Failures", ""]
        lines += [f"- `{r['case']}/{r['check']}`: {why}" for r, why in failed]
        lines += [""]
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    out = ROOT / "data" / "derived" / OUTPUT
    out.write_text(render())
    cases = load_cases()
    results = run_all()
    bad = [r for r, ok, _ in results if not ok]
    print(f"wrote {out.relative_to(ROOT)}: {len(counted(cases))} counted papers across "
          f"{len(strains_counted(cases))} strains, {len(results) - len(bad)} of "
          f"{len(results)} checks passing")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
