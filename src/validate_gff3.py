#!/usr/bin/env python3
"""Check every GFF3 this repository ships with independent external parsers.

ROADMAP M2: "parses with two independent GFF3 parsers, zero warnings". Three run, each
an unrelated implementation:

1. **GenomeTools `gt gff3validator`** (C), run with `-strict` and `-typecheck so`, so
   every feature type is checked against the Sequence Ontology and every `Parent`
   relationship against SO's part-of graph. CI installs it with `apt-get install
   genometools`. Locally, any `gt` on PATH, or the binary named by `$GT`.
2. **The `gff3` validator from PyPI** (pure Python, unrelated code), which parses each
   file and then runs its own checks: unresolved parents, children outside their
   parent's span, every coordinate inside `##sequence-region` and inside the
   sequence itself (read from the GenBank record), N bases inside a CDS, and phase.
3. **bcbio-gff** (`BCBio.GFF`, built on Biopython), the parser most Python pipelines
   load GFF3 with. It must emit no warning, and the feature tree it builds must hold
   every line of the file, with every `Parent` resolved into a nested feature.

The two Python tools run under the interpreter named by `$GFF3_PYTHON`, by default
`.venv/bin/python` (`python3 -m venv .venv && .venv/bin/pip install gff3==1.0.1
bcbio-gff==0.7.1`; `.venv/` is gitignored). This script itself, like the library, uses
the standard library only.

Every line either tool prints is shown verbatim. The run fails if either tool is
missing (it never reports a tool it could not run as passing), if `gt` prints anything
but "input is valid GFF3", if bcbio-gff warns or loses a line, or if `gff3` logs
anything at all with one stated exception.

**The one exception, and why it is not a flag.** `gff3`'s phase check groups CDS lines
by `Parent` and assumes each parent holds one CDS. GFF3 does not say that: a CDS is
identified by its ID, and the GTF puts two ORFs on one transcript in two places
(`TRL2_CDS` and `TRL2A_ORF` on `TRL2_RNA`, and their `IRL2` copies), so the tool reads
two CDS as one and reports their phases as wrong. A PHASE finding is therefore excused
only where its parent holds more than one CDS ID, and it is still printed. Every phase
in every file is checked against the sequence, by translation, in
`tests/test_m2_gff3.py`, and those four CDS translate to their record's protein.

Run: python3 -I src/validate_gff3.py
"""
from __future__ import annotations

import collections
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
sys.path.insert(0, str(ROOT / "src"))
from build_gff3 import GFF3_FILES, RAW, SEQUENCE_RECORD  # noqa: E402
from build_crosswalk import genbank_sequence  # noqa: E402

GT_SUCCESS = "input is valid GFF3"

DRIVER = r'''
import importlib.metadata, json, logging, sys
records = []
class Keep(logging.Handler):
    def emit(self, record):
        records.append((record.levelname, record.getMessage()))
logger = logging.getLogger("gff3-validation")
logger.setLevel(logging.DEBUG)
logger.addHandler(Keep())
from gff3 import Gff3
gff_path, fasta_path = sys.argv[1], sys.argv[2]
gff = Gff3(logger=logger)
gff.parse(gff_path)
# A file object, not the path: given a path, gff3 1.0.1 opens it in binary mode and
# compares bytes with ">" under Python 3, so it reads no sequence and crashes.
with open(fasta_path) as fasta:
    gff.parse_fasta_external(fasta)
gff.check_unresolved_parents()
gff.check_parent_boundary()
gff.check_reference(sequence_region=True, fasta_external=True)
gff.check_phase()
errors = []
for line in gff.lines:
    for error in line.get("line_errors") or []:
        parents = line.get("attributes", {}).get("Parent") or []
        errors.append({"line": line["line_index"] + 1, "type": error["error_type"],
                       "message": error["message"], "parents": parents})
print(json.dumps({"version": importlib.metadata.version("gff3"),
                  "lines": len(gff.lines), "log": records, "errors": errors}))
'''


BCBIO_DRIVER = r'''
import collections, importlib.metadata, json, logging, sys, warnings
records = []
class Keep(logging.Handler):
    def emit(self, record):
        records.append((record.levelname, record.getMessage()))
logging.getLogger().addHandler(Keep())
logging.getLogger().setLevel(logging.WARNING)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    from BCBio import GFF
    with open(sys.argv[1]) as fh:
        parsed = list(GFF.parse(fh))
def walk(features, depth):
    for feature in features:
        yield feature, depth
        yield from walk(feature.sub_features, depth + 1)
types = collections.Counter()
orphans = 0
for record in parsed:
    for feature, depth in walk(record.features, 0):
        types[feature.type] += 1
        if depth == 0 and "Parent" in feature.qualifiers:
            orphans += 1
print(json.dumps({"version": importlib.metadata.version("bcbio-gff"),
                  "biopython": importlib.metadata.version("biopython"),
                  "records": [r.id for r in parsed], "types": types,
                  "orphans": orphans,
                  "warnings": [str(w.message) for w in caught] + [
                      f"{level} {message}" for level, message in records]}))
'''


def gt_command() -> str | None:
    return os.environ.get("GT") or shutil.which("gt")


def gff3_python() -> pathlib.Path:
    return pathlib.Path(os.environ.get("GFF3_PYTHON") or ROOT / ".venv" / "bin" / "python")


def cds_ids_by_parent(path: pathlib.Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = collections.defaultdict(set)
    for text in path.read_text().splitlines():
        if text.startswith("#") or not text.strip():
            continue
        fields = text.split("\t")
        if fields[2] != "CDS":
            continue
        # Tolerant of a malformed pair: the parsers report that, and this helper only
        # decides which of their phase findings to excuse.
        attrs = dict(pair.split("=", 1) for pair in fields[8].split(";") if "=" in pair)
        for parent in attrs.get("Parent", "").split(","):
            if parent:
                out[parent].add(attrs.get("ID", ""))
    return out


def run_gt(files: list[pathlib.Path]) -> list[str]:
    gt = gt_command()
    if not gt:
        return ["genometools `gt` was not found on PATH or in $GT, so the first parser "
                "did not run. Install it (`apt-get install genometools`, or build it "
                "from source) and run this again."]
    version = subprocess.run([gt, "--version"], capture_output=True, text=True)
    print(f"[gt] {version.stdout.splitlines()[0] if version.stdout else version.stderr}")
    failures = []
    for path in files:
        proc = subprocess.run([gt, "gff3validator", "-strict", "-typecheck", "so",
                               str(path)], capture_output=True, text=True)
        output = (proc.stdout + proc.stderr).strip()
        print(f"[gt] gt gff3validator -strict -typecheck so {path.name}")
        for line in output.splitlines():
            print(f"[gt]   {line}")
        if proc.returncode != 0 or output != GT_SUCCESS:
            failures.append(f"gt: {path.name}: exit {proc.returncode}: {output[:300]}")
    return failures


def run_gff3(files: list[pathlib.Path]) -> list[str]:
    python = gff3_python()
    if not python.exists():
        return [f"{python} does not exist, so the second parser did not run. Create it "
                f"with `python3 -m venv .venv && .venv/bin/pip install gff3==1.0.1`, "
                f"or name another interpreter in $GFF3_PYTHON."]
    failures = []
    with tempfile.TemporaryDirectory() as scratch:
        for path in files:
            source = next(s for s, name in GFF3_FILES.items() if name == path.name)
            seqid = next(line.split()[1] for line in path.read_text().splitlines()
                         if line.startswith("##sequence-region"))
            fasta = pathlib.Path(scratch) / f"{source}.fa"
            sequence = genbank_sequence(RAW / SEQUENCE_RECORD[source])
            fasta.write_text(f">{seqid}\n" + "\n".join(
                sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
            proc = subprocess.run([str(python), "-I", "-c", DRIVER, str(path),
                                   str(fasta)], capture_output=True, text=True)
            if proc.returncode != 0:
                failures.append(f"gff3: {path.name}: the validator crashed: "
                                f"{proc.stderr[-500:]}")
                continue
            report = json.loads(proc.stdout)
            print(f"[gff3 {report['version']}] {path.name}: {report['lines']} lines "
                  f"parsed, {len(report['log'])} log records, "
                  f"{len(report['errors'])} line errors")
            for level, message in report["log"]:
                for line in message.splitlines():
                    print(f"[gff3]   {level} {line}")
            if proc.stderr.strip():
                for line in proc.stderr.strip().splitlines():
                    print(f"[gff3]   stderr {line}")
                failures.append(f"gff3: {path.name}: wrote to stderr")
            shared = cds_ids_by_parent(path)
            excused = [e for e in report["errors"]
                       if e["type"] == "PHASE" and e["parents"]
                       and all(len(shared[p]) > 1 for p in e["parents"])]
            unexcused = [e for e in report["errors"] if e not in excused]
            for e in excused:
                print(f"[gff3]   excused: line {e['line']} PHASE under "
                      f"{e['parents']}, which holds {len(shared[e['parents'][0]])} "
                      f"CDS IDs ({sorted(shared[e['parents'][0]])}); the tool assumes "
                      f"one CDS per parent")
            if unexcused:
                failures.append(f"gff3: {path.name}: {len(unexcused)} findings: "
                                f"{unexcused[:3]}")
            # Every log record must be one of the excused errors, so nothing the tool
            # says at any level (it logs unresolved parents at INFO) passes unread.
            if len(report["log"]) != len(excused):
                failures.append(f"gff3: {path.name}: {len(report['log'])} log records "
                                f"against {len(excused)} excused findings")
    return failures


def run_bcbio(files: list[pathlib.Path]) -> list[str]:
    python = gff3_python()
    if not python.exists():
        return [f"{python} does not exist, so bcbio-gff did not run. Create it with "
                f"`python3 -m venv .venv && .venv/bin/pip install gff3==1.0.1 "
                f"bcbio-gff==0.7.1`, or name another interpreter in $GFF3_PYTHON."]
    failures = []
    for path in files:
        proc = subprocess.run([str(python), "-I", "-c", BCBIO_DRIVER, str(path)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            failures.append(f"bcbio-gff: {path.name}: crashed: {proc.stderr[-500:]}")
            continue
        report = json.loads(proc.stdout)
        lines = collections.Counter(
            text.split("\t")[2] for text in path.read_text().splitlines()
            if text and not text.startswith("#"))
        seqid = next(line.split()[1] for line in path.read_text().splitlines()
                     if line.startswith("##sequence-region"))
        print(f"[bcbio-gff {report['version']}, biopython {report['biopython']}] "
              f"{path.name}: records {report['records']}, "
              f"{sum(report['types'].values())} features in the tree, "
              f"{report['orphans']} unresolved parents, "
              f"{len(report['warnings'])} warnings")
        for warning in report["warnings"]:
            print(f"[bcbio-gff]   {warning}")
        if proc.stderr.strip():
            for line in proc.stderr.strip().splitlines():
                print(f"[bcbio-gff]   stderr {line}")
            failures.append(f"bcbio-gff: {path.name}: wrote to stderr")
        if report["warnings"]:
            failures.append(f"bcbio-gff: {path.name}: {len(report['warnings'])} "
                            f"warnings")
        if report["records"] != [seqid]:
            failures.append(f"bcbio-gff: {path.name}: records {report['records']}")
        if report["orphans"]:
            failures.append(f"bcbio-gff: {path.name}: {report['orphans']} features "
                            f"name a Parent that was not resolved")
        if collections.Counter(report["types"]) != lines:
            failures.append(f"bcbio-gff: {path.name}: the tree holds "
                            f"{report['types']} where the file has {dict(lines)}")
    return failures


def main() -> int:
    files = [DERIVED / name for name in sorted(GFF3_FILES.values())]
    absent = [f.name for f in files if not f.exists()]
    if absent:
        print(f"not built: {absent}; run `make gff3` first", file=sys.stderr)
        return 1
    failures = run_gt(files) + run_bcbio(files) + run_gff3(files)
    if failures:
        print("\nGFF3 validation FAILED:\n  " + "\n  ".join(failures), file=sys.stderr)
        return 1
    print(f"\n{len(files)} GFF3 files: gt and bcbio-gff report nothing; gff3 reports "
          f"nothing beyond the excused phase findings listed above")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
