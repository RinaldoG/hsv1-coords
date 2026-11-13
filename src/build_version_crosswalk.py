#!/usr/bin/env python3
"""Map coordinates between the two versions of the reference accession (roadmap M6).

`NC_001806.1` is 152,261 bp and `NC_001806.2` is 152,222 bp. Published work is defined
against both, and the difference is **not** a constant shift, so a coordinate can be
wrong purely because two papers cited two versions of one accession, with no repeat
involved at all. Until this existed, nothing in this repository could detect that.

Built the same way as the strain crosswalk: align with minimap2, reduce to ungapped
blocks, chain to monotonic, and decline where there is no unique answer.

Output:
  data/derived/version_crosswalk.tsv   ungapped blocks, .1 query against .2 target
  data/derived/version_crosswalk.md    where the two versions diverge, and by how much

Run: python3 -I src/build_version_crosswalk.py
"""
from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"

sys.path.insert(0, str(ROOT / "src"))
from build_crosswalk import (  # noqa: E402
    PRESET, blocks_from_cs, drop_ambiguous, drop_non_monotonic, genbank_sequence)

OLD, NEW = "NC_001806.1", "NC_001806.2"


def main() -> int:
    if not shutil.which("minimap2"):
        print("minimap2 not found; it is a declared dependency", file=sys.stderr)
        return 1
    old_path, new_path = RAW / "ref_NC_001806.1.gb", RAW / "ref_NC_001806.2.gb"
    if not old_path.exists():
        print(f"missing {old_path}, run `make inputs` first", file=sys.stderr)
        return 1

    old_seq, new_seq = genbank_sequence(old_path), genbank_sequence(new_path)
    print(f"  {OLD}: {len(old_seq):,} bp")
    print(f"  {NEW}: {len(new_seq):,} bp  (difference {len(old_seq) - len(new_seq):+d})")

    tmp = pathlib.Path(tempfile.mkdtemp())
    (tmp / "old.fa").write_text(f">{OLD}\n{old_seq}\n")
    (tmp / "new.fa").write_text(f">{NEW}\n{new_seq}\n")
    proc = subprocess.run(
        ["minimap2", "-cx", PRESET, "--cs", "-t", "4", "--secondary=no",
         str(tmp / "new.fa"), str(tmp / "old.fa")],
        capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"minimap2 failed: {proc.stderr[-300:]}", file=sys.stderr)
        return 1

    blocks: list[tuple[int, int, int]] = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if f[4] != "+":
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs:
            blocks.extend(blocks_from_cs(int(f[2]), int(f[7]), cs))
    blocks.sort()
    blocks, dropped_amb = drop_ambiguous(blocks)
    blocks, dropped_mono = drop_non_monotonic(blocks)

    out = DERIVED / "version_crosswalk.tsv"
    with out.open("w") as fh:
        fh.write("query_start\ttarget_start\tlength\n")
        for qs, ts, ln in blocks:
            fh.write(f"{qs}\t{ts}\t{ln}\n")

    # The offsets, which is the point: a reader needs to see that there is no one number.
    offsets: list[tuple[int, int, int]] = []   # (.1 position, .2 position, offset)
    for qs, ts, ln in blocks:
        offsets.append((qs + 1, ts + 1, ts - qs))
    distinct = sorted({o for _, _, o in offsets})

    report = DERIVED / "version_crosswalk.md"
    with report.open("w") as fh:
        fh.write(f"# `{OLD}` against `{NEW}`\n\n")
        fh.write(f"Built by `src/build_version_crosswalk.py`. The two versions of the "
                 f"same accession differ by {len(old_seq) - len(new_seq)} bp in total "
                 f"length, and **the offset between them is not constant**.\n\n")
        fh.write(f"- Ungapped blocks: **{len(blocks)}**\n")
        fh.write(f"- Mapped: **{sum(b[2] for b in blocks):,}** of {len(old_seq):,} bp "
                 f"({100 * sum(b[2] for b in blocks) / len(old_seq):.2f} per cent)\n")
        fh.write(f"- Dropped as ambiguous or non-monotonic: "
                 f"{dropped_amb + dropped_mono:,} bp\n")
        fh.write(f"- **Distinct offsets: {len(distinct)}** ({', '.join(f'{o:+d}' for o in distinct)})\n\n")
        fh.write("A coordinate quoted against one version and read against the other is "
                 "wrong by whichever of these applies at that position, and nothing in "
                 "either record says which version a paper used unless the paper does.\n\n")
        fh.write("## Blocks\n\n")
        fh.write(f"| `{OLD}` from | to | `{NEW}` from | length | offset |\n")
        fh.write("|---|---|---|---|---|\n")
        for qs, ts, ln in blocks:
            fh.write(f"| {qs + 1:,} | {qs + ln:,} | {ts + 1:,} | {ln:,} | "
                     f"{ts - qs:+d} |\n")

    print(f"  {len(blocks)} blocks, {len(distinct)} distinct offsets: "
          f"{', '.join(f'{o:+d}' for o in distinct)}")
    print(f"  wrote {out.name}, {report.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
