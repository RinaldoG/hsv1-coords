"""Run the README's own examples, because round 8 found they did not run at all.

Three faults in five lines of quickstart, after seven rounds of review, none of which
tried the first line the way a new user would:

- `from hsv1coords import ...` raised `ModuleNotFoundError` on a fresh clone. There was
  no `pyproject.toml` and no `PYTHONPATH` note anywhere in the file.
- The flagship `lift` call was annotated with a number; it returns `None`, because the
  position is inside a repeat and the lift declines, which the same file says it does.
- One number, 118777, was used as a coordinate in two different genomes four lines
  apart, so the counterpart printed beside it belonged to the other frame.

Every reviewer read the examples. None executed them. So the suite does, and it
executes them as a subprocess with only the package path set, which is the state a
reader is in.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
DOCUMENTS = ("README.md", "ROADMAP.md", "FINDINGS.md", "CHANGELOG.md")
"""Every document, not only the README.

Round 10: only README python blocks were executed, so ROADMAP.md still showed two
flagship examples handing out an alternative with no strand note, nine commits after the
strand note was added and the README was fixed. That is PITFALLS root cause 14, fixing
the file the blocker quoted.
"""


def python_blocks(body: str) -> list[str]:
    return re.findall(r"```python\n(.*?)```", body, re.S)


def output_groups(block: str) -> list[str]:
    """The quoted-output comment groups of a block, in order.

    A group is a run of `#` comment lines. It is quoted output when the line before the
    run is a `print(`, and prose otherwise. Round 9: the earlier rule asked whether the
    PREVIOUS line was a print or another comment, so the first line of a multi-line
    prose comment was correctly skipped and every line after it was read as output. A
    three-line note about the example was compared against the program's stdout.
    """
    groups: list[str] = []
    lines = block.splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip().startswith("#"):
            i += 1
            continue
        start = i
        run = []
        while i < len(lines) and lines[i].strip().startswith("#"):
            run.append(lines[i].strip().lstrip("#").strip())
            i += 1
        prior = lines[start - 1].strip() if start else ""
        if prior.startswith("print("):
            groups.append(" ".join(x for x in run if x))
    return groups


class TheReadmeExamplesRun(unittest.TestCase):
    source_of: dict[str, str] = {}

    @classmethod
    def setUpClass(cls):
        cls.blocks = []
        for name in DOCUMENTS:
            path = ROOT / name
            if not path.exists():
                continue
            for block in python_blocks(path.read_text()):
                if not re.search(r"^(from|import) ", block, re.M):
                    # A fragment with no import cannot be executed. It is also not
                    # allowed to quote output: that is how an example dodges this gate.
                    continue
                cls.blocks.append(block)
                cls.source_of[block] = name
        if not cls.blocks:
            raise unittest.SkipTest("the documents state no Python examples")
        for need in ("data/derived/crosswalk", "data/derived/repeat_map.tsv"):
            if not (ROOT / need).exists():
                raise unittest.SkipTest(f"{need} is built, not committed; run `make all`")

    def _run(self, code: str) -> str:
        # Not `-I`: isolated mode ignores PYTHONPATH, so this gate could not import
        # the package it exists to prove a reader can import. `-E` would do the same.
        proc = subprocess.run(
            [sys.executable, "-c", code], cwd=ROOT, capture_output=True,
            text=True, env={"PYTHONPATH": str(ROOT / "src"), "PATH": "/usr/bin:/bin",
                            "PYTHONDONTWRITEBYTECODE": "1"})
        self.assertEqual(
            proc.returncode, 0,
            f"a README example failed to run:\n{proc.stderr[-2000:]}")
        return proc.stdout

    def test_every_python_block_executes(self):
        """A fresh reader runs this first. It must not raise."""
        for i, block in enumerate(self.blocks):
            with self.subTest(block=i):
                self._run(block)

    def test_the_quoted_outputs_are_what_the_code_prints(self):
        """Each quoted output line must appear in the real output.

        Kept beside the equality check below: this one names the offending fragment,
        which is the more useful failure when a single line has drifted.
        """
        checked = 0
        for i, block in enumerate(self.blocks):
            if "print(" not in block:
                continue
            flat = re.sub(r"\s+", " ", self._run(block))
            for group in output_groups(block):
                fragment = re.sub(r"\s+", " ", group).strip()
                if len(fragment) < 8:
                    continue
                with self.subTest(block=i, quoted=fragment[:60]):
                    self.assertIn(
                        fragment, flat,
                        f"the README quotes output the code does not print:\n"
                        f"  quoted: {fragment}\n  actual: {flat[:600]}")
                checked += 1
        self.assertGreater(checked, 3,
                           f"only {checked} quoted output groups were checked, so this "
                           f"proves little")

    def test_the_quoted_output_is_the_whole_output_in_order(self):
        """Equality, not containment, because containment missed three things.

        The check above asserts each quoted fragment appears somewhere in the real
        output. Round 9 showed what that cannot see: the README quoted six of the seven
        lines `report` printed and dropped the seventh, which was the ambiguity caveat,
        the most load-bearing line in the file; reordering three quoted lines passed;
        and deleting two quoted lines entirely passed, because a `findall` over zero
        occurrences runs zero subtests and reports success. All three were under a
        heading that says "Every output below is copied from a real run".

        So this reassembles every quoted output comment in the block, in order, and
        requires it to equal what the block actually prints. An example may still be
        written without quoting its output; it may not quote part of it and leave the
        rest unsaid.
        """
        checked, unquoted = 0, []
        for i, block in enumerate(self.blocks):
            if "print(" not in block:
                continue
            groups = output_groups(block)
            if not groups:
                # A block with a print and no quoted output: run it, so a claim in its
                # prose cannot rest on code that does not work, but there is nothing to
                # compare. Round 10 added such a block whose comment said every position
                # in the inverted repeats lifts cleanly, and nothing looked at it.
                self._run(block)
                unquoted.append(self.source_of.get(block, "?"))
                continue
            expected = re.sub(r"\s+", " ", " ".join(groups)).strip()
            actual = re.sub(r"\s+", " ", self._run(block)).strip()
            with self.subTest(block=i):
                self.assertEqual(
                    actual, expected,
                    "the README's quoted output is not what the block prints, in "
                    "order and in full. A line the code prints and the README does "
                    "not quote is the failure this catches.")
            checked += 1
        self.assertGreater(checked, 0,
                           "no block with quoted output was compared, so this gate ran "
                           "on nothing")
        self.assertEqual(
            unquoted, [],
            f"{len(unquoted)} example blocks print something and quote none of it, so "
            f"whatever they claim in prose is unchecked: {unquoted}")

    def test_an_inline_expectation_comment_matches_the_value(self):
        """`expr  # value` comments, such as `resolve("VP16")  # 'UL48'`."""
        checked = 0
        for i, block in enumerate(self.blocks):
            for line in block.splitlines():
                m = re.match(r"^([A-Za-z_][^#]*?)\s{2,}#\s*(\S.*)$", line)
                if not m:
                    continue
                if m.group(1).strip().startswith("print("):
                    # A trailing comment on a print line was read by nothing: the
                    # inline check skipped it and `output_groups` only collects
                    # whole-line comments. Round 10 turned "# the LAT TATA box" into a
                    # false claim about the position with the suite green. Such a
                    # comment is a label, so it must not look like an assertion about a
                    # value: no digits, and no claim words.
                    label = m.group(2).strip()
                    bad = re.search(r"\d|\bis\b|\bequals\b|\breturns\b|"
                                    r"\bwell defined\b|\bnone\b", label, re.I)
                    with self.subTest(label=label[:50]):
                        self.assertIsNone(
                            bad,
                            f"a trailing comment on a print line is a label, and "
                            f"nothing checks it, so it must not read as a statement "
                            f"about the output: {label!r}")
                    continue
                expr, expected = m.group(1).strip(), m.group(2).strip()
                # Round 10: skipping any expectation starting "or " or "the " let a
                # wrong value hide behind one word. `# 59722` was changed to `# the
                # answer is 99999` and the suite stayed green. An expectation that is
                # prose rather than a value has to say so with a marker, so that the
                # skip is a decision somebody wrote down rather than a word they
                # happened to start with.
                if expr.endswith(("(", ",")):
                    continue
                if expected.startswith("prose:"):
                    continue
                header = "\n".join(l for l in block.splitlines()
                                   if l.startswith(("from ", "import ")))
                out = self._run(f"{header}\nprint(repr({expr}))")
                with self.subTest(block=i, expression=expr):
                    self.assertEqual(
                        out.strip(), expected,
                        f"the README says `{expr}` gives {expected}, and it gives "
                        f"{out.strip()}")
                checked += 1
        self.assertGreater(checked, 2,
                           f"only {checked} inline expectations were checked")

    def test_the_readme_says_how_to_install_the_package(self):
        """Round 8's blocker: the examples import a package nothing installs."""
        body = README.read_text()
        self.assertTrue(
            re.search(r"pip install|PYTHONPATH", body),
            "the README's examples import `hsv1coords` but the file never says how to "
            "make that import work. It fails on a fresh clone.")
        self.assertTrue((ROOT / "pyproject.toml").exists(),
                        "there is no pyproject.toml, so `pip install -e .` cannot work")


class EveryQuotedExampleIsExecutable(unittest.TestCase):
    """A block that quotes output must be runnable, or the quote is unchecked.

    Round 10: only README blocks were executed, and two ROADMAP blocks quoted output
    while importing nothing, so their text could say anything. One of them did: it
    showed an alternative with no strand note, nine commits after the strand note was
    added.
    """

    def test_a_block_that_quotes_output_imports_what_it_uses(self):
        offenders = []
        for name in DOCUMENTS:
            path = ROOT / name
            if not path.exists():
                continue
            for i, block in enumerate(python_blocks(path.read_text())):
                if re.search(r"^(from|import) ", block, re.M):
                    continue
                lines = block.splitlines()
                quotes = [n for n, line in enumerate(lines)
                          if line.strip().startswith("#") and n
                          and not lines[n - 1].strip().startswith("#")
                          and lines[n - 1].strip()
                          and not lines[n - 1].strip().startswith("#")]
                if quotes:
                    offenders.append((name, i, lines[quotes[0]].strip()[:60]))
        self.assertEqual(
            offenders, [],
            f"{len(offenders)} example blocks quote output without importing anything, "
            f"so nothing can run them and the quoted text is unchecked: {offenders}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
