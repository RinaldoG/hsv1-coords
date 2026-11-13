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


def python_blocks(body: str) -> list[str]:
    return re.findall(r"```python\n(.*?)```", body, re.S)


class TheReadmeExamplesRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocks = python_blocks(README.read_text())
        if not cls.blocks:
            raise unittest.SkipTest("the README states no Python examples")
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
        """Each `# ...` comment after a print must appear in the real output.

        This is what makes the examples a gate rather than decoration: an output
        comment that drifts from the behaviour fails here.
        """
        checked = 0
        for i, block in enumerate(self.blocks):
            if "print(" not in block:
                continue
            out = self._run(block)
            flat = re.sub(r"\s+", " ", out)
            lines = block.splitlines()
            for n, line in enumerate(lines):
                stripped = line.strip()
                # An output comment is a `#` line that follows a print or another
                # output comment, and is not prose about the code.
                if not stripped.startswith("# "):
                    continue
                # Directly after a print, or after another output comment, with no
                # blank line between. A prose comment separated by a blank line is a
                # note about the code, not quoted output.
                prior = lines[n - 1].strip() if n else ""
                if not (prior.startswith("print(") or prior.startswith("# ")):
                    continue
                fragment = re.sub(r"\s+", " ", stripped[2:]).strip()
                if len(fragment) < 8:
                    continue
                with self.subTest(block=i, quoted=fragment[:60]):
                    self.assertIn(
                        fragment, flat,
                        f"the README quotes output the code does not print:\n"
                        f"  quoted: {fragment}\n  actual: {flat[:600]}")
                checked += 1
        self.assertGreater(checked, 5,
                           f"only {checked} quoted output lines were checked, so this "
                           f"proves little")

    def test_an_inline_expectation_comment_matches_the_value(self):
        """`expr  # value` comments, such as `resolve("VP16")  # 'UL48'`."""
        checked = 0
        for i, block in enumerate(self.blocks):
            for line in block.splitlines():
                m = re.match(r"^([A-Za-z_][^#]*?)\s{2,}#\s*(\S.*)$", line)
                if not m or m.group(1).strip().startswith("print("):
                    continue
                expr, expected = m.group(1).strip(), m.group(2).strip()
                if expr.endswith(("(", ",")) or expected.startswith(("or ", "the ")):
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
