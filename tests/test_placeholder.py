"""M0 holds no logic to test yet. This asserts the repo layout the gate expects."""
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]


class Layout(unittest.TestCase):
    def test_expected_paths_exist(self):
        for rel in ("README.md", "ROADMAP.md", "LICENSE", "Makefile",
                    "src/fetch_inputs.py", "src/check_inputs.py"):
            self.assertTrue((ROOT / rel).exists(), f"missing {rel}")
