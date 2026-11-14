"""Where this package's data lives, whether installed or run from a checkout.

Round 11: every module computed its data path as `parents[2] / "data"`, which is the
repo layout. In a non-editable install that resolves inside site-packages' parent, so
`pip install .` produced an importable library in which every loader raised
FileNotFoundError and `Crosswalk.available()` returned an empty list, making a sweep a
silent no-op. `pip install -e .` worked, which is what the README recommends, so the
breakage only hit the user who installed it as a dependency.

The data ships inside the package now (`hsv1coords/data`, wired in pyproject.toml), and
this prefers the installed copy, falls back to the checkout, and raises with both paths
named if neither is there.
"""
from __future__ import annotations

import pathlib

_PACKAGE_DATA = pathlib.Path(__file__).resolve().parent / "data"
_REPO_DATA = pathlib.Path(__file__).resolve().parents[2] / "data"


def data_dir() -> pathlib.Path:
    """The data directory: the checkout's if there is one, else the installed copy.

    The checkout wins deliberately. `make package-data` stages a copy inside the package
    for wheel building, and if that copy were preferred then development in a checkout
    would silently read a stale snapshot of its own artefacts. In an installed copy there
    is no checkout beside it, so the staged data is what gets used.
    """
    for candidate in (_REPO_DATA, _PACKAGE_DATA):
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        f"hsv1coords cannot find its data. Looked in {_REPO_DATA} (a repository "
        f"checkout) and {_PACKAGE_DATA} (an installed copy). If you installed from a "
        f"source tree without the data built, run `make all` in the checkout and "
        f"install again, or install with `pip install -e .`")


def data_path(*parts: str) -> pathlib.Path:
    return data_dir().joinpath(*parts)
