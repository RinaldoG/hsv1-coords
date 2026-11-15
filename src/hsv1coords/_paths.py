"""Where this package's data lives, whether installed or run from a checkout.

Round 11: every module computed its data path as `parents[2] / "data"`, which is the
repo layout. In a non-editable install that resolves inside site-packages' parent, so
`pip install .` produced an importable library in which every loader raised
FileNotFoundError and `Crosswalk.available()` returned an empty list, making a sweep a
silent no-op. `pip install -e .` worked, which is what the README recommends, so the
breakage only hit the user who installed it as a dependency.

The data ships inside the package now (`hsv1coords/data`, wired in pyproject.toml), and
this prefers the checkout, falls back to the installed copy, and says so when neither is
there.

Round 12: `data_path` raised, and five modules call it at module scope, so a build from
a tree with no staged data made `import hsv1coords` itself fail. That is worse than what
round 11 fixed, not better: before, the library imported and failed per call; after, it
could not be imported at all, for anyone doing `pip install git+...` or installing an
sdist. Resolution is therefore split in two. `data_path` never raises, so a module-level
constant is safe. `require()` is what a loader calls when it is about to read, and that
is where the explanation lives.

`requested()` records every path the library has asked for, which is how
`src/stage_package_data.py` knows what a wheel has to carry: the staged list is derived
from what the loaders actually read rather than typed into a tuple beside them.
"""
from __future__ import annotations

import pathlib

_PACKAGE_DATA = pathlib.Path(__file__).resolve().parent / "data"
_REPO_DATA = pathlib.Path(__file__).resolve().parents[2] / "data"

_REQUESTED: set[str] = set()


def _explain() -> str:
    return (f"hsv1coords cannot find its data. Looked in {_REPO_DATA} (a repository "
            f"checkout) and {_PACKAGE_DATA} (an installed copy). If you installed from "
            f"a source tree without the data built, run `make all` in the checkout and "
            f"install again, or install with `pip install -e .`")


def data_dir(strict: bool = True) -> pathlib.Path:
    """The data directory: the checkout's if there is one, else the installed copy.

    The checkout wins deliberately. `make package-data` stages a copy inside the package
    for wheel building, and if that copy were preferred then development in a checkout
    would silently read a stale snapshot of its own artefacts. In an installed copy there
    is no checkout beside it, so the staged data is what gets used.

    `strict=False` returns the installed location even when it does not exist, for
    callers that only need to name a path. Nothing is read, so nothing can be wrong yet.
    """
    for candidate in (_REPO_DATA, _PACKAGE_DATA):
        if candidate.is_dir():
            return candidate
    if strict:
        raise FileNotFoundError(_explain())
    return _PACKAGE_DATA


def data_path(*parts: str) -> pathlib.Path:
    """Name a data file. Never raises: module-level constants are built with this."""
    _REQUESTED.add("/".join(parts))
    return data_dir(strict=False).joinpath(*parts)


def require(path: pathlib.Path) -> pathlib.Path:
    """Assert a data file is there, with the two locations named if it is not.

    Called by a loader at the moment it reads, so the diagnosis arrives at the call that
    needs the file rather than at import.
    """
    if path.exists():
        return path
    if not (_REPO_DATA.is_dir() or _PACKAGE_DATA.is_dir()):
        raise FileNotFoundError(_explain())
    raise FileNotFoundError(
        f"{path} is missing. It is a built artefact: run `make all` in a checkout, or "
        f"reinstall a wheel built with `make package-data` first.")


def requested() -> tuple[str, ...]:
    """Every data file the library has asked for in this process, as a relative path.

    Used by the staging script, so what a wheel carries is derived from what the loaders
    read. Round 12: the staged list was a tuple typed beside the copying code and the
    test over it was membership-only against five names, so dropping the annotation
    table or the block files from it passed.
    """
    return tuple(sorted(_REQUESTED))
