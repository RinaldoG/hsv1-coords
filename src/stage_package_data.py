"""Copy the built artefacts into the package, so a wheel carries its own data.

Round 11: `pip install .` produced an importable library in which every loader raised
FileNotFoundError, because each module resolved its data as `parents[2] / "data"`, which
is the repository layout and not the installed one. `pip install -e .` worked, which is
what README.md recommends, so the breakage hit only the user who installed this as a
dependency, which is the user least able to diagnose it.

`hsv1coords._paths` prefers the checkout and falls back to `hsv1coords/data`, so this
script is what makes the second of those exist. The staged copy is gitignored: it is a
build product, not a second copy of the data to keep in step.

Round 12 found three things wrong with the round 11 version.

The staged list was a tuple typed beside the copying code, and the test over it was
membership-only against five names, so dropping the annotation table or the block files
from it passed. The list is now derived: the loaders are exercised, `_paths.requested()`
reports every data file they actually read, and that is what gets staged. A new loader
reading a new file is staged without anyone remembering to add it.

`accessions.txt` was staged and then not shipped, because no `pyproject` glob matched
`.txt`, and this script printed success. Every staged file is now checked against the
`package-data` globs, and a file that would not ship fails the run.

And nothing compared the staged copy with its source. It was stale on the machine this
was found on: the 15 alias rows round 11 corrected were still asserted as attested in
`src/hsv1coords/data`, invisible because `_paths` prefers the checkout and the staged
tree is gitignored, so neither the suite nor CI ever read it. Every file is compared
byte for byte after copying.
"""
from __future__ import annotations

import fnmatch
import pathlib
import shutil
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data"
TARGET = ROOT / "src" / "hsv1coords" / "data"

sys.path.insert(0, str(ROOT / "src"))


def what_the_loaders_read() -> list[str]:
    """Relative data paths, derived by running every loader this library offers."""
    import hsv1coords
    from hsv1coords import _paths

    hsv1coords.Dictionary.load()
    hsv1coords.ConventionSet.load()
    hsv1coords.RepeatMap.load()
    hsv1coords.VersionCrosswalk.load()
    available = hsv1coords.Crosswalk.available()
    if available:
        hsv1coords.Crosswalk.load(available[0])
    return list(_paths.requested())


def shipped_globs() -> list[str]:
    with (ROOT / "pyproject.toml").open("rb") as fh:
        config = tomllib.load(fh)
    return config["tool"]["setuptools"]["package-data"]["hsv1coords"]


def main() -> int:
    if not SOURCE.is_dir():
        print(f"no {SOURCE}; run the build first", file=sys.stderr)
        return 1

    wanted = what_the_loaders_read()
    if TARGET.exists():
        shutil.rmtree(TARGET)

    staged: list[str] = []
    missing: list[str] = []
    for relative in wanted:
        src = SOURCE.joinpath(*relative.split("/"))
        if src.is_dir():
            # A loader that names a directory (the crosswalk's 183 block files) rather
            # than a file. Every file directly inside it is staged, so this does not
            # depend on guessing the naming pattern.
            files = sorted(p for p in src.iterdir() if p.is_file())
            if not files:
                missing.append(f"{relative}/ (the directory is empty)")
                continue
            for one in files:
                dest = TARGET.joinpath(*relative.split("/"), one.name)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(one, dest)
                staged.append(f"{relative}/{one.name}")
        elif src.exists():
            dest = TARGET.joinpath(*relative.split("/"))
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            staged.append(relative)
        else:
            missing.append(relative)

    # Compared, not assumed. A staged copy that differs from its source ships a
    # different library from the one the suite tests.
    differs = [r for r in staged
               if TARGET.joinpath(*r.split("/")).read_bytes()
               != SOURCE.joinpath(*r.split("/")).read_bytes()]

    # And a staged file no `package-data` glob matches is staged and then not shipped.
    globs = shipped_globs()
    unshipped = [r for r in staged
                 if not any(fnmatch.fnmatch(f"data/{r}", g) for g in globs)]

    print(f"staged {len(staged)} files into {TARGET.relative_to(ROOT)}, derived from "
          f"{len(wanted)} paths the loaders read")
    problems = []
    if missing:
        problems.append(f"{len(missing)} files the loaders read are not built: "
                        f"{missing}. Run `make all` first, or a wheel built from this "
                        f"tree will fail on the first call.")
    if differs:
        problems.append(f"{len(differs)} staged files differ from their source: "
                        f"{differs[:5]}")
    if unshipped:
        problems.append(f"{len(unshipped)} staged files match no package-data glob in "
                        f"pyproject.toml, so they are copied here and then left out of "
                        f"the wheel: {unshipped[:5]}. Globs: {globs}")
    if problems:
        print("\n" + "\n".join(problems), file=sys.stderr)
        return 1
    print("every staged file is identical to its source and matched by a shipped glob")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
