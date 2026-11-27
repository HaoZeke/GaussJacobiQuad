#!/usr/bin/env python3
# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
"""Split a project version into major, minor, and patch.

The compact integer matches tblite: 10000 * major + 100 * minor + patch.
Each component must be in 0..99. A larger component collides
(1.100.0 == 2.0.0, 1.2.100 == 1.3.0).

A hyphen or plus suffix is SemVer pre-release or build metadata. A dotted
or glued dev/post/rc tail (1.2.0.dev1, 1.2.0rc.1) is the same release triple.
A fourth numeric component and junk in the last component are rejected.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Longest tags first so "preview" is not parsed as "pre".
_TAIL = re.compile(
    r"(?i)(?:[._]?(?:alpha|beta|preview|pre|post|rev|rc|dev|a|b|c|r)[._]?\d*)+$"
)
_META = re.compile(r"^[0-9A-Za-z.-]+$")
_ROOT = Path(__file__).resolve().parents[1]


class VersionError(ValueError):
    """The version cannot be packed into GJP_VERSION."""


def parse_version(version: str) -> tuple[int, int, int]:
    text = version.strip()
    if not text:
        raise VersionError("empty version")
    core, _meta = _split_meta(text)
    parts = core.split(".")
    if not 1 <= len(parts) <= 3:
        raise VersionError(f"{version!r} has more than 3 numeric components")
    nums: list[int] = []
    for part in parts:
        if not part.isdigit():
            raise VersionError(f"component {part!r} in {version!r} is not an integer")
        nums.append(int(part))
    while len(nums) < 3:
        nums.append(0)
    major, minor, patch = nums
    for name, value in (("major", major), ("minor", minor), ("patch", patch)):
        if value > 99:
            raise VersionError(
                f"{name} {value} in {version!r} exceeds 99; "
                "GJP_VERSION = 10000*major + 100*minor + patch collides past 99"
            )
    return major, minor, patch


def _split_meta(text: str) -> tuple[str, str]:
    core = text
    meta = ""
    if "+" in core:
        core, meta = core.split("+", 1)
        if not meta or _META.fullmatch(meta) is None:
            raise VersionError(f"bad build metadata in {text!r}")
    if "-" in core:
        core, pre = core.split("-", 1)
        if not pre or _META.fullmatch(pre) is None:
            raise VersionError(f"bad pre-release in {text!r}")
    else:
        stripped = _TAIL.sub("", core)
        if stripped != core:
            if not stripped or stripped.endswith("."):
                raise VersionError(f"bad pre-release tail in {text!r}")
            core = stripped
    if not core:
        raise VersionError(f"no numeric version in {text!r}")
    return core, meta


def fortran_constants(path: Path) -> tuple[int, int, int, str]:
    text = path.read_text()

    def grab(name: str) -> int:
        match = re.search(rf"gjp_version_{name}\s*=\s*(\d+)", text)
        if not match:
            raise VersionError(f"missing gjp_version_{name} in {path}")
        return int(match.group(1))

    string_match = re.search(r'gjp_version_string\s*=\s*"([^"]*)"', text)
    if not string_match:
        raise VersionError(f"missing gjp_version_string in {path}")
    return grab("major"), grab("minor"), grab("patch"), string_match.group(1)


def _replace_one(text: str, pattern: str, repl: str, path: Path) -> str:
    new, count = re.subn(pattern, lambda _match: repl, text, count=1, flags=re.M)
    if count != 1:
        raise VersionError(f"{path}: pattern did not match once")
    return new


def apply_version(version: str, root: Path | None = None) -> tuple[int, int, int]:
    """Rewrite release files. fpm.toml gets only major.minor.patch."""
    major, minor, patch = parse_version(version)
    triple = f"{major}.{minor}.{patch}"
    root = _ROOT if root is None else root
    fpm = root / "fpm.toml"
    meson = root / "meson.build"
    doxy = root / "apidocs" / "Doxygen-GaussJacobiQuad.cfg"
    citation = root / "CITATION.cff"
    source = root / "src" / "gjp_version.f90"
    updated = {
        fpm: _replace_one(
            fpm.read_text(),
            r'^version = "[^"]*"',
            f'version = "{triple}"',
            fpm,
        ),
        meson: _replace_one(
            meson.read_text(),
            r"version: '[^']*'",
            f"version: '{version}'",
            meson,
        ),
        doxy: _replace_one(
            doxy.read_text(),
            r"^PROJECT_NUMBER[ \t]*=.*",
            f'PROJECT_NUMBER          = "{version}"',
            doxy,
        ),
        citation: _replace_one(
            citation.read_text(),
            r"^version: .*",
            f"version: v{version}",
            citation,
        ),
        source: source.read_text(),
    }
    fortran = updated[source]
    for name, value in (("major", major), ("minor", minor), ("patch", patch)):
        fortran = _replace_one(
            fortran,
            rf"^integer, parameter :: gjp_version_{name} = .*",
            f"integer, parameter :: gjp_version_{name} = {value}",
            source,
        )
    fortran = _replace_one(
        fortran,
        r'^character\(len=\*\), parameter :: gjp_version_string = ".*"',
        f'character(len=*), parameter :: gjp_version_string = "{version}"',
        source,
    )
    updated[source] = fortran
    for path, text in updated.items():
        path.write_text(text)
    return major, minor, patch


def main(argv: list[str]) -> int:
    if len(argv) == 3 and argv[1] == "--apply":
        major, minor, patch = apply_version(argv[2])
        print(f"{major} {minor} {patch}")
        return 0
    if len(argv) == 4 and argv[1] == "--check-fortran":
        major, minor, patch = parse_version(argv[3])
        got_major, got_minor, got_patch, got_string = fortran_constants(Path(argv[2]))
        if (got_major, got_minor, got_patch) != (
            major,
            minor,
            patch,
        ) or got_string != argv[3]:
            raise VersionError(
                f"{argv[2]} has {got_major}.{got_minor}.{got_patch} ({got_string!r}) "
                f"but the project version is {argv[3]!r} -> {major}.{minor}.{patch}"
            )
        print(f"{major} {minor} {patch}")
        return 0
    if len(argv) == 2:
        major, minor, patch = parse_version(argv[1])
        print(f"{major} {minor} {patch}")
        return 0
    sys.stderr.write(
        "usage: parse_gjp_version.py VERSION\n"
        "       parse_gjp_version.py --apply VERSION\n"
        "       parse_gjp_version.py --check-fortran SRC.f90 VERSION\n"
    )
    return 2


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except VersionError as exc:
        sys.stderr.write(f"error: {exc}\n")
        sys.exit(1)
