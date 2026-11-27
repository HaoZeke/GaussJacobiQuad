#!/usr/bin/env python3
"""Split a project version into major, minor, and patch.

The compact integer matches tblite: 10000 * major + 100 * minor + patch.
Each component must be in 0..99. A larger component collides
(1.100.0 == 2.0.0, 1.2.100 == 1.3.0).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_LEADING_INT = re.compile(r"^(\d+)")


class VersionError(ValueError):
    """The version cannot be packed into GJP_VERSION."""


def parse_version(version: str) -> tuple[int, int, int]:
    text = version.strip()
    if not text:
        raise VersionError("empty version")
    core = text.split("+", 1)[0]
    core = core.split("-", 1)[0]
    if not core:
        raise VersionError(f"no numeric version in {version!r}")
    nums: list[int] = []
    for index, part in enumerate(core.split(".")):
        match = _LEADING_INT.match(part)
        if not match:
            raise VersionError(
                f"component {part!r} in {version!r} has no leading integer"
            )
        rest = part[match.end() :]
        if rest and index != len(core.split(".")) - 1:
            raise VersionError(
                f"non-numeric {rest!r} in component {part!r} of {version!r}"
            )
        nums.append(int(match.group(1)))
        if rest:
            break
    if len(nums) > 3:
        raise VersionError(f"{version!r} has more than 3 numeric components")
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


def main(argv: list[str]) -> int:
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
        "       parse_gjp_version.py --check-fortran SRC.f90 VERSION\n"
    )
    return 2


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except VersionError as exc:
        sys.stderr.write(f"error: {exc}\n")
        sys.exit(1)
