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
"""Insert the MIT banner into Fortran sources for the pre-commit hook.

The f2py mode line counts only on the first two lines. An existing banner
counts wherever it already sits, so a later run does not stack another one.
A mode line with no trailing newline is broken before the banner is written.
"""

import os
import re
import sys

_MODE_LINE = re.compile(r"^!\s*-\*-\s*f90\s*-\*-$")


def _prefix(path):
    with open(path, encoding="utf-8", newline="") as handle:
        plain = handle.readlines()
    prefixed = []
    for line in plain:
        extra = " " if line.strip() else ""
        prefixed.append("!%s%s" % (extra, line))
    if prefixed and not prefixed[-1].endswith("\n"):
        prefixed[-1] += "\n"
    return prefixed


def _read(path):
    last_error = None
    for encoding in ("utf-8", "latin-1"):
        try:
            with open(path, "r", encoding=encoding, newline="") as handle:
                return handle.readlines(), encoding
        except UnicodeDecodeError as exc:
            last_error = exc
    print("Error while processing: %s" % path, file=sys.stderr)
    if last_error is not None:
        raise last_error
    raise RuntimeError("could not read %s" % path)


def _has_banner(lines, banners):
    stripped = [line.strip() for line in lines]
    for banner in banners:
        want = [line.strip() for line in banner]
        width = len(want)
        if width == 0 or width > len(stripped):
            continue
        for start in range(len(stripped) - width + 1):
            if stripped[start : start + width] == want:
                return True
    return False


def _insertion_index(lines):
    """Mode line only on line 1 or line 2; otherwise the banner goes first."""
    for index in range(min(2, len(lines))):
        if _MODE_LINE.match(lines[index].strip()):
            return index + 1
    return 0


def _line_ending(lines):
    for line in lines:
        if line.endswith("\r\n"):
            return "\r\n"
        if line.endswith("\n"):
            return "\n"
    return "\n"


def _with_break(lines, index):
    if index <= 0:
        return lines
    previous = lines[index - 1]
    if previous.endswith("\n"):
        return lines
    updated = list(lines)
    updated[index - 1] = previous + _line_ending(lines)
    return updated


def _process(path, banners):
    lines, encoding = _read(path)
    if _has_banner(lines, banners):
        return False
    index = _insertion_index(lines)
    lines = _with_break(lines, index)
    updated = lines[:index] + banners[0] + lines[index:]
    with open(path, "w", encoding=encoding, newline="") as handle:
        handle.write("".join(updated))
    return True


def main(argv):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    scripts = os.path.join(root, "scripts")
    banners = [
        _prefix(os.path.join(scripts, name))
        for name in ("license_header.txt", "license_header_legacy.txt")
    ]
    changed = False
    for path in argv[1:]:
        if _process(path, banners):
            changed = True
    return 1 if changed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
