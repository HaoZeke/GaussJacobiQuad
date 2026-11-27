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
"""A headerless gjquadpy.pyf keeps the f2py mode line on line 1."""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYF = ROOT / "interfaces" / "PyInterface" / "gjquadpy.pyf"
COPY = ROOT / "interfaces" / "PyInterface" / "gjquadpy_mode_check.pyf"


def _strip_banner(text):
    lines = text.splitlines(keepends=True)
    begin = next(i for i, line in enumerate(lines) if line.strip() == "! BEGIN_HEADER")
    end = next(
        i
        for i, line in enumerate(lines)
        if i > begin and line.strip() == "! END_HEADER"
    )
    return "".join(lines[:begin] + lines[end + 1 :])


def _pre_commit():
    found = shutil.which("pre-commit")
    if found:
        return [found]
    return [sys.executable, "-m", "pre_commit"]


def main():
    original = PYF.read_text()
    mode = original.splitlines()[0]
    COPY.write_text(_strip_banner(original))
    try:
        result = subprocess.run(
            _pre_commit()
            + ["run", "insert-license", "--files", str(COPY.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        lines = COPY.read_text().splitlines()
    finally:
        COPY.unlink(missing_ok=True)
    if not lines or lines[0] != mode:
        print("FAIL mode line: %r" % (lines[0] if lines else ""), file=sys.stderr)
        return 1
    if len(lines) < 2 or lines[1].strip() != "! BEGIN_HEADER":
        print("FAIL banner was not inserted under the mode line", file=sys.stderr)
        return 1
    if sum(1 for line in lines if line.strip() == "! BEGIN_HEADER") != 1:
        print("FAIL banner count", file=sys.stderr)
        return 1
    if result.returncode not in (0, 1):
        sys.stderr.write(result.stderr)
        print("FAIL insert-license EXIT:%s" % result.returncode, file=sys.stderr)
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
