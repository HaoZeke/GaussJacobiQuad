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
"""The Fortran insert-license hook inserts one banner and then stays put.

Copies put the f2py mode line at index 0..5, with and without a trailing
newline. A banner that already starts past the old 5-line window is kept.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODE = "!    -*- f90 -*-"
BODY = "module keep_me\n  integer :: x = 1\nend module\n"


def _pre_commit():
    found = shutil.which("pre-commit")
    if found:
        return [found]
    return [sys.executable, "-m", "pre_commit"]


def _banner_lines():
    raw = (ROOT / "scripts" / "license_header.txt").read_text()
    prefixed = []
    for line in raw.splitlines(keepends=True):
        extra = " " if line.strip() else ""
        prefixed.append("!%s%s" % (extra, line))
    if prefixed and not prefixed[-1].endswith("\n"):
        prefixed[-1] += "\n"
    return "".join(prefixed).splitlines()


def _fresh(index, newline):
    parts = ["! pad %d\n" % i for i in range(index)]
    if newline:
        parts.append(MODE + "\n")
        parts.append(BODY)
    else:
        parts.append(MODE)
    return "".join(parts)


def _already(banner):
    parts = ["! pad %d\n" % i for i in range(4)]
    parts.append(MODE + "\n")
    parts.append(banner)
    parts.append(BODY)
    return "".join(parts)


def _run(paths):
    return subprocess.run(
        _pre_commit() + ["run", "insert-license", "--files"] + paths,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


def _one_banner(text):
    return sum(1 for line in text.splitlines() if line.strip() == "! BEGIN_HEADER")


def _check_fresh(index, newline, text, banner_lines):
    lines = text.splitlines()
    if _one_banner(text) != 1:
        return "banner count %s" % _one_banner(text)
    if (MODE + "!") in text:
        return "mode line glued to the banner"
    width = len(banner_lines)
    if index <= 1:
        if lines[index] != MODE:
            return "mode line moved to %r" % (lines[index] if lines else "")
        if lines[index + 1 : index + 1 + width] != banner_lines:
            return "banner is not directly under the mode line"
    else:
        if lines[:width] != banner_lines:
            return "banner is not at the start of the file"
        mode_at = width + index
        if mode_at >= len(lines) or lines[mode_at] != MODE:
            return "mode line was not kept after the leading banner"
    if newline:
        for line in ("module keep_me", "  integer :: x = 1", "end module"):
            if line not in lines:
                return "lost %r" % line
    else:
        for pad in range(index):
            if "! pad %d" % pad not in lines:
                return "lost pad %d" % pad
        if index >= 2 and text.endswith("\n"):
            return "unterminated mode line gained a newline"
    return None


def main():
    banner_lines = _banner_lines()
    banner = "\n".join(banner_lines) + "\n"
    failed = 0
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        folder = Path(tmp)
        specs = []
        for index in range(6):
            for newline in (True, False):
                for suffix in ("f90", "pyf"):
                    name = "mode_%d_%s.%s" % (
                        index,
                        "nl" if newline else "nonew",
                        suffix,
                    )
                    path = folder / name
                    path.write_text(_fresh(index, newline))
                    specs.append((path, index, newline, None))
        for suffix in ("f90", "pyf"):
            path = folder / ("already_deep.%s" % suffix)
            original = _already(banner)
            path.write_text(original)
            specs.append((path, 4, True, original))

        rels = [str(path.relative_to(ROOT)) for path, _i, _n, _o in specs]
        snapshots = []
        for run in (1, 2, 3):
            result = _run(rels)
            if run == 1 and result.returncode not in (0, 1):
                sys.stderr.write(result.stderr)
                print(
                    "FAIL insert-license run 1 EXIT:%s" % result.returncode,
                    file=sys.stderr,
                )
                return 1
            if run > 1 and result.returncode != 0:
                sys.stderr.write(result.stdout)
                sys.stderr.write(result.stderr)
                print(
                    "FAIL insert-license run %s EXIT:%s" % (run, result.returncode),
                    file=sys.stderr,
                )
                failed += 1
            current = []
            for path, index, newline, original in specs:
                text = path.read_text()
                current.append(text)
                if original is not None:
                    if text != original or _one_banner(text) != 1:
                        print(
                            "FAIL %s: deep banner changed" % path.name,
                            file=sys.stderr,
                        )
                        failed += 1
                    continue
                if run == 1:
                    reason = _check_fresh(index, newline, text, banner_lines)
                    if reason:
                        print("FAIL %s: %s" % (path.name, reason), file=sys.stderr)
                        failed += 1
            if run > 1 and current != snapshots:
                print(
                    "FAIL insert-license run %s rewrote a file" % run, file=sys.stderr
                )
                failed += 1
            snapshots = current
    if failed:
        print("failed %s" % failed)
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
