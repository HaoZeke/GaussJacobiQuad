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
"""The header updater must leave the source line after END_HEADER intact."""

import importlib.util
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = "Authors: Rohit Goswami <rgoswami[at]ieee.org>"
SAMPLES = (
    (ROOT / "src" / "gjp_imtqlx.f90", "!", "f90"),
    (ROOT / "interfaces" / "CInterface" / "c_cli_gjpq.c", "//", "c"),
    (ROOT / "pytests" / "test_gjp_rec.py", "#", "py"),
)


def load_add_headers():
    path = ROOT / "scripts" / "add_headers.py"
    spec = importlib.util.spec_from_file_location("add_headers_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def suffix_after(text, cchar):
    marker = f"{cchar} END_HEADER\n"
    idx = text.find(marker)
    if idx == -1:
        return None
    return text[idx + len(marker) :]


def first_line(text):
    if not text:
        return None
    return text.splitlines()[0]


def check_syntax(dest, suffix, src):
    if suffix == "c":
        header = src.parent / "GaussJacobiQuadCInterp.h"
        shutil.copyfile(header, dest.parent / header.name)
        cmd = ["gcc", "-fsyntax-only", str(dest)]
    elif suffix == "py":
        cmd = [sys.executable, "-m", "py_compile", str(dest)]
    else:
        return 0
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        print("FAIL %s EXIT:%s" % (" ".join(cmd), result.returncode), file=sys.stderr)
    return result.returncode


def main():
    updater = load_add_headers()
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for src, cchar, suffix in SAMPLES:
            original = src.read_text()
            expected = suffix_after(original, cchar)
            if expected is None:
                print("FAIL %s: sample has no END_HEADER line" % src.name, file=sys.stderr)
                failed += 1
                continue
            dest_dir = tmp / suffix
            dest_dir.mkdir()
            dest = dest_dir / src.name
            dest.write_text(original.replace(AUTHOR, "Authors: SENTINEL", 1))
            updater.add_headers([str(dest_dir)], [suffix], cchar=cchar)
            updated = dest.read_text()
            got = suffix_after(updated, cchar)
            if "SENTINEL" in updated or got != expected:
                print(
                    "FAIL %s: next line %r != %r"
                    % (src.name, first_line(got), first_line(expected)),
                    file=sys.stderr,
                )
                failed += 1
                continue
            updater.add_headers([str(dest_dir)], [suffix], cchar=cchar)
            if dest.read_text() != updated:
                print("FAIL %s: second run changed the file" % src.name, file=sys.stderr)
                failed += 1
                continue
            if check_syntax(dest, suffix, src) != 0:
                failed += 1
    if failed:
        print("failed %s" % failed)
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
