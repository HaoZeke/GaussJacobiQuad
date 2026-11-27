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
import os
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


def check_crlf_and_shebang(updater, tmp):
    """A CRLF file stays CRLF, and a shebang stays above the banner."""
    failed = 0
    crlf_dir = tmp / "crlf"
    crlf_dir.mkdir()
    sample = crlf_dir / "sample.py"
    sample.write_bytes(b"# BEGIN_HEADER\r\n# old\r\n# END_HEADER\r\nimport stay\r\n")
    updater.add_headers([str(crlf_dir)], ["py"], cchar="#")
    data = sample.read_bytes()
    bare = data.replace(b"\r\n", b"")
    if b"\n" in bare or not data.endswith(b"import stay\r\n"):
        print("FAIL crlf: line endings or body changed", file=sys.stderr)
        failed += 1
    elif b"# END_HEADER\r\nimport stay\r\n" not in data:
        print("FAIL crlf: source was glued or dropped", file=sys.stderr)
        failed += 1
    else:
        updater.add_headers([str(crlf_dir)], ["py"], cchar="#")
        if sample.read_bytes() != data:
            print("FAIL crlf: second run changed the file", file=sys.stderr)
            failed += 1

    sh_dir = tmp / "shebang"
    sh_dir.mkdir()
    tool = sh_dir / "tool.py"
    tool.write_bytes(b"#!/usr/bin/env python3\r\nprint('hi')\r\n")
    updater.add_headers([str(sh_dir)], ["py"], cchar="#")
    sh_data = tool.read_bytes()
    prefix = b"#!/usr/bin/env python3\r\n# BEGIN_HEADER\r\n"
    if not sh_data.startswith(prefix):
        print("FAIL shebang: banner was placed above #!", file=sys.stderr)
        failed += 1
    elif b"print('hi')\r\n" not in sh_data or b"\n" in sh_data.replace(b"\r\n", b""):
        print("FAIL shebang: body or CRLF lost", file=sys.stderr)
        failed += 1
    else:
        updater.add_headers([str(sh_dir)], ["py"], cchar="#")
        if tool.read_bytes() != sh_data:
            print("FAIL shebang: second run changed the file", file=sys.stderr)
            failed += 1
    return failed


def _git_env():
    env = os.environ.copy()
    for key in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"):
        env.pop(key, None)
    return env


def check_shell_updater():
    """addHeaders.sh keeps a glued closer and skips a missing one."""
    script = ROOT / "scripts" / "addHeaders.sh"
    failed = 0
    env = _git_env()
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp)
        src = repo / "src"
        src.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=repo, check=True, env=env)
        glued = src / "glued.f90"
        glued.write_text(
            "! BEGIN_HEADER\n"
            "! old\n"
            "! END_HEADERmodule glued_should_survive\n"
            "SECOND_LINE\n"
        )
        missing = src / "partial.f90"
        original_missing = "! BEGIN_HEADER\nmodule survive\n"
        missing.write_text(original_missing)
        result = subprocess.run(
            ["bash", str(script)],
            cwd=repo,
            capture_output=True,
            text=True,
            env=env,
        )
        glued_lines = glued.read_text().splitlines()
        if "module glued_should_survive" not in glued_lines:
            print("FAIL addHeaders.sh glued: module line removed", file=sys.stderr)
            failed += 1
        elif "SECOND_LINE" not in glued_lines:
            print("FAIL addHeaders.sh glued: following line removed", file=sys.stderr)
            failed += 1
        elif glued.read_text().count("END_HEADER") != 1:
            print(
                "FAIL addHeaders.sh glued: END_HEADER was duplicated",
                file=sys.stderr,
            )
            failed += 1
        else:
            again = subprocess.run(
                ["bash", str(script)],
                cwd=repo,
                capture_output=True,
                text=True,
                env=env,
            )
            if glued.read_text().splitlines() != glued_lines or again.returncode != 0:
                print(
                    "FAIL addHeaders.sh glued: second run changed it",
                    file=sys.stderr,
                )
                failed += 1
        if missing.read_text() != original_missing:
            print(
                "FAIL addHeaders.sh missing END_HEADER: file changed",
                file=sys.stderr,
            )
            failed += 1
        elif "no END_HEADER" not in result.stderr:
            print(
                "FAIL addHeaders.sh missing END_HEADER: no warning",
                file=sys.stderr,
            )
            failed += 1
        if result.returncode != 0:
            print("FAIL addHeaders.sh EXIT:%s" % result.returncode, file=sys.stderr)
            failed += 1
    return failed


def main():
    updater = load_add_headers()
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for src, cchar, suffix in SAMPLES:
            original = src.read_text()
            expected = suffix_after(original, cchar)
            if expected is None:
                print(
                    "FAIL %s: sample has no END_HEADER line" % src.name,
                    file=sys.stderr,
                )
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
                print(
                    "FAIL %s: second run changed the file" % src.name,
                    file=sys.stderr,
                )
                failed += 1
                continue
            if check_syntax(dest, suffix, src) != 0:
                failed += 1
        failed += check_crlf_and_shebang(updater, tmp)
    failed += check_shell_updater()
    if failed:
        print("failed %s" % failed)
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
