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
"""Accept and reject cases for the release-tag SemVer check.

Default mode calls scripts/semver_tag.py, the same entry point as
.github/workflows/release.yml. Pass --regex to evaluate another pattern
(python3 re) instead, which is how the previous grep -E tag check is compared.
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "semver_tag.py"
ROOT = SCRIPT.parents[1]

ACCEPT = [
    "v1.0.0",
    "v1.0.0-rc.1",
    "v1.0.0+build.5",
    "v1.0.0-alpha.1+sha.abc",
    "v10.20.30",
]
REJECT = [
    "v01.0.0",
    "v1.0.0-01",
    "v1.0",
    "v1.0.0-",
    "v1.0.0+",
    "1.0.0",
]


def cli_matches(tag: str) -> bool:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), tag],
        capture_output=True,
        text=True,
    )
    if proc.returncode not in (0, 1):
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode or 1)
    return proc.returncode == 0


def run_tags(matches) -> int:
    failed = 0
    for tag in ACCEPT:
        if matches(tag):
            print(f"ok accept {tag}")
        else:
            print(f"FAIL accept {tag}")
            failed += 1
    for tag in REJECT:
        if matches(tag):
            print(f"FAIL reject {tag}")
            failed += 1
        else:
            print(f"ok reject {tag}")
    total = len(ACCEPT) + len(REJECT)
    print(f"RESULT fail={failed} total={total}")
    return 1 if failed else 0


def run_notes() -> int:
    sample = (
        "# 1.2.3 (day)\n"
        "alpha\n"
        "# 01.0.0 (day)\n"
        "still-leading-zero\n"
        "# 1.0.0-01 (day)\n"
        "still-pre-zero\n"
        "# 1.0.0+build.5 (day)\n"
        "beta\n"
        "0.1.0 (day)\n"
        "old\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "NEWS.md"
        path.write_text(sample, encoding="utf-8")
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "--notes", "1.2.3", str(path)],
            capture_output=True,
            text=True,
        )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        print("FAIL notes exit")
        print("NOTES fail=1 total=1")
        return 1
    expected = (
        "# 1.2.3 (day)\n"
        "alpha\n"
        "# 01.0.0 (day)\n"
        "still-leading-zero\n"
        "# 1.0.0-01 (day)\n"
        "still-pre-zero\n"
    )
    if proc.stdout != expected:
        print("FAIL notes boundary")
        print(proc.stdout)
        print("NOTES fail=1 total=1")
        return 1
    print("ok notes boundary")
    print("NOTES fail=0 total=1")
    return 0


def run_workflow_wiring() -> int:
    text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    failed = 0
    if "python3 scripts/semver_tag.py" not in text:
        print("FAIL workflow does not call scripts/semver_tag.py")
        failed += 1
    old_tag = r"^v[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z-]+(\.[0-9A-Za-z-]+)*)?$"
    if old_tag in text:
        print("FAIL workflow still uses the old tag grep")
        failed += 1
    # Prerelease flagging stays the hyphen check it already had.
    prerelease = r"grep -Eq '^v[0-9]+\.[0-9]+\.[0-9]+-'"
    if prerelease not in text:
        print("FAIL prerelease flag regex changed")
        failed += 1
    old_boundary = r"[0-9]+\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z-]+(\\.[0-9A-Za-z-]+)*)?"
    if old_boundary in text:
        print("FAIL workflow still uses the old NEWS boundary")
        failed += 1
    if failed == 0:
        print("ok workflow wiring")
    print(f"WIRING fail={failed} total=4")
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--regex",
        help="python3 re used instead of scripts/semver_tag.py (tag cases only)",
    )
    args = parser.parse_args()
    if args.regex is not None:
        compiled = re.compile(args.regex)
        return run_tags(lambda tag: compiled.fullmatch(tag) is not None)
    status = run_tags(cli_matches)
    status |= run_notes()
    status |= run_workflow_wiring()
    return status


if __name__ == "__main__":
    sys.exit(main())
