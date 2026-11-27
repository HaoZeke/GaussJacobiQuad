# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
"""Checked-in versions: fpm triple, and meson, Fortran, and towncrier text."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import parse_gjp_version as parser  # noqa: E402


def quoted_version(path: Path, pattern: str) -> str:
    match = re.search(pattern, path.read_text())
    if not match:
        raise SystemExit(f"no version in {path}")
    return match.group(1)


def main() -> int:
    source = ROOT / "src" / "gjp_version.f90"
    if not source.is_file():
        raise SystemExit(f"missing {source}")
    fpm_version = quoted_version(ROOT / "fpm.toml", r'(?m)^version\s*=\s*"([^"]+)"')
    meson_version = quoted_version(ROOT / "meson.build", r"version:\s*'([^']+)'")
    town_version = quoted_version(
        ROOT / "towncrier.toml", r'(?m)^version\s*=\s*"([^"]+)"'
    )
    major, minor, patch, text = parser.fortran_constants(source)
    fpm_parts = parser.parse_version(fpm_version)
    meson_parts = parser.parse_version(meson_version)
    town_parts = parser.parse_version(town_version)
    triple = f"{major}.{minor}.{patch}"
    if (
        fpm_parts != (major, minor, patch)
        or meson_parts != (major, minor, patch)
        or town_parts != (major, minor, patch)
    ):
        raise SystemExit(
            f"components fpm {fpm_parts} meson {meson_parts} "
            f"towncrier {town_parts} module {(major, minor, patch)}"
        )
    # fpm rejects pre-release text, so fpm.toml stores only the release triple.
    if fpm_version != triple:
        raise SystemExit(f"fpm.toml {fpm_version!r} != release triple {triple}")
    if text != meson_version:
        raise SystemExit(f"{source} string {text!r} != meson {meson_version!r}")
    # towncrier build without --version uses this string, suffix included.
    if town_version != meson_version:
        raise SystemExit(f"towncrier.toml {town_version!r} != meson {meson_version!r}")
    print(
        f"agree fpm {fpm_version} meson {meson_version} "
        f"towncrier {town_version} -> {triple}"
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except parser.VersionError as exc:
        sys.stderr.write(f"error: {exc}\n")
        sys.exit(1)
