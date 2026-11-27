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
"""Stage gjp_version.mod into the build root for the uninstalled .pc file."""

from __future__ import annotations

import sys
from pathlib import Path


def _same_mod(hits: list[Path]) -> Path:
    if not hits:
        raise SystemExit("no gjp_version.mod")
    data = hits[0].read_bytes()
    for hit in hits[1:]:
        if hit.read_bytes() != data:
            raise SystemExit(f"gjp_version.mod copies differ: {hits}")
    return hits[0]


def find_mod(libs: list[Path]) -> Path:
    """The library private dir only. Do not rglob the whole build tree."""
    hits: list[Path] = []
    seen: set[Path] = set()
    for lib in libs:
        parent = lib.parent
        patterns = (
            lib.name + ".p/gjp_version.mod",
            "libGaussJacobiQuad*.p/gjp_version.mod",
        )
        for pattern in patterns:
            for hit in parent.glob(pattern):
                if hit.name.lower() != "gjp_version.mod":
                    continue
                resolved = hit.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                hits.append(hit)
    if not hits:
        raise SystemExit(f"no gjp_version.mod beside {libs}")
    return _same_mod(hits)


def stage(libs: list[Path], output: Path) -> int:
    mod = find_mod(libs)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(mod.read_bytes())
    print(f"Staging {mod} to {output}")
    return 0


def main(argv: list[str]) -> int:
    if len(argv) >= 4 and argv[1] == "--stage":
        return stage([Path(item) for item in argv[2:-1]], Path(argv[-1]))
    sys.stderr.write("usage: install_gjp_mod.py --stage LIB [LIB...] OUTPUT\n")
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
