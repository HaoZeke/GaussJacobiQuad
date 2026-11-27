#!/usr/bin/env python3
"""Stage or install gjp_version.mod into the compiler-specific include directory."""

from __future__ import annotations

import os
import sys
from pathlib import Path, PurePath


def install_destination(prefix: Path, arg: str) -> Path:
    """Join prefix and includedir. An absolute includedir keeps prefix.

    Meson rebases an absolute install path onto DESTDIR with
    PurePath(destdir, *PurePath(path).parts[1:]). pathlib's `/` drops the
    left operand when the right one is absolute, which is what used to
    ignore MESON_INSTALL_DESTDIR_PREFIX.
    """
    raw = PurePath(arg)
    if raw.is_absolute():
        return Path(PurePath(prefix, *raw.parts[1:]))
    return prefix / raw


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


def install_tree(includedir: str) -> int:
    build = Path(os.environ["MESON_BUILD_ROOT"])
    prefix = Path(
        os.environ.get(
            "MESON_INSTALL_DESTDIR_PREFIX", os.environ["MESON_INSTALL_PREFIX"]
        )
    )
    dest = install_destination(prefix, includedir)
    hits: list[Path] = []
    seen: set[Path] = set()
    for hit in build.glob("libGaussJacobiQuad*.p/gjp_version.mod"):
        if hit.name.lower() != "gjp_version.mod":
            continue
        resolved = hit.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        hits.append(hit)
    if not hits:
        sys.stderr.write(f"expected a gjp_version.mod under {build}, found {hits}\n")
        return 1
    try:
        mod = _same_mod(hits)
    except SystemExit as exc:
        sys.stderr.write(f"{exc}\n")
        return 1
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / mod.name
    target.write_bytes(mod.read_bytes())
    print(f"Installing {mod} to {target}")
    return 0


def main(argv: list[str]) -> int:
    if len(argv) >= 4 and argv[1] == "--stage":
        return stage([Path(item) for item in argv[2:-1]], Path(argv[-1]))
    if len(argv) == 2:
        return install_tree(argv[1])
    sys.stderr.write(
        "usage: install_gjp_mod.py INCLUDE/REL/DIR\n"
        "       install_gjp_mod.py --stage LIB [LIB...] OUTPUT\n"
    )
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
