#!/usr/bin/env python3
"""Install gjp_version.mod into the compiler-specific include directory."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        sys.stderr.write("usage: install_gjp_mod.py INCLUDE/REL/DIR\n")
        return 2
    build = Path(os.environ["MESON_BUILD_ROOT"])
    prefix = Path(os.environ.get("MESON_INSTALL_DESTDIR_PREFIX", os.environ["MESON_INSTALL_PREFIX"]))
    dest = prefix / sys.argv[1]
    found = [path for path in build.rglob("*.mod") if path.name.lower() == "gjp_version.mod"]
    if len(found) != 1:
        sys.stderr.write(f"expected one gjp_version.mod under {build}, found {found}\n")
        return 1
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / found[0].name
    target.write_bytes(found[0].read_bytes())
    print(f"Installing {found[0]} to {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
