"""Finding 6: the checked-in module matches fpm.toml and the meson project version."""

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
    if fpm_version != meson_version:
        raise SystemExit(f"fpm.toml {fpm_version} != meson.build {meson_version}")
    major, minor, patch, text = parser.fortran_constants(source)
    if text != fpm_version:
        raise SystemExit(f"{source} string {text!r} != project version {fpm_version!r}")
    parsed = parser.parse_version(fpm_version)
    if (major, minor, patch) != parsed:
        raise SystemExit(f"{source} components {(major, minor, patch)} != {parsed}")
    print(f"agree {fpm_version} -> {major}.{minor}.{patch}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except parser.VersionError as exc:
        sys.stderr.write(f"error: {exc}\n")
        sys.exit(1)
