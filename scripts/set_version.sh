#!/usr/bin/env bash
# Bump version strings across release-managed files. Invoked by cog pre_bump_hooks.
set -euo pipefail

version="${1:?usage: set_version.sh <semver>}"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

# Python rewrites the files. fpm.toml receives only major.minor.patch because
# fpm rejects pre-release text. towncrier.toml receives the full version,
# including any rc or pre-release suffix. The same path is what BSD sed
# would need instead of GNU sed -i.
python3 "$root/scripts/parse_gjp_version.py" --apply "$version"

echo "set_version.sh: updated project files to ${version}"
