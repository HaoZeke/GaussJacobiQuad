#!/usr/bin/env bash
# Bump version strings across release-managed files. Invoked by cog pre_bump_hooks.
set -euo pipefail

version="${1:?usage: set_version.sh <semver>}"

if ! [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-].*)?$ ]]; then
  echo "error: invalid semver: $version" >&2
  exit 1
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

# fpm.toml: version = "X.Y.Z"
sed -i "s/^version = \".*\"/version = \"${version}\"/" fpm.toml

# meson.build: version: 'X.Y.Z'
sed -i "s/version: '[^']*'/version: '${version}'/" meson.build

# Doxygen project number
sed -i "s/^PROJECT_NUMBER[[:space:]]*=.*/PROJECT_NUMBER          = \"${version}\"/" \
  apidocs/Doxygen-GaussJacobiQuad.cfg

# CITATION.cff tracks the release tag form
sed -i "s/^version: .*/version: v${version}/" CITATION.cff

echo "set_version.sh: updated project files to ${version}"
