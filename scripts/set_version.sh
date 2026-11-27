#!/usr/bin/env bash
# Bump version strings across release-managed files. Invoked by cog pre_bump_hooks.
set -euo pipefail

version="${1:?usage: set_version.sh <semver>}"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

# Same rules as the generated macros: rc suffixes, a missing patch, and
# the 0..99 compact bound. Fails before any file is rewritten.
parsed="$(python3 "$root/scripts/parse_gjp_version.py" "$version")"
read -r major minor patch <<<"$parsed"

# fpm.toml: version = "X.Y.Z"
sed -i "s/^version = \".*\"/version = \"${version}\"/" fpm.toml

# meson.build: version: 'X.Y.Z'
sed -i "s/version: '[^']*'/version: '${version}'/" meson.build

# Doxygen project number
sed -i "s/^PROJECT_NUMBER[[:space:]]*=.*/PROJECT_NUMBER          = \"${version}\"/" \
  apidocs/Doxygen-GaussJacobiQuad.cfg

# CITATION.cff tracks the release tag form
sed -i "s/^version: .*/version: v${version}/" CITATION.cff

# fpm compiles this file; there is no configure_file on that path.
sed -i \
  -e "s/^integer, parameter :: gjp_version_major = .*/integer, parameter :: gjp_version_major = ${major}/" \
  -e "s/^integer, parameter :: gjp_version_minor = .*/integer, parameter :: gjp_version_minor = ${minor}/" \
  -e "s/^integer, parameter :: gjp_version_patch = .*/integer, parameter :: gjp_version_patch = ${patch}/" \
  -e "s/^character(len=\\*), parameter :: gjp_version_string = .*/character(len=*), parameter :: gjp_version_string = \"${version}\"/" \
  src/gjp_version.f90

echo "set_version.sh: updated project files to ${version}"
