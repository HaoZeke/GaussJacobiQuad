#!/bin/sh
# pyproject.toml configures the wheel with -Dpython_package=true. That
# configuration must compile, including src/gjp_version.f90.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$(mktemp -d)
trap 'rm -rf "$build"' EXIT
meson setup "$build" "$root" -Dpython_package=true -Dwarning_level=1 -Dbuildtype=release
meson compile -C "$build"
