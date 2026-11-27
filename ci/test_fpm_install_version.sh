#!/bin/sh
# A clean fpm install must ship gjp_version.mod and link a program that uses it.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
prefix=$(mktemp -d)
build=$(mktemp -d)
trap 'rm -rf "$prefix" "$build"' EXIT
(
    cd "$root"
    fpm install --prefix "$prefix" --build-dir "$build" --flag "-fcoarray=single"
)
mod=$(find "$prefix" -name gjp_version.mod | head -n 1)
if [ -z "$mod" ]; then
    echo "gjp_version.mod missing under $prefix" >&2
    find "$prefix" -type f >&2
    exit 1
fi
lib=$(find "$prefix" -name libGaussJacobiQuad.a | head -n 1)
if [ -z "$lib" ]; then
    echo "libGaussJacobiQuad.a missing under $prefix" >&2
    exit 1
fi
gfortran -o "$prefix/use_gjp_version" "$root/ci/use_gjp_version.f90" \
    -I"$(dirname "$mod")" "$lib"
"$prefix/use_gjp_version"
