#!/bin/sh
# The uninstalled pc must publish the module directory and find gjp_version.mod.
set -eu
build=$1
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
export PKG_CONFIG_PATH="$build/meson-uninstalled${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
moddir=$(pkg-config --variable=moduledir gjp_cinterp-uninstalled)
test -n "$moddir"
test -f "$moddir/gjp_version.mod"
gfortran -o "$build/use_gjp_version" "$root/ci/use_gjp_version.f90" \
    $(pkg-config --cflags --libs gjp_cinterp-uninstalled)
