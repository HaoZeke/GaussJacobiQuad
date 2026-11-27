#!/bin/sh
# Link a C consumer with the shared pkg-config line, not --static.
set -eu
build=$1
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
export PKG_CONFIG_PATH="$build/meson-uninstalled${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
cc -o "$build/c_pkgconfig_consumer" "$root/ci/c_pkgconfig_consumer.c" \
    $(pkg-config --cflags --libs gjp_cinterp-uninstalled)
