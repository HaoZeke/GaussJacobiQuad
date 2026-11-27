#!/bin/bash
# src/*.f90 goes through add_headers.py. Only a whole-line banner in the
# leading comment block is rewritten; marker text elsewhere is left as-is.
GITROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$SCRIPT_DIR/add_headers.py" --dirs "$GITROOT/src" --ftypes f90 --cchar '!'
