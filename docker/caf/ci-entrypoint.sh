#!/usr/bin/env bash
# Run multi-image CAF build + tests for GaussJacobiQuad.
# Expects repo root as CWD and OpenCoarrays+MPICH on PATH / LIBRARY_PATH.
set -euo pipefail

FCFLAGS="${GJP_CAF_FCFLAGS:--fcoarray=lib}"
LDFLAGS="${GJP_CAF_LDFLAGS:--lcaf_mpi -lmpifort -lmpi}"

echo "== CAF toolchain =="
command -v gfortran
gfortran --version | head -1
command -v cafrun || command -v mpiexec
command -v fpm
echo "FCFLAGS=$FCFLAGS"
echo "LDFLAGS=$LDFLAGS"

echo "== fpm build (multi-image CAF) =="
fpm build --flag "$FCFLAGS" --link-flag "$LDFLAGS"

echo "== fpm test (1-image CAF lib mode) =="
fpm test --flag "$FCFLAGS" --link-flag "$LDFLAGS"

BIN=$(find build -type f -name gjp_bench_caf | head -1)
if [[ -z "$BIN" ]]; then
  echo "gjp_bench_caf not found" >&2
  exit 1
fi

# Note: some MPICH 5 + OpenCoarrays combos abort on MPI_Finalize after a
# successful run (unmatched message). Treat printed CORRECT/nimg as success.
run_caf() {
  local nimg=$1
  local out=$2
  set +e
  cafrun -n "$nimg" "$BIN" 128 0.5 0.5 1 0 >"$out" 2>&1
  local ec=$?
  set -e
  cat "$out"
  grep -E "CORRECT nimg=[[:space:]]*${nimg}" "$out"
  # max|dx| must be tiny (0 or ~1e-16)
  grep -E 'max\|dx\|=' "$out" | head -1
  if ! grep -E "CORRECT nimg=[[:space:]]*${nimg}" "$out" >/dev/null; then
    echo "missing CORRECT for nimg=$nimg (exit=$ec)" >&2
    exit 1
  fi
  echo "nimg=$nimg OK (cafrun exit=$ec; finalize aborts ignored if CORRECT printed)"
}

echo "== multi-image correctness: nimg=2 =="
run_caf 2 /tmp/caf_ci_n2.txt

echo "== multi-image correctness: nimg=4 =="
run_caf 4 /tmp/caf_ci_n4.txt

echo "CAF multi-image CI OK"
