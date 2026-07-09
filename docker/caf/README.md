# Multi-image CAF toolchain image (nix + dockerTools)

Reproducible **gfortran + MPICH + OpenCoarrays** stack for GaussJacobiQuad
`rec_caf`. OpenMPI is **not** used: OpenCoarrays refuses OpenMPI with
gfortran ≥14.

## Packages (flake)

| Attr | What |
|------|------|
| `.#opencoarrays` | `libcaf_mpi` + `caf` / `cafrun` against nixpkgs MPICH |
| `.#tools` | Full env: gfortran, mpich, OpenCoarrays, fpm, openblas, pytest |
| `.#oci` | OCI image (`docker load < result`) |
| `.#oci-stream` | Stream script: `./result | docker load` |
| `.#devShells.default` | Interactive multi-image CAF shell |

## Build locally (needs Nix flakes)

```bash
cd docker/caf
# library only
nix build .#opencoarrays -L

# multi-image smoke (num_images=2)
nix build .#checks.x86_64-linux.opencoarrays-probe -L   # or: nix flake check

# Docker image
nix build .#oci -o gjp-caf-image
docker load < gjp-caf-image
docker run --rm -v "$PWD/../..":/work -w /work \
  ghcr.io/haozeke/gaussjacobiquad-caf:latest \
  bash -lc 'fpm test --flag "-fcoarray=lib" --link-flag "-lcaf_mpi -lmpifort -lmpi"'
```

Dev shell (no Docker):

```bash
cd docker/caf && nix develop
# then from repo root:
fpm build --flag "$GJP_CAF_FCFLAGS" --link-flag "$GJP_CAF_LDFLAGS"
cafrun -n 4 ./build/*/app/gjp_bench_caf 4096 0.5 0.5 3 1
```

## CI

Workflow `.github/workflows/caf_multiimage.yml`:

1. Installs Nix on `ubuntu-latest`
2. Builds `docker/caf#oci` and loads it into Docker **or** uses `nix develop` /
   `nix shell` for the same toolchain
3. Builds GaussJacobiQuad with `-fcoarray=lib` + OpenCoarrays link flags
4. Runs `fpm test` and a short `cafrun -n 2` / `cafrun -n 4` bench correctness pass

Optional GHCR publish (on `main` when `docker/caf/**` changes) uses
`docker/caf#oci` and `docker push ghcr.io/haozeke/gaussjacobiquad-caf`.

## Why this layout

- **MPICH**, not OpenMPI, for gfortran CAF multi-image
- Recipe lives in-tree so CI and developers share one toolchain definition
- `dockerTools.buildLayeredImage` = pure Nix OCI (no host Dockerfile stages)

## Quirks

- OpenCoarrays **refuses OpenMPI** with gfortran ≥14; this flake hard-wires **MPICH**.
- Release **2.10.2** is too old for gfortran ≥15 CAF ABI; `opencoarrays.nix` pins
  commit `74a5d0a` (post-2.10.2) which provides `caf_get_from_remote` / accessor registration.
- nixpkgs package name `fpm` is the **Ruby** packaging tool — this flake ships
  **fortran-lang fpm** as `fortran-fpm` / `fpm` on `PATH` via `fortran-fpm.nix`.
- OCI `buildEnv` must not list packages that re-export the same files: do **not**
  add `gcc` or `binutils` next to `gfortran` (shared `ld.bfd`), and do **not**
  add a separate `lapack` next to `openblas` (shared `liblapack.so.3`).
