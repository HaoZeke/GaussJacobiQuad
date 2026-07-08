# 0.2.0 (2026-07-08)

## Features

- Add Bogaert Legendre asymptotics (`bogaert`), Glaser–Liu–Rokhlin-style Prüfer+Newton (`glr`/`glr_caf`), and a single public entry `gauss_jacobi_rule` with optional method / `auto` selection via `select_method_auto`. (bogaert_glr_auto)
- Improve C/C++ ergonomics (`gauss_jacobi_rule_c` with status codes and auto method) and Python ctypes bindings (`gauss_jacobi` returning numpy arrays). (c_py_bindings)
- Add Coarray batch API `gauss_jacobi_batch_caf` so independent rules for rec, gw, and algo665 are partitioned across images (multi-image speedup for eigensolve methods). (caf_batch)
- Add Coarray Fortran path `rec_caf`: partition independent recurrence Newton nodes across images and gather with coarrays (`this_image` / `num_images` / `sync all`). Build with `-fcoarray=single` (or multi-image CAF when available). (caf_rec)
- Add crates.io Rust crate `gauss-jacobi-quad` and PyPI package `gauss-jacobi-quad` over the C ABI, with OIDC publish workflow. (rust_pypi_packages)
- Add `sturm` / `sturm_caf` Golub–Welsch paths: Sturm bisection of the Jacobi tridiagonal plus inverse iteration, with Coarray ownership over eigenvalue indices (CAF-friendly alternative to sequential `imtqlx`). (sturm)
- Add SymPy-based accuracy matrix (`scripts/sympy_accuracy.py`) covering all single-image methods across mild, high-α/β, Legendre, and auto regimes. (sympy_accuracy)
- Added a variant of algorithm 665, refactored into gjp_common to
  accentuate similarities between 665 and standard GW (7.algo665)


0.1.0 (2023-08-26)

### Features

- Added a CLI for validating against SciPy
- Added an ISO_C_BINDING compatible interface and a C header
- Added an f2py generated wrapper for Python interoperability
- Implemented the Golub-Welsch algorithm with newer LAPACK routines
- Wrote a CLI for obtaining formatted results from SymPy

### Improved Documentation

- Added benchmarks and a discussion on the interfaces (#8)

### Misc

-
