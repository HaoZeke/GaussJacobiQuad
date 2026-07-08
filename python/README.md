# gauss-jacobi-quad

CPython **extension** for [GaussJacobiQuad](https://github.com/HaoZeke/GaussJacobiQuad).

- Calls **ISO_C_BINDING** entry `gauss_jacobi_rule_c` (same ABI as the C library)
- **Stable ABI** (`Py_LIMITED_API` 3.9 → one `cp39-abi3` manylinux wheel for 3.9+ GIL builds)
- **Free-threaded** wheels (`cp313t`, full C API + `Py_MOD_GIL_NOT_USED`) — separate ABI, cannot share abi3
- **Multi-phase** module init (PEP 489)
- Heap exception type via `PyErr_NewException` (dynamic type)
- Wheel embeds the Fortran/C kernels (meson-python + gfortran)

## Install

```bash
pip install gauss-jacobi-quad
```

Building from source requires `gfortran`, `openblas`/`lapack`, `meson`, `ninja`.

## Usage

```python
from gauss_jacobi_quad import gauss_jacobi, GaussJacobiError

x, w = gauss_jacobi(32, 0.0, 0.0)                 # auto
x, w = gauss_jacobi(40, 0.5, 0.5, method="glr")
```

No `GJP_CINTERP` needed for the extension wheel.

## License

MIT
