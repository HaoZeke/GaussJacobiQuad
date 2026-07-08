# gauss-jacobi-quad (PyPI)

Python bindings for [GaussJacobiQuad](https://github.com/HaoZeke/GaussJacobiQuad).
Calls the Fortran library through **`libgjp_cinterp`** (ctypes). The quadrature
math is not reimplemented in Python.

## Install

```bash
pip install gauss-jacobi-quad
```

You must also provide a built shared library:

```bash
# from GaussJacobiQuad sources
meson setup bbdir -Dfcoarray=single && meson compile -C bbdir
export GJP_CINTERP=$PWD/bbdir/libgjp_cinterp.so
export LD_LIBRARY_PATH=$PWD/bbdir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
```

## Usage

```python
from gauss_jacobi_quad import gauss_jacobi

x, w = gauss_jacobi(32, 0.0, 0.0)  # auto
x, w = gauss_jacobi(40, 0.5, 0.5, method="glr")
```

CLI:

```bash
gauss-jacobi-quad --npts 32 --alpha 0 --beta 0 --meth auto
```

## License

MIT
