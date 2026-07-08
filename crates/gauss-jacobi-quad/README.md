# gauss-jacobi-quad

Rust bindings for [GaussJacobiQuad](https://github.com/HaoZeke/GaussJacobiQuad):
Gauss–Jacobi quadrature nodes and weights via the **C ABI** of `libgjp_cinterp`
(Fortran kernels). This crate does **not** reimplement the math.

## Install

```toml
[dependencies]
gauss-jacobi-quad = "0.2"
```

Build the native library first (from the GaussJacobiQuad repo):

```bash
meson setup bbdir -Dfcoarray=single
meson compile -C bbdir
export GJP_CINTERP=$PWD/bbdir/libgjp_cinterp.so
export LD_LIBRARY_PATH=$PWD/bbdir:$LD_LIBRARY_PATH
```

## Usage

```rust
use gauss_jacobi_quad::{gauss_jacobi, Method};

let (x, w) = gauss_jacobi(32, 0.0, 0.0, Method::Auto)?;
let (x, w) = gauss_jacobi(40, 0.5, 0.5, Method::Glr)?;
```

## License

MIT
