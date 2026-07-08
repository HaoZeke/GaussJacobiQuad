//! Integration tests require a built libgjp_cinterp (set GJP_CINTERP).
use gauss_jacobi_quad::{gauss_jacobi, gauss_jacobi_method_str, Method};

#[test]
fn bogaert_and_auto_agree_legendre() {
    let (xa, wa) = gauss_jacobi(32, 0.0, 0.0, Method::Auto).unwrap();
    let (xb, wb) = gauss_jacobi(32, 0.0, 0.0, Method::Bogaert).unwrap();
    for i in 0..32 {
        assert!((xa[i] - xb[i]).abs() < 1e-14);
        assert!((wa[i] - wb[i]).abs() < 1e-14);
    }
}

#[test]
fn method_str_sturm() {
    let (x, w) = gauss_jacobi_method_str(24, 0.5, 0.5, "sturm").unwrap();
    assert_eq!(x.len(), 24);
    assert!(w.iter().all(|&wi| wi > 0.0));
}
