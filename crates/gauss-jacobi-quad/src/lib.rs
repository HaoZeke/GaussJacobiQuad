// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// -----------------------------------------------------------------------------
// END_HEADER
//! Safe Rust bindings to **GaussJacobiQuad** via the C ABI
//! (`gauss_jacobi_rule_c` in `libgjp_cinterp`).
//!
//! The quadrature math stays in Fortran; this crate only marshals arrays and
//! status codes. Build `libgjp_cinterp` with meson first, then point
//! `GJP_CINTERP` at the `.so` (or install a `gjp_cinterp.pc`).
//!
//! # Example
//!
//! ```no_run
//! # #[cfg(feature = "native")]
//! # {
//! use gauss_jacobi_quad::{gauss_jacobi, Method};
//!
//! let (x, w) = gauss_jacobi(32, 0.0, 0.0, Method::Auto).unwrap();
//! assert_eq!(x.len(), 32);
//! assert!(w.iter().all(|&wi| wi > 0.0));
//! # }
//! ```

#[cfg(feature = "native")]
use std::ffi::{CStr, CString};
#[cfg(feature = "native")]
use std::os::raw::{c_char, c_double, c_int};
#[cfg(feature = "native")]
use std::ptr;

#[cfg(feature = "spectral")]
mod spectral;
#[cfg(feature = "spectral")]
pub use spectral::{
    cc_project, cc_synthesize, chebyshev_antiderivative, clenshaw, clenshaw_bound,
    clenshaw_curtis_panel, direct_cosine_synthesize, direct_project, direct_sine_synthesize,
    dst_project, dst_synthesize, endpoint_trapezoid_array, forward_inplace,
    gauss_legendre_two_point, half_period_grid, log_composite_trapezoid,
    log_composite_trapezoid_by, log_stable_trapezoid, midpoint_integral,
    midpoint_trapezoid_enclosure, monotone_exponential_enclosure, piecewise_exponential_integral,
    piecewise_exponential_integral_by, piecewise_exponential_panel, ChirpGrid, ChirpMap,
    ClenshawCurtisPanel, EvaluationError, ExponentialSample, HalfPeriodGrid, MapPoint,
    MonotoneExponentialEnclosure, SpectralError,
};

/// Status codes matching `GaussJacobiQuadCInterp.h`.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Ok = 0,
    ErrNpts = 1,
    ErrAlpha = 2,
    ErrBeta = 3,
    ErrMethod = 4,
    ErrBogaertAb = 5,
    ErrBogaertN = 6,
    ErrBounds = 7,
    Unknown = -1,
}

#[cfg(feature = "native")]
impl From<c_int> for Status {
    fn from(v: c_int) -> Self {
        match v {
            0 => Status::Ok,
            1 => Status::ErrNpts,
            2 => Status::ErrAlpha,
            3 => Status::ErrBeta,
            4 => Status::ErrMethod,
            5 => Status::ErrBogaertAb,
            6 => Status::ErrBogaertN,
            7 => Status::ErrBounds,
            _ => Status::Unknown,
        }
    }
}

impl Status {
    pub fn message(self) -> &'static str {
        match self {
            Status::Ok => "ok",
            Status::ErrNpts => "npts must be positive",
            Status::ErrAlpha => "alpha must be greater than -1",
            Status::ErrBeta => "beta must be greater than -1",
            Status::ErrMethod => {
                "unknown method (use auto|rec|gw|algo665|algo665_dc|sturm|glr|bogaert)"
            }
            Status::ErrBogaertAb => "bogaert requires alpha=beta=0 (Gauss-Legendre only)",
            Status::ErrBogaertN => "bogaert requires n >= 21",
            Status::ErrBounds => "invalid interval or non-finite scaled node or weight",
            Status::Unknown => "unknown status",
        }
    }
}

/// Error from the C library or argument validation.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("gauss_jacobi_rule_c: {status:?} ({msg}) method={method:?}")]
    Library {
        status: Status,
        msg: &'static str,
        method: Option<String>,
    },
    /// `gauss_jacobi_scale_c` returned a non-zero status.
    #[error("gauss_jacobi_scale_c: {status:?} ({msg})")]
    Scale { status: Status, msg: &'static str },
    /// `nodes` and `weights` differ in length.
    #[error("gauss_jacobi_scale: nodes and weights must have the same length")]
    ScaleLength,
    /// Empty input, or a length that does not fit in `i32`.
    #[error("gauss_jacobi_scale: npts must be positive and fit in a C int")]
    ScaleNpts,
    #[error("invalid method string (interior NUL)")]
    Nul(#[from] std::ffi::NulError),
}

/// Method selection for [`gauss_jacobi`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    /// `select_method_auto` policy (empty / `"auto"`).
    #[default]
    Auto,
    Rec,
    Gw,
    Algo665,
    Algo665Dc,
    Sturm,
    Glr,
    Bogaert,
}

impl Method {
    pub fn as_c_str(self) -> Option<&'static str> {
        match self {
            Method::Auto => None,
            Method::Rec => Some("rec"),
            Method::Gw => Some("gw"),
            Method::Algo665 => Some("algo665"),
            Method::Algo665Dc => Some("algo665_dc"),
            Method::Sturm => Some("sturm"),
            Method::Glr => Some("glr"),
            Method::Bogaert => Some("bogaert"),
        }
    }
}

#[cfg(feature = "native")]
#[link(name = "gjp_cinterp")]
extern "C" {
    fn gauss_jacobi_rule_c(
        npts: c_int,
        alpha: c_double,
        beta: c_double,
        x: *mut c_double,
        wts: *mut c_double,
        method: *const c_char,
    ) -> c_int;

    fn gauss_jacobi_scale_c(
        npts: c_int,
        alpha: c_double,
        beta: c_double,
        a: c_double,
        b: c_double,
        x: *mut c_double,
        wts: *mut c_double,
    ) -> c_int;

    fn gjp_status_string(status: c_int) -> *const c_char;
}

/// Human-readable status from the C library (static C string).
#[cfg(feature = "native")]
pub fn status_string(status: Status) -> String {
    unsafe {
        let p = gjp_status_string(status as c_int);
        if p.is_null() {
            return status.message().to_string();
        }
        CStr::from_ptr(p).to_string_lossy().into_owned()
    }
}

/// Compute Gauss–Jacobi nodes and weights.
///
/// Returns `(nodes, weights)` of length `npts`, ordered as produced by the
/// library (typically ascending nodes).
#[cfg(feature = "native")]
pub fn gauss_jacobi(
    npts: usize,
    alpha: f64,
    beta: f64,
    method: Method,
) -> Result<(Vec<f64>, Vec<f64>), Error> {
    if npts == 0 || npts > i32::MAX as usize {
        return Err(Error::Library {
            status: Status::ErrNpts,
            msg: Status::ErrNpts.message(),
            method: method.as_c_str().map(str::to_string),
        });
    }
    let mut x = vec![0.0_f64; npts];
    let mut w = vec![0.0_f64; npts];
    let c_method: Option<CString> = match method.as_c_str() {
        None => None,
        Some(s) => Some(CString::new(s)?),
    };
    let method_ptr = match c_method.as_ref() {
        Some(value) => value.as_ptr(),
        None => ptr::null(),
    };

    let st = unsafe {
        gauss_jacobi_rule_c(
            npts as c_int,
            alpha,
            beta,
            x.as_mut_ptr(),
            w.as_mut_ptr(),
            method_ptr,
        )
    };
    let status = Status::from(st);
    if status != Status::Ok {
        return Err(Error::Library {
            status,
            msg: status.message(),
            method: method.as_c_str().map(str::to_string),
        });
    }
    Ok((x, w))
}

/// Force a method by raw C string (e.g. `"sturm_caf"` is rejected by C surface;
/// use Fortran for CAF). Empty / `"auto"` select the policy.
#[cfg(feature = "native")]
pub fn gauss_jacobi_method_str(
    npts: usize,
    alpha: f64,
    beta: f64,
    method: &str,
) -> Result<(Vec<f64>, Vec<f64>), Error> {
    if npts == 0 || npts > i32::MAX as usize {
        return Err(Error::Library {
            status: Status::ErrNpts,
            msg: Status::ErrNpts.message(),
            method: Some(method.to_string()),
        });
    }
    let mut x = vec![0.0_f64; npts];
    let mut w = vec![0.0_f64; npts];
    let c_method = if method.is_empty() || method == "auto" {
        None
    } else {
        Some(CString::new(method)?)
    };
    let method_ptr = match c_method.as_ref() {
        Some(value) => value.as_ptr(),
        None => ptr::null(),
    };
    let st = unsafe {
        gauss_jacobi_rule_c(
            npts as c_int,
            alpha,
            beta,
            x.as_mut_ptr(),
            w.as_mut_ptr(),
            method_ptr,
        )
    };
    let status = Status::from(st);
    if status != Status::Ok {
        return Err(Error::Library {
            status,
            msg: status.message(),
            method: Some(method.to_string()),
        });
    }
    Ok((x, w))
}

/// Map nodes and weights from `[-1, 1]` onto `[a, b]`.
///
/// Nodes are one ties-to-even binary64 rounding of the real expression
/// `(b - a) / 2 * x + (a + b) / 2`. Weights are one ties-to-even binary64
/// rounding of `w * H**(alpha + beta + 1)`, where `H` is the binary64
/// half-width `fl((b - a) / 2)`. A half-width that flushes to zero is
/// treated as `2**-1075`. `x = -1` is exactly `a` and
/// `x = 1` is exactly `b`. A non-finite scaled node or weight is an error.
/// A weight factor that overflows on its own is not an error when the
/// product with the weight is finite.
///
/// `nodes` and `weights` must have the same length, and that length must fit
/// in `i32`. The C call is not made when the length is rejected.
#[cfg(feature = "native")]
pub fn gauss_jacobi_scale(
    nodes: &[f64],
    weights: &[f64],
    alpha: f64,
    beta: f64,
    a: f64,
    b: f64,
) -> Result<(Vec<f64>, Vec<f64>), Error> {
    if nodes.len() != weights.len() {
        return Err(Error::ScaleLength);
    }
    if nodes.is_empty() || nodes.len() > i32::MAX as usize {
        return Err(Error::ScaleNpts);
    }
    let mut x = nodes.to_vec();
    let mut w = weights.to_vec();
    let st = unsafe {
        gauss_jacobi_scale_c(
            nodes.len() as c_int,
            alpha,
            beta,
            a,
            b,
            x.as_mut_ptr(),
            w.as_mut_ptr(),
        )
    };
    let status = Status::from(st);
    if status != Status::Ok {
        return Err(Error::Scale {
            status,
            msg: status.message(),
        });
    }
    Ok((x, w))
}

#[cfg(test)]
mod status_messages {
    use super::Status;

    #[test]
    fn method_and_bogaert_match_c_strings() {
        assert_eq!(
            Status::ErrMethod.message(),
            "unknown method (use auto|rec|gw|algo665|algo665_dc|sturm|glr|bogaert)"
        );
        assert_eq!(
            Status::ErrBogaertAb.message(),
            "bogaert requires alpha=beta=0 (Gauss-Legendre only)"
        );
        assert_eq!(
            Status::ErrBounds.message(),
            "invalid interval or non-finite scaled node or weight"
        );
    }
}

#[cfg(all(test, feature = "native"))]
mod tests {
    use super::*;

    #[test]
    fn auto_legendre_positive_weights() {
        let (x, w) = gauss_jacobi(32, 0.0, 0.0, Method::Auto).expect("auto");
        assert_eq!(x.len(), 32);
        assert!(w.iter().all(|&wi| wi > 0.0 && wi.is_finite()));
        assert!(x.iter().all(|&xi| xi > -1.0 && xi < 1.0 && xi.is_finite()));
        let s: f64 = w.iter().sum();
        assert!((s - 2.0).abs() < 1e-12);
    }

    #[test]
    fn glr_matches_algo665() {
        let (xg, wg) = gauss_jacobi(40, 0.5, 0.5, Method::Glr).unwrap();
        let (xa, wa) = gauss_jacobi(40, 0.5, 0.5, Method::Algo665).unwrap();
        // sort for comparison
        let mut ig: Vec<_> = (0..xg.len()).collect();
        ig.sort_by(|&i, &j| xg[i].partial_cmp(&xg[j]).unwrap());
        let mut ia: Vec<_> = (0..xa.len()).collect();
        ia.sort_by(|&i, &j| xa[i].partial_cmp(&xa[j]).unwrap());
        for k in 0..xg.len() {
            assert!((xg[ig[k]] - xa[ia[k]]).abs() < 1e-12);
            assert!((wg[ig[k]] - wa[ia[k]]).abs() < 1e-11);
        }
    }

    #[test]
    fn bogaert_non_legendre_errors() {
        let err = gauss_jacobi(32, 0.5, 0.0, Method::Bogaert).unwrap_err();
        match err {
            Error::Library { status, .. } => assert_eq!(status, Status::ErrBogaertAb),
            _ => panic!("unexpected"),
        }
    }

    #[test]
    fn scale_rejects_bad_lengths() {
        let err = gauss_jacobi_scale(&[0.0, 1.0], &[1.0], 0.0, 0.0, 0.0, 1.0).unwrap_err();
        let text = err.to_string();
        assert!(matches!(err, Error::ScaleLength), "{text}");
        assert!(text.contains("gauss_jacobi_scale"));
        assert!(text.contains("same length"));
        assert!(!text.contains("npts must be positive"), "{text}");
        let err = gauss_jacobi_scale(&[], &[], 0.0, 0.0, 0.0, 1.0).unwrap_err();
        let text = err.to_string();
        assert!(matches!(err, Error::ScaleNpts), "{text}");
        assert!(text.contains("gauss_jacobi_scale"));
        assert!(text.contains("npts must be positive"));
        assert!(!text.contains("gauss_jacobi_rule_c"), "{text}");
    }

    #[test]
    fn scale_alpha1023_weight_one_is_bounds() {
        let err = gauss_jacobi_scale(&[0.25], &[1.0], 1023.0, 0.0, -2.0, 2.0).unwrap_err();
        let text = err.to_string();
        match err {
            Error::Scale { status, .. } => assert_eq!(status, Status::ErrBounds),
            other => panic!("scale overflow reported as {other}: {text}"),
        }
        assert!(text.starts_with("gauss_jacobi_scale_c"), "{text}");
        assert!(text.contains("non-finite scaled node or weight"), "{text}");
        assert!(!text.contains("method="), "{text}");
        let pred = f64::from_bits(1.0f64.to_bits() - 1);
        let (_x, w) = gauss_jacobi_scale(&[0.25], &[pred], 1023.0, 0.0, -2.0, 2.0).unwrap();
        assert_eq!(w[0], f64::MAX);
    }

    #[test]
    fn scale_nan_node_is_bounds() {
        let err = gauss_jacobi_scale(&[f64::NAN], &[1.0], 0.0, 0.0, 0.0, 1.0).unwrap_err();
        let text = err.to_string();
        match err {
            Error::Scale { status, .. } => assert_eq!(status, Status::ErrBounds),
            other => panic!("NaN node reported as {other}: {text}"),
        }
        assert!(text.contains("non-finite scaled node or weight"), "{text}");
    }

    #[test]
    fn scale_wide_endpoint_is_exact() {
        let (x, w) = gauss_jacobi_scale(&[-1.0], &[0.0], 0.0, 0.0, -1.0, 1.0e16).unwrap();
        assert_eq!(x[0], -1.0);
        assert_eq!(w[0], 0.0);
        let (x, _) = gauss_jacobi_scale(&[1.0], &[0.0], 0.0, 0.0, -1.0e16, 1.0).unwrap();
        assert_eq!(x[0], 1.0);
    }
}
