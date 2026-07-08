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
//! use gauss_jacobi_quad::{gauss_jacobi, Method};
//!
//! let (x, w) = gauss_jacobi(32, 0.0, 0.0, Method::Auto).unwrap();
//! assert_eq!(x.len(), 32);
//! assert!(w.iter().all(|&wi| wi > 0.0));
//! ```

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_int};
use std::ptr;

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
    Unknown = -1,
}

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
            Status::ErrMethod => "unknown method",
            Status::ErrBogaertAb => "bogaert requires alpha=beta=0",
            Status::ErrBogaertN => "bogaert requires n >= 21",
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

    fn gjp_status_string(status: c_int) -> *const c_char;
}

/// Human-readable status from the C library (static C string).
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
    let method_ptr = c_method
        .as_ref()
        .map(|c| c.as_ptr())
        .unwrap_or(ptr::null());

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
    let method_ptr = c_method
        .as_ref()
        .map(|c| c.as_ptr())
        .unwrap_or(ptr::null());
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

#[cfg(test)]
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
        let mut xg = xg;
        let mut xa = xa;
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
}
