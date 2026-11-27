// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// -----------------------------------------------------------------------------
// END_HEADER
//! Endpoint-halved cosine quadrature, DCT-I/DST-I transforms, and mapped grids.
//!
//! These rules are pure Rust and do not require the native Gauss--Jacobi
//! library. Enable the `spectral` feature with default features disabled when
//! only these routines are needed.

use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::{Add, Div, Mul, Sub};
use std::sync::Arc;

use rustfft::num_complex::Complex;
use rustfft::{Fft, FftPlanner};

const PI: f64 = std::f64::consts::PI;

/// Invalid input to a spectral rule or mapped grid.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum SpectralError {
    #[error("the chirp centre {value:.6e} must be finite and inside (0, pi)")]
    ChirpCentre { value: f64 },
    #[error("the chirp width {value:.6e} must produce a finite positive map span")]
    ChirpWidth { value: f64 },
    #[error("{name} must be finite and inside [0, pi], got {value:.6e}")]
    MapCoordinate { name: &'static str, value: f64 },
    #[error("{name} needs at least {minimum} panel(s), got {panels}")]
    PanelCount {
        name: &'static str,
        panels: usize,
        minimum: usize,
    },
    #[error("{name}={value:.6e} must be finite and positive")]
    Positive { name: &'static str, value: f64 },
    #[error("{name} interval [{lower:.6e}, {upper:.6e}] must be finite and increasing")]
    Interval {
        name: &'static str,
        lower: f64,
        upper: f64,
    },
    #[error("{name} mode count {modes} exceeds the {panels}-panel limit")]
    ModeCount {
        name: &'static str,
        modes: usize,
        panels: usize,
    },
    #[error("{name} dimensions must be positive")]
    Dimension { name: &'static str },
    #[error("{name} has length {actual}, expected {expected}")]
    Shape {
        name: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("{name} size overflows usize arithmetic")]
    SizeOverflow { name: &'static str },
    #[error("{name} size {size} cannot be allocated")]
    Allocation { name: &'static str, size: usize },
    #[error("{name}[{index}] is not finite: {value:.6e}")]
    NonFinite {
        name: &'static str,
        index: usize,
        value: f64,
    },
    #[error("{name}[{index}]={value:.6e} is outside the required range")]
    OutOfRange {
        name: &'static str,
        index: usize,
        value: f64,
    },
    #[error("{name} series did not converge to a finite value")]
    NonConverged { name: &'static str },
}

/// A checked-rule failure or an error returned by the caller's evaluator.
#[derive(Debug, thiserror::Error)]
pub enum EvaluationError<E> {
    #[error(transparent)]
    Input(#[from] SpectralError),
    #[error("spectral-rule evaluator failed: {0}")]
    Evaluator(E),
}

const MAX_EXACT_INDEX: u128 = 1_u128 << f64::MANTISSA_DIGITS;

/// Angle, imaginary-time, and endpoint-halved weight vectors.
pub type HalfPeriodGrid = (Vec<f64>, Vec<f64>, Vec<f64>);

fn checked_nodes(
    panels: usize,
    minimum: usize,
    name: &'static str,
) -> Result<usize, SpectralError> {
    if panels < minimum {
        return Err(SpectralError::PanelCount {
            name,
            panels,
            minimum,
        });
    }
    if panels as u128 > MAX_EXACT_INDEX {
        return Err(SpectralError::SizeOverflow { name });
    }
    panels
        .checked_add(1)
        .ok_or(SpectralError::SizeOverflow { name })
}

fn checked_product(left: usize, right: usize, name: &'static str) -> Result<usize, SpectralError> {
    left.checked_mul(right)
        .ok_or(SpectralError::SizeOverflow { name })
}

fn filled<T: Clone>(size: usize, value: T, name: &'static str) -> Result<Vec<T>, SpectralError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(size)
        .map_err(|_| SpectralError::Allocation { name, size })?;
    output.resize(size, value);
    Ok(output)
}

fn validate_interval(lower: f64, upper: f64, name: &'static str) -> Result<(), SpectralError> {
    if !(lower.is_finite() && upper.is_finite() && upper > lower) {
        return Err(SpectralError::Interval { name, lower, upper });
    }
    Ok(())
}

fn validate_finite(value: f64, name: &'static str, index: usize) -> Result<(), SpectralError> {
    if !value.is_finite() {
        return Err(SpectralError::NonFinite { name, index, value });
    }
    Ok(())
}

/// Endpoint-halved uniform cosine grid on a half period.
///
/// The returned vectors are `(theta, tau, weights)`, with
/// `theta[j] = pi*j/m`, `tau[j] = beta*theta[j]/(2*pi)`, and trapezoid
/// weights `beta/m`, halved at both endpoints. This operation order is kept
/// stable because downstream spectral calculations depend on bitwise
/// reproducibility.
pub fn half_period_grid(m: usize, beta: f64) -> Result<HalfPeriodGrid, SpectralError> {
    let nodes = checked_nodes(m, 1, "half-period grid")?;
    if !(beta > 0.0 && beta.is_finite()) {
        return Err(SpectralError::Positive {
            name: "half-period beta",
            value: beta,
        });
    }
    let mut theta = filled(nodes, 0.0, "half-period theta")?;
    let mut tau = filled(nodes, 0.0, "half-period tau")?;
    let mut weights = filled(nodes, beta / m as f64, "half-period weights")?;
    for j in 0..=m {
        theta[j] = PI * j as f64 / m as f64;
        tau[j] = beta * theta[j] / (2.0 * PI);
    }
    weights[0] *= 0.5;
    weights[m] *= 0.5;
    Ok((theta, tau, weights))
}

/// Two-point Gauss--Legendre rule on an interval described by centre and width.
pub fn gauss_legendre_two_point(centre: f64, width: f64) -> Result<[(f64, f64); 2], SpectralError> {
    validate_finite(centre, "Gauss-Legendre centre", 0)?;
    if !(width > 0.0 && width.is_finite()) {
        return Err(SpectralError::Positive {
            name: "Gauss-Legendre width",
            value: width,
        });
    }
    let offset = width / (2.0 * 3.0_f64.sqrt());
    let weight = 0.5 * width;
    Ok([(centre - offset, weight), (centre + offset, weight)])
}

/// Endpoint-halved trapezoid for a fixed-size vector-valued integrand.
pub fn endpoint_trapezoid_array<const N: usize, F>(
    lower: f64,
    upper: f64,
    panels: usize,
    mut value: F,
) -> Result<[f64; N], SpectralError>
where
    F: FnMut(f64) -> [f64; N],
{
    validate_interval(lower, upper, "endpoint trapezoid")?;
    checked_nodes(panels, 1, "endpoint trapezoid")?;
    let span = upper - lower;
    let mut sum = [0.0; N];
    for index in 0..=panels {
        let weight_factor = if index == 0 || index == panels {
            0.5
        } else {
            1.0
        };
        let weight = weight_factor * span / panels as f64;
        let sampled = value(lower + span * index as f64 / panels as f64);
        for component in 0..N {
            validate_finite(sampled[component], "endpoint sample", index * N + component)?;
            sum[component] += weight * sampled[component];
        }
    }
    Ok(sum)
}

thread_local! {
    static FORWARD: RefCell<PlanCache> = RefCell::new(PlanCache::forward());
    static INVERSE: RefCell<PlanCache> = RefCell::new(PlanCache::inverse());
}

struct PlanCache {
    planner: FftPlanner<f64>,
    plans: HashMap<usize, Arc<dyn Fft<f64>>>,
    inverse: bool,
}

impl PlanCache {
    fn forward() -> Self {
        Self {
            planner: FftPlanner::new(),
            plans: HashMap::new(),
            inverse: false,
        }
    }

    fn inverse() -> Self {
        Self {
            planner: FftPlanner::new(),
            plans: HashMap::new(),
            inverse: true,
        }
    }

    fn plan(&mut self, len: usize) -> Arc<dyn Fft<f64>> {
        if let Some(found) = self.plans.get(&len) {
            return Arc::clone(found);
        }
        let plan = if self.inverse {
            self.planner.plan_fft_inverse(len)
        } else {
            self.planner.plan_fft_forward(len)
        };
        self.plans.insert(len, Arc::clone(&plan));
        plan
    }
}

fn validate_modes(
    m: usize,
    k_count: usize,
    d: usize,
    name: &'static str,
) -> Result<usize, SpectralError> {
    let nodes = checked_nodes(m, 1, name)?;
    if d == 0 {
        return Err(SpectralError::Dimension { name });
    }
    if k_count > nodes {
        return Err(SpectralError::ModeCount {
            name,
            modes: k_count,
            panels: m,
        });
    }
    Ok(nodes)
}

fn validate_transform(
    actual: usize,
    m: usize,
    k_count: usize,
    d: usize,
    name: &'static str,
) -> Result<usize, SpectralError> {
    let nodes = validate_modes(m, k_count, d, name)?;
    let expected = checked_product(nodes, d, name)?;
    if actual != expected {
        return Err(SpectralError::Shape {
            name,
            expected,
            actual,
        });
    }
    Ok(nodes)
}

fn validate_coefficients(
    actual: usize,
    m: usize,
    k_count: usize,
    d: usize,
    name: &'static str,
) -> Result<usize, SpectralError> {
    let nodes = validate_modes(m, k_count, d, name)?;
    let expected = checked_product(k_count, d, name)?;
    if actual != expected {
        return Err(SpectralError::Shape {
            name,
            expected,
            actual,
        });
    }
    Ok(nodes)
}

/// Unnormalized forward FFT with a plan cached on the calling thread.
pub fn forward_inplace(buf: &mut [Complex<f64>]) -> Result<(), SpectralError> {
    let len = buf.len();
    checked_nodes(len, 1, "FFT buffer")?;
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        plan.process(buf);
    });
    Ok(())
}

/// Weighted cosine projection on an endpoint-halved grid of `m` panels.
pub fn cc_project(
    samples: &[f64],
    m: usize,
    k_count: usize,
    d: usize,
    beta: f64,
) -> Result<Vec<f64>, SpectralError> {
    let nodes = validate_transform(samples.len(), m, k_count, d, "cosine samples")?;
    if !(beta > 0.0 && beta.is_finite()) {
        return Err(SpectralError::Positive {
            name: "cosine projection beta",
            value: beta,
        });
    }
    if m < 2 {
        return direct_project(samples, m, k_count, d, beta);
    }
    let output_len = checked_product(k_count, d, "cosine projection")?;
    let mut out = filled(output_len, 0.0, "cosine projection")?;
    let scale = beta / (2.0 * m as f64);
    validate_finite(scale, "cosine projection scale", 0)?;
    let wave_len = checked_product(2, m, "cosine FFT")?;
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(wave_len);
        let mut wave = filled(wave_len, Complex::new(0.0, 0.0), "cosine FFT workspace")?;
        for a in 0..d {
            wave[0] = Complex::new(samples[a], 0.0);
            wave[m] = Complex::new(samples[m * d + a], 0.0);
            for t in 1..m {
                let value = Complex::new(samples[t * d + a], 0.0);
                wave[t] = value;
                wave[wave_len - t] = value;
            }
            plan.process(&mut wave);
            for k in 0..k_count {
                out[k * d + a] = wave[k].re * scale;
            }
        }
        Ok::<(), SpectralError>(())
    })?;
    let _ = nodes;
    Ok(out)
}

/// Cosine-series synthesis on an endpoint-halved grid.
pub fn cc_synthesize(
    coeff: &[f64],
    k_count: usize,
    d: usize,
    m: usize,
) -> Result<Vec<f64>, SpectralError> {
    let nodes = validate_coefficients(coeff.len(), m, k_count, d, "cosine coefficients")?;
    if m < 2 {
        return direct_cosine_synthesize(coeff, k_count, d, m);
    }
    let output_len = checked_product(nodes, d, "cosine synthesis")?;
    let mut out = filled(output_len, 0.0, "cosine synthesis")?;
    let len = checked_product(2, m, "cosine inverse FFT")?;
    INVERSE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        let mut spectrum = filled(len, Complex::new(0.0, 0.0), "cosine inverse FFT workspace")?;
        for a in 0..d {
            spectrum.fill(Complex::new(0.0, 0.0));
            spectrum[0] = Complex::new(coeff[a], 0.0);
            for k in 1..k_count {
                let half = Complex::new(0.5 * coeff[k * d + a], 0.0);
                spectrum[k] = half;
                spectrum[len - k] = half;
            }
            plan.process(&mut spectrum);
            for t in 0..=m {
                out[t * d + a] = spectrum[t].re;
            }
        }
        Ok::<(), SpectralError>(())
    })?;
    Ok(out)
}

/// Weighted sine projection on an endpoint-halved grid.
pub fn dst_project(
    samples: &[f64],
    m: usize,
    k_count: usize,
    d: usize,
    beta: f64,
) -> Result<Vec<f64>, SpectralError> {
    validate_transform(samples.len(), m, k_count, d, "sine samples")?;
    if !(beta > 0.0 && beta.is_finite()) {
        return Err(SpectralError::Positive {
            name: "sine projection beta",
            value: beta,
        });
    }
    let output_len = checked_product(k_count, d, "sine projection")?;
    let mut out = filled(output_len, 0.0, "sine projection")?;
    let scale = -beta / (2.0 * m as f64);
    validate_finite(scale, "sine projection scale", 0)?;
    let wave_len = checked_product(2, m, "sine FFT")?;
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(wave_len);
        let mut wave = filled(wave_len, Complex::new(0.0, 0.0), "sine FFT workspace")?;
        for a in 0..d {
            wave.fill(Complex::new(0.0, 0.0));
            for t in 1..m {
                let value = samples[t * d + a];
                wave[t] = Complex::new(value, 0.0);
                wave[wave_len - t] = Complex::new(-value, 0.0);
            }
            plan.process(&mut wave);
            for k in 1..k_count {
                out[k * d + a] = wave[k].im * scale;
            }
        }
        Ok::<(), SpectralError>(())
    })?;
    Ok(out)
}

/// Sine-series synthesis on an endpoint-halved grid.
pub fn dst_synthesize(
    coeff: &[f64],
    k_count: usize,
    d: usize,
    m: usize,
) -> Result<Vec<f64>, SpectralError> {
    let nodes = validate_coefficients(coeff.len(), m, k_count, d, "sine coefficients")?;
    if m < 2 {
        return direct_sine_synthesize(coeff, k_count, d, m);
    }
    let output_len = checked_product(nodes, d, "sine synthesis")?;
    let mut out = filled(output_len, 0.0, "sine synthesis")?;
    let len = checked_product(2, m, "sine inverse FFT")?;
    INVERSE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        let mut spectrum = filled(len, Complex::new(0.0, 0.0), "sine inverse FFT workspace")?;
        for a in 0..d {
            spectrum.fill(Complex::new(0.0, 0.0));
            for k in 1..k_count {
                let half = 0.5 * coeff[k * d + a];
                spectrum[k] = Complex::new(0.0, -half);
                spectrum[len - k] = Complex::new(0.0, half);
            }
            plan.process(&mut spectrum);
            for t in 1..m {
                out[t * d + a] = spectrum[t].re;
            }
        }
        Ok::<(), SpectralError>(())
    })?;
    Ok(out)
}

/// Direct weighted cosine projection, provided as a small-grid reference.
pub fn direct_project(
    samples: &[f64],
    m: usize,
    k_count: usize,
    d: usize,
    beta: f64,
) -> Result<Vec<f64>, SpectralError> {
    validate_transform(samples.len(), m, k_count, d, "direct cosine samples")?;
    let (theta, _, weights) = half_period_grid(m, beta)?;
    let output_len = checked_product(k_count, d, "direct cosine projection")?;
    let mut out = filled(output_len, 0.0, "direct cosine projection")?;
    for k in 0..k_count {
        for t in 0..=m {
            let weighted_cosine = weights[t] * (theta[t] * k as f64).cos();
            for a in 0..d {
                out[k * d + a] += weighted_cosine * samples[t * d + a];
            }
        }
    }
    Ok(out)
}

/// Kahan-summed midpoint rule on `[0, upper]`.
///
/// The callback receives each panel midpoint. Its error type is preserved.
pub fn midpoint_integral<E, F>(
    upper: f64,
    panels: usize,
    mut value: F,
) -> Result<f64, EvaluationError<E>>
where
    F: FnMut(f64) -> Result<f64, E>,
{
    if !(upper > 0.0 && upper.is_finite()) {
        return Err(SpectralError::Positive {
            name: "midpoint upper bound",
            value: upper,
        }
        .into());
    }
    checked_nodes(panels, 1, "midpoint rule")?;
    let step = upper / panels as f64;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let point = (panel as f64 + 0.5) * step;
        let term = value(point).map_err(EvaluationError::Evaluator)?;
        validate_finite(term, "midpoint sample", panel)?;
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    let integral = step * sum;
    validate_finite(integral, "midpoint integral", 0)?;
    Ok(integral)
}

/// Endpoint enclosure for `integral_0^infinity f(s) exp(-s) ds`.
#[derive(Clone, Copy, Debug)]
pub struct MonotoneExponentialEnclosure {
    pub lower: f64,
    pub upper: f64,
    pub monotone: bool,
}

/// Enclose a nonnegative, nondecreasing factor under an exponential weight.
///
/// The upper result includes the conservative tail `exp(-upper)` obtained
/// from the additional assumption `f <= 1`.
pub fn monotone_exponential_enclosure<E, F>(
    upper: f64,
    panels: usize,
    mut value: F,
) -> Result<MonotoneExponentialEnclosure, EvaluationError<E>>
where
    F: FnMut(f64) -> Result<f64, E>,
{
    if !(upper > 0.0 && upper.is_finite()) {
        return Err(SpectralError::Positive {
            name: "monotone enclosure upper bound",
            value: upper,
        }
        .into());
    }
    checked_nodes(panels, 1, "monotone enclosure")?;
    let step = upper / panels as f64;
    let panel_weight_factor = -(-step).exp_m1();
    let mut lower = 0.0;
    let mut lower_correction = 0.0;
    let mut upper_sum = 0.0;
    let mut upper_correction = 0.0;
    let mut previous = 0.0;
    let mut monotone = true;
    for panel in 0..panels {
        let left = panel as f64 * step;
        let right = (panel + 1) as f64 * step;
        let current = value(right).map_err(EvaluationError::Evaluator)?;
        validate_finite(current, "monotone enclosure sample", panel)?;
        if !(0.0..=1.0).contains(&current) {
            return Err(SpectralError::OutOfRange {
                name: "monotone enclosure sample",
                index: panel,
                value: current,
            }
            .into());
        }
        if current < previous {
            monotone = false;
            break;
        }
        let weight = (-left).exp() * panel_weight_factor;
        kahan_add(&mut lower, &mut lower_correction, previous * weight);
        kahan_add(&mut upper_sum, &mut upper_correction, current * weight);
        previous = current;
    }
    let upper_sum = upper_sum + (-upper).exp();
    validate_finite(lower, "monotone enclosure lower bound", 0)?;
    validate_finite(upper_sum, "monotone enclosure upper bound", 0)?;
    Ok(MonotoneExponentialEnclosure {
        lower,
        upper: upper_sum,
        monotone,
    })
}

/// Midpoint lower sum and trapezoid upper sum over segmented uniform grids.
///
/// `exact` converts binary64 nodes and constants into the caller's arithmetic
/// type, allowing interval and extended-precision callers without coupling
/// this crate to one arithmetic implementation.
pub fn midpoint_trapezoid_enclosure<T, L, U, X>(
    segments: &[(f64, f64, usize)],
    mut lower_value: L,
    mut upper_value: U,
    exact: X,
) -> Result<(T, T, usize), SpectralError>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
    L: FnMut(T) -> T,
    U: FnMut(T) -> T,
    X: Fn(f64) -> T,
{
    let mut lower = exact(0.0);
    let mut upper = exact(0.0);
    let mut total_panels: usize = 0;
    for &(start, end, count) in segments {
        validate_interval(start, end, "midpoint-trapezoid segment")?;
        checked_nodes(count, 1, "midpoint-trapezoid segment")?;
        total_panels = total_panels
            .checked_add(count)
            .ok_or(SpectralError::SizeOverflow {
                name: "midpoint-trapezoid panels",
            })?;
        let width = (end - start) / count as f64;
        let mut left = start;
        let mut left_value = upper_value(exact(left));
        for _ in 0..count {
            let right = left + width;
            let left_exact = exact(left);
            let right_exact = exact(right);
            let step = right_exact - left_exact;
            let right_value = upper_value(right_exact);
            upper = upper + step * (left_value + right_value) / exact(2.0);
            let middle = (left_exact + right_exact) / exact(2.0);
            lower = lower + step * lower_value(middle);
            left = right;
            left_value = right_value;
        }
    }
    Ok((lower, upper, total_panels))
}

fn kahan_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

/// Composite trapezoid in log space over equally refined segments.
///
/// The endpoint weight is applied independently on every segment, so a shared
/// endpoint contributes once from each adjacent segment.
pub fn log_composite_trapezoid<E, F>(
    segments: &[(f64, f64)],
    panels: usize,
    mut log_value: F,
) -> Result<Option<f64>, EvaluationError<E>>
where
    F: FnMut(f64) -> Result<Option<f64>, E>,
{
    log_composite_trapezoid_by(segments, panels, |point, log_spacing, log_weight| {
        log_value(point).map(|term| term.map(|value| log_spacing + log_weight + value))
    })
}

/// Composite trapezoid whose callback controls the floating-point association.
///
/// The callback receives the node, `ln(spacing)`, and `ln(endpoint weight)`.
/// This form lets callers preserve a reference expression's exact operation
/// order while sharing validation, storage, and log-sum-exp reduction.
pub fn log_composite_trapezoid_by<E, F>(
    segments: &[(f64, f64)],
    panels: usize,
    mut log_term: F,
) -> Result<Option<f64>, EvaluationError<E>>
where
    F: FnMut(f64, f64, f64) -> Result<Option<f64>, E>,
{
    let nodes = checked_nodes(panels, 1, "log composite trapezoid")?;
    let capacity = checked_product(segments.len(), nodes, "log composite terms")?;
    let mut terms = Vec::new();
    terms
        .try_reserve_exact(capacity)
        .map_err(|_| SpectralError::Allocation {
            name: "log composite terms",
            size: capacity,
        })?;
    for &(lower, upper) in segments {
        validate_interval(lower, upper, "log composite segment")?;
        let spacing = (upper - lower) / panels as f64;
        let log_spacing = spacing.ln();
        for index in 0..=panels {
            let point = lower + spacing * index as f64;
            let weight: f64 = if index == 0 || index == panels {
                0.5
            } else {
                1.0
            };
            if let Some(term) =
                log_term(point, log_spacing, weight.ln()).map_err(EvaluationError::Evaluator)?
            {
                validate_finite(term, "log composite term", terms.len())?;
                terms.push(term);
            }
        }
    }
    let Some(largest) = terms.iter().copied().max_by(f64::total_cmp) else {
        return Ok(None);
    };
    let scaled_sum: f64 = terms.iter().map(|term| (*term - largest).exp()).sum();
    let result = largest + scaled_sum.ln();
    validate_finite(result, "log composite result", 0)?;
    Ok(Some(result))
}

/// Composite trapezoid for a positive log integrand.
///
/// Each panel is scaled by its larger endpoint before exponentiation.
pub fn log_stable_trapezoid<F>(
    lower: f64,
    upper: f64,
    panels: usize,
    mut log_value: F,
) -> Result<f64, SpectralError>
where
    F: FnMut(f64) -> f64,
{
    validate_interval(lower, upper, "log-stable trapezoid")?;
    checked_nodes(panels, 1, "log-stable trapezoid")?;
    let mut sum = 0.0;
    let mut previous = log_value(lower);
    validate_finite(previous, "log-stable sample", 0)?;
    for index in 0..panels {
        let left = lower + (upper - lower) * index as f64 / panels as f64;
        let right = lower + (upper - lower) * (index + 1) as f64 / panels as f64;
        let next = log_value(right);
        validate_finite(next, "log-stable sample", index + 1)?;
        let width = right - left;
        // `top` is the exact log-sum-exp scale bound for both endpoints, not
        // a floor: subtracting it makes each exponential no larger than one.
        let top = if previous >= next { previous } else { next };
        sum += width * top.exp() * 0.5 * ((previous - top).exp() + (next - top).exp());
        previous = next;
    }
    validate_finite(sum, "log-stable result", 0)?;
    Ok(sum)
}

/// One sample for piecewise-linear exponential quadrature.
#[derive(Clone, Copy, Debug)]
pub struct ExponentialSample {
    pub coordinate: f64,
    pub exponent: f64,
    pub amplitude: f64,
}

/// Trapezoid and exact integral of piecewise-linear amplitude and exponent.
pub fn piecewise_exponential_integral(
    samples: &[ExponentialSample],
) -> Result<(f64, f64), SpectralError> {
    piecewise_exponential_integral_by(samples, |sample| *sample)
}

/// Generic form of [`piecewise_exponential_integral`] without a copy.
pub fn piecewise_exponential_integral_by<T, F>(
    samples: &[T],
    sample: F,
) -> Result<(f64, f64), SpectralError>
where
    F: Fn(&T) -> ExponentialSample,
{
    let mut trapezoid = 0.0;
    let mut linear = 0.0;
    for index in 0..samples.len().saturating_sub(1) {
        let left = sample(&samples[index]);
        let right = sample(&samples[index + 1]);
        validate_exponential_sample(left, index)?;
        validate_exponential_sample(right, index + 1)?;
        let width = right.coordinate - left.coordinate;
        if !(width > 0.0 && width.is_finite()) {
            return Err(SpectralError::Interval {
                name: "piecewise exponential panel",
                lower: left.coordinate,
                upper: right.coordinate,
            });
        }
        let amplitude0 = left.amplitude;
        let amplitude1 = right.amplitude;
        let exponent0 = left.exponent;
        let exponent1 = right.exponent;
        let value0 = (-exponent0).exp();
        let value1 = (-exponent1).exp();
        trapezoid += 0.5 * (amplitude0 * value0 + amplitude1 * value1) * width;
        linear += piecewise_exponential_panel(width, amplitude0, amplitude1, exponent0, exponent1)?;
    }
    validate_finite(trapezoid, "piecewise trapezoid", 0)?;
    validate_finite(linear, "piecewise linear integral", 0)?;
    Ok((trapezoid, linear))
}

fn validate_exponential_sample(
    sample: ExponentialSample,
    index: usize,
) -> Result<(), SpectralError> {
    validate_finite(sample.coordinate, "exponential coordinate", index)?;
    validate_finite(sample.exponent, "exponential exponent", index)?;
    validate_finite(sample.amplitude, "exponential amplitude", index)
}

/// Exact panel integral for linear amplitude and exponent interpolants.
pub fn piecewise_exponential_panel(
    width: f64,
    amplitude0: f64,
    amplitude1: f64,
    exponent0: f64,
    exponent1: f64,
) -> Result<f64, SpectralError> {
    if !(width > 0.0 && width.is_finite()) {
        return Err(SpectralError::Positive {
            name: "piecewise panel width",
            value: width,
        });
    }
    for (index, value) in [amplitude0, amplitude1, exponent0, exponent1]
        .into_iter()
        .enumerate()
    {
        validate_finite(value, "piecewise panel argument", index)?;
    }
    let exponent_delta = exponent1 - exponent0;
    let amplitude_delta = amplitude1 - amplitude0;
    let value0 = (-exponent0).exp();
    let integral = if exponent_delta.abs() < 1.0e-10 {
        width * value0 * (amplitude0 + 0.5 * amplitude_delta)
    } else {
        let decay = (-exponent_delta).exp();
        let constant_integral = (1.0 - decay) / exponent_delta;
        let linear_integral =
            (1.0 - decay * (1.0 + exponent_delta)) / (exponent_delta * exponent_delta);
        width * value0 * (amplitude0 * constant_integral + amplitude_delta * linear_integral)
    };
    validate_finite(integral, "piecewise panel integral", 0)?;
    Ok(integral)
}

/// One fixed-degree Clenshaw--Curtis panel and its interpolation data.
#[derive(Clone, Debug)]
pub struct ClenshawCurtisPanel {
    pub value: f64,
    pub rounding: f64,
    pub monotone: bool,
    pub antiderivative: Vec<f64>,
    pub at_plus: f64,
}

/// Clenshaw--Curtis interpolation and integration on `[lower, upper]`.
///
/// The callback returns the sampled value and an absolute sample error.
pub fn clenshaw_curtis_panel<E, F>(
    lower: f64,
    upper: f64,
    degree: usize,
    mut sample: F,
) -> Result<ClenshawCurtisPanel, EvaluationError<E>>
where
    F: FnMut(f64) -> Result<(f64, f64), E>,
{
    validate_interval(lower, upper, "Clenshaw-Curtis panel")?;
    let nodes = checked_nodes(degree, 2, "Clenshaw-Curtis panel")?;
    let middle = 0.5 * (lower + upper);
    let half = 0.5 * (upper - lower);
    let mut samples = filled(nodes, 0.0, "Clenshaw-Curtis samples")?;
    let mut errors = filled(nodes, 0.0, "Clenshaw-Curtis sample errors")?;
    for index in 0..=degree {
        let coordinate = (PI * index as f64 / degree as f64).cos();
        let (value, error) =
            sample(middle + half * coordinate).map_err(EvaluationError::Evaluator)?;
        validate_finite(value, "Clenshaw-Curtis sample", index)?;
        validate_finite(error, "Clenshaw-Curtis sample error", index)?;
        if error < 0.0 {
            return Err(SpectralError::OutOfRange {
                name: "Clenshaw-Curtis sample error",
                index,
                value: error,
            }
            .into());
        }
        samples[index] = value;
        errors[index] = error;
    }
    let mut coefficients = filled(nodes, 0.0, "Clenshaw-Curtis coefficients")?;
    let mut rounding = 0.0_f64;
    for (mode, coefficient) in coefficients.iter_mut().enumerate() {
        let mut sum = 0.5 * samples[0];
        let mut error = 0.5 * errors[0];
        for index in 1..degree {
            let cosine = (PI * (index * mode) as f64 / degree as f64).cos();
            sum += samples[index] * cosine;
            error += errors[index] * cosine.abs();
        }
        let endpoint = if mode % 2 == 0 { 0.5 } else { -0.5 };
        sum += endpoint * samples[degree];
        error += 0.5 * errors[degree];
        let scale = if mode == 0 || mode == degree {
            1.0 / degree as f64
        } else {
            2.0 / degree as f64
        };
        *coefficient = scale * sum;
        if mode % 2 == 0 {
            let integral = if mode == 0 {
                2.0
            } else {
                2.0 / (1.0 - (mode * mode) as f64)
            };
            rounding += half * scale * error * integral.abs();
        }
    }
    let antiderivative = chebyshev_antiderivative(&coefficients);
    let at_plus = clenshaw(&antiderivative, 1.0);
    let value = half * (at_plus - clenshaw(&antiderivative, -1.0));
    validate_finite(value, "Clenshaw-Curtis integral", 0)?;
    validate_finite(rounding, "Clenshaw-Curtis rounding", 0)?;
    // Monotonicity is pairwise only when samples are exact. With the certified
    // accumulated error `rounding`, compare every partial integral against the
    // greatest earlier partial integral: if any earlier value exceeds it by
    // more than that bound, no monotone function is consistent with the data.
    // `previous` is therefore a history statistic, not a numerical floor.
    let mut previous = 0.0_f64;
    let mut monotone = true;
    for index in 0..=degree {
        let coordinate = (PI * index as f64 / degree as f64).cos();
        let above = half * (at_plus - clenshaw(&antiderivative, coordinate));
        validate_finite(above, "Clenshaw-Curtis partial integral", index)?;
        if above + rounding < previous {
            monotone = false;
            break;
        }
        if above > previous {
            previous = above;
        }
    }
    Ok(ClenshawCurtisPanel {
        value,
        rounding,
        monotone,
        antiderivative,
        at_plus,
    })
}

/// Evaluate a Chebyshev series by Clenshaw recurrence.
pub fn clenshaw(coefficients: &[f64], coordinate: f64) -> f64 {
    clenshaw_bound(coefficients, coordinate).0
}

/// Chebyshev-series value and a forward roundoff estimate.
pub fn clenshaw_bound(coefficients: &[f64], coordinate: f64) -> (f64, f64) {
    let mut next = 0.0;
    let mut current = 0.0;
    let mut l1 = 0.0;
    for &term in coefficients.iter().rev() {
        l1 += term.abs();
        let updated = term + 2.0 * coordinate * current - next;
        next = current;
        current = updated;
    }
    let value = current - coordinate * next;
    let operations = coefficients.len() as f64 * 3.0 + 2.0;
    // `l1` bounds the exact series on [-1, 1], while `abs(value)` covers
    // evaluation outside that interval and accumulated rounding. Their larger
    // value is therefore the derived forward-error scale, not a value floor.
    let scale = if l1 >= value.abs() { l1 } else { value.abs() };
    (value, operations * f64::EPSILON * scale)
}

/// Coefficients of an indefinite Chebyshev antiderivative.
pub fn chebyshev_antiderivative(coefficients: &[f64]) -> Vec<f64> {
    let count = coefficients.len();
    let mut integrated = vec![0.0; count + 2];
    let coefficient = |index: usize| {
        if index < count {
            coefficients[index]
        } else {
            0.0
        }
    };
    if count > 0 {
        integrated[1] = coefficient(0) - coefficient(2) / 2.0;
    }
    for (mode, value) in integrated.iter_mut().enumerate().skip(2) {
        *value = (coefficient(mode - 1) - coefficient(mode + 1)) / (2.0 * mode as f64);
    }
    integrated
}

/// Direct cosine-series synthesis, provided as a small-grid reference.
pub fn direct_cosine_synthesize(
    coeff: &[f64],
    k_count: usize,
    d: usize,
    m: usize,
) -> Result<Vec<f64>, SpectralError> {
    direct_synthesize(coeff, k_count, d, m, false)
}

/// Direct sine-series synthesis, provided as a small-grid reference.
pub fn direct_sine_synthesize(
    coeff: &[f64],
    k_count: usize,
    d: usize,
    m: usize,
) -> Result<Vec<f64>, SpectralError> {
    direct_synthesize(coeff, k_count, d, m, true)
}

fn direct_synthesize(
    coeff: &[f64],
    k_count: usize,
    d: usize,
    m: usize,
    sine: bool,
) -> Result<Vec<f64>, SpectralError> {
    let nodes = validate_coefficients(coeff.len(), m, k_count, d, "direct coefficients")?;
    let (theta, _, _) = half_period_grid(m, 1.0)?;
    let output_len = checked_product(nodes, d, "direct synthesis")?;
    let mut out = filled(output_len, 0.0, "direct synthesis")?;
    for t in 0..=m {
        for k in 0..k_count {
            let phase = theta[t] * k as f64;
            let basis = if sine { phase.sin() } else { phase.cos() };
            for a in 0..d {
                out[t * d + a] += basis * coeff[k * d + a];
            }
        }
    }
    Ok(out)
}

/// `sinhc h = sinh(h)/h` and `E` with `d sinhc/dh = h E`, from their series.
///
/// Every term is nonnegative. The finite iteration limit is larger than the
/// binary64 exponent range, so a finite series must settle or underflow before
/// the limit; a non-finite intermediate is refused immediately.
fn sinhc_pair(h: f64) -> Result<(f64, f64), SpectralError> {
    validate_finite(h, "sinhc argument", 0)?;
    let h2 = h * h;
    validate_finite(h2, "sinhc squared argument", 0)?;
    let mut value = 1.0;
    let mut slope = 0.0;
    let mut e = 1.0 / 6.0;
    let mut n = 1.0_f64;
    for _ in 0..4096 {
        let term = e * h2;
        let slope_term = 2.0 * n * e;
        validate_finite(term, "sinhc term", n as usize)?;
        validate_finite(slope_term, "sinhc slope term", n as usize)?;
        let settled = term <= f64::EPSILON * value && slope_term <= f64::EPSILON * slope;
        value += term;
        slope += slope_term;
        validate_finite(value, "sinhc value", n as usize)?;
        validate_finite(slope, "sinhc slope", n as usize)?;
        if settled {
            return Ok((value, slope));
        }
        n += 1.0;
        e *= h2 / ((2.0 * n) * (2.0 * n + 1.0));
    }
    Err(SpectralError::NonConverged { name: "sinhc" })
}

/// `asinh(hi) - asinh(lo)` from the exact `gap = hi - lo`.
///
/// The `ln_1p` form keeps relative precision when the gap is small.
fn asinh_gap(hi: f64, lo: f64, gap: f64) -> Result<f64, SpectralError> {
    for (index, value) in [hi, lo, gap].into_iter().enumerate() {
        validate_finite(value, "asinh-gap argument", index)?;
    }
    let roots = (1.0 + hi * hi).sqrt() + (1.0 + lo * lo).sqrt();
    let bracket = 1.0 + (hi + lo) / roots;
    let result = (gap * bracket * (-lo.asinh()).exp()).ln_1p();
    validate_finite(result, "asinh gap", 0)?;
    Ok(result)
}

/// The chirp `theta -> phi`, fixed in units of `theta` for every period.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpMap {
    /// Kink centre `theta*` inside `(0, pi)`.
    pub theta_star: f64,
    /// Kink half-width in `u = cos(theta)`.
    pub sigma: f64,
    /// `asinh((1 - u*) / sigma)`.
    s_plus: f64,
    /// `asinh((-1 - u*) / sigma)`, which is negative.
    s_minus: f64,
}

/// One map node: `theta(phi)`, `d theta / d phi`, and its derivative.
#[derive(Debug, Clone, Copy)]
pub struct MapPoint {
    /// Physical half-period angle.
    pub theta: f64,
    /// `d theta / d phi`.
    pub g: f64,
    /// `d² theta / d phi²`.
    pub g_phi: f64,
}

impl ChirpMap {
    /// Construct a finite monotone map centred at `theta_star`.
    pub fn new(theta_star: f64, sigma: f64) -> Result<Self, SpectralError> {
        if !(theta_star > 0.0 && theta_star < PI && theta_star.is_finite()) {
            return Err(SpectralError::ChirpCentre { value: theta_star });
        }
        if !(sigma > 0.0 && sigma.is_finite()) {
            return Err(SpectralError::ChirpWidth { value: sigma });
        }
        let half = 0.5 * theta_star;
        let above = 2.0 * half.sin().powi(2) / sigma;
        let below = 2.0 * half.cos().powi(2) / sigma;
        if !(above.is_finite() && below.is_finite()) {
            return Err(SpectralError::ChirpWidth { value: sigma });
        }
        let s_plus = above.asinh();
        let s_minus = -below.asinh();
        let span = s_plus - s_minus;
        if !(span > 0.0 && span.is_finite()) {
            return Err(SpectralError::ChirpWidth { value: sigma });
        }
        Ok(Self {
            theta_star,
            sigma,
            s_plus,
            s_minus,
        })
    }

    /// Finite map span `s_plus - s_minus`.
    pub fn span(&self) -> f64 {
        self.s_plus - self.s_minus
    }

    /// Evaluate `theta`, `d theta/d phi`, and `d² theta/d phi²`.
    ///
    /// The half-angle form cancels the common endpoint factor analytically,
    /// retaining relative precision next to `phi = 0` and `phi = pi`.
    pub fn inverse(&self, phi: f64) -> Result<MapPoint, SpectralError> {
        if !(phi.is_finite() && (0.0..=PI).contains(&phi)) {
            return Err(SpectralError::MapCoordinate {
                name: "phi",
                value: phi,
            });
        }
        let span = self.span();
        let q_plus = (0.5 * phi).sin().powi(2);
        let q_minus = (0.5 * phi).cos().powi(2);
        let h_plus = 0.5 * span * q_plus;
        let h_minus = 0.5 * span * q_minus;
        let s = self.s_plus - 2.0 * h_plus;
        let a_plus = self.s_plus - h_plus;
        let a_minus = self.s_minus + h_minus;
        let one_minus_u = 2.0 * self.sigma * a_plus.cosh() * h_plus.sinh();
        let one_plus_u = 2.0 * self.sigma * a_minus.cosh() * h_minus.sinh();
        let theta = 2.0 * one_minus_u.sqrt().atan2(one_plus_u.sqrt());
        let (c_plus, e_plus) = sinhc_pair(h_plus)?;
        let (c_minus, e_minus) = sinhc_pair(h_minus)?;
        let product = a_plus.cosh() * c_plus * a_minus.cosh() * c_minus;
        let g = s.cosh() / product.sqrt();
        let sin_phi = phi.sin();
        let ds = -0.5 * span * sin_phi;
        let dh_plus = 0.25 * span * sin_phi;
        let dh_minus = -dh_plus;
        let dlog = s.tanh() * ds
            - 0.5
                * ((h_plus * e_plus / c_plus - a_plus.tanh()) * dh_plus
                    + (h_minus * e_minus / c_minus + a_minus.tanh()) * dh_minus);
        let point = MapPoint {
            theta,
            g,
            g_phi: g * dlog,
        };
        for (index, value) in [point.theta, point.g, point.g_phi].into_iter().enumerate() {
            validate_finite(value, "inverse chirp result", index)?;
        }
        Ok(point)
    }

    /// Evaluate `phi(theta)` using cancellation-free inverse-hyperbolic gaps.
    pub fn forward(&self, theta: f64) -> Result<f64, SpectralError> {
        if !theta.is_finite() {
            return Err(SpectralError::MapCoordinate {
                name: "theta",
                value: theta,
            });
        }
        // A folded half-period coordinate is formed by a multiply and a
        // subtraction from `2*pi`; their two rounded operations can place an
        // exact endpoint at most `2 eps*pi` outside the closed interval.
        let endpoint_roundoff = 2.0 * f64::EPSILON * PI;
        if theta < 0.0 {
            if theta >= -endpoint_roundoff {
                return Ok(0.0);
            }
            return Err(SpectralError::MapCoordinate {
                name: "theta",
                value: theta,
            });
        }
        if theta > PI {
            if theta <= PI + endpoint_roundoff {
                return Ok(PI);
            }
            return Err(SpectralError::MapCoordinate {
                name: "theta",
                value: theta,
            });
        }
        let above = 2.0 * (0.5 * theta).sin().powi(2) / self.sigma;
        let below = 2.0 * (0.5 * theta).cos().powi(2) / self.sigma;
        let a = self.s_plus.sinh();
        let c = -self.s_minus.sinh();
        let b = a - above;
        let upper = asinh_gap(a, b, above)?;
        let lower = asinh_gap(c, -b, below)?;
        let phi = 2.0 * upper.sqrt().atan2(lower.sqrt());
        validate_finite(phi, "forward chirp result", 0)?;
        Ok(phi)
    }
}

/// Uniform phase nodes with physical times and map rates.
#[derive(Debug, Clone)]
pub struct ChirpGrid {
    /// Number of panels.
    pub m: usize,
    /// Imaginary-time period.
    pub beta: f64,
    /// Coordinate map shared by all nodes.
    pub map: ChirpMap,
    /// Uniform phase nodes.
    pub phi: Vec<f64>,
    /// Physical imaginary-time nodes.
    pub tau: Vec<f64>,
    /// `d phi / d tau` at each node.
    pub rate: Vec<f64>,
    /// Endpoint-halved phase weights.
    pub weights: Vec<f64>,
}

impl ChirpGrid {
    /// Build a checked endpoint grid for a finite positive period.
    pub fn new(map: ChirpMap, beta: f64, m: usize) -> Result<Self, SpectralError> {
        let nodes = checked_nodes(m, 2, "chirp grid")?;
        if !(beta > 0.0 && beta.is_finite()) {
            return Err(SpectralError::Positive {
                name: "chirp-grid beta",
                value: beta,
            });
        }
        let omega = 2.0 * PI / beta;
        let mut phi = filled(nodes, 0.0, "chirp phase")?;
        let mut tau = filled(nodes, 0.0, "chirp time")?;
        let mut rate = filled(nodes, 0.0, "chirp rate")?;
        let mut weights = filled(nodes, PI / m as f64, "chirp weights")?;
        weights[0] *= 0.5;
        weights[m] *= 0.5;
        for t in 0..=m {
            let phase = PI * t as f64 / m as f64;
            let point = map.inverse(phase)?;
            let mapped_rate = omega / point.g;
            validate_finite(mapped_rate, "chirp rate", t)?;
            phi[t] = phase;
            tau[t] = point.theta / omega;
            rate[t] = mapped_rate;
        }
        Ok(Self {
            m,
            beta,
            map,
            phi,
            tau,
            rate,
            weights,
        })
    }

    /// Number of nodes, including both endpoints.
    pub fn nodes(&self) -> usize {
        self.m + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(condition: bool, message: &'static str) -> Result<(), String> {
        if condition {
            Ok(())
        } else {
            Err(message.to_string())
        }
    }

    fn input_error<E>(error: EvaluationError<E>) -> Result<SpectralError, String> {
        match error {
            EvaluationError::Input(input) => Ok(input),
            EvaluationError::Evaluator(_) => Err("unexpected evaluator error".to_string()),
        }
    }

    #[test]
    fn half_period_grid_checks_size_and_endpoints() -> Result<(), String> {
        let (theta, tau, weights) = half_period_grid(4, 2.0).map_err(|e| e.to_string())?;
        check(
            theta.len() == 5 && tau.len() == 5 && weights.len() == 5,
            "grid length",
        )?;
        check(theta[0].to_bits() == 0.0_f64.to_bits(), "theta start")?;
        check(theta[4].to_bits() == PI.to_bits(), "theta end")?;
        check(tau[4].to_bits() == 1.0_f64.to_bits(), "tau end")?;
        check(
            weights[0].to_bits() == 0.25_f64.to_bits(),
            "endpoint weight",
        )?;
        match half_period_grid(usize::MAX, 1.0) {
            Err(SpectralError::SizeOverflow { .. }) => Ok(()),
            _ => Err("huge half-period grid was not refused".to_string()),
        }
    }

    #[test]
    fn gauss_legendre_two_point_integrates_a_quadratic() -> Result<(), String> {
        let rule = gauss_legendre_two_point(0.0, 2.0).map_err(|e| e.to_string())?;
        let integral: f64 = rule
            .iter()
            .map(|(point, weight)| point * point * weight)
            .sum();
        check((integral - 2.0 / 3.0).abs() < 1.0e-15, "two-point rule")?;
        match gauss_legendre_two_point(0.0, f64::NAN) {
            Err(SpectralError::Positive { .. }) => Ok(()),
            _ => Err("non-finite Gauss-Legendre width was not refused".to_string()),
        }
    }

    #[test]
    fn endpoint_trapezoid_array_checks_panels_and_samples() -> Result<(), String> {
        let value = endpoint_trapezoid_array::<2, _>(0.0, 1.0, 4, |x| [x, 1.0])
            .map_err(|e| e.to_string())?;
        check(value[0].to_bits() == 0.5_f64.to_bits(), "linear trapezoid")?;
        check(
            value[1].to_bits() == 1.0_f64.to_bits(),
            "constant trapezoid",
        )?;
        match endpoint_trapezoid_array::<1, _>(0.0, 1.0, 0, |x| [x]) {
            Err(SpectralError::PanelCount { .. }) => Ok(()),
            _ => Err("zero-panel endpoint trapezoid was not refused".to_string()),
        }
    }

    #[test]
    fn forward_inplace_checks_empty_and_transforms() -> Result<(), String> {
        let mut values = [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
        forward_inplace(&mut values).map_err(|e| e.to_string())?;
        check(values[0].re.to_bits() == 1.0_f64.to_bits(), "FFT zero mode")?;
        check(
            values[1].re.to_bits() == 1.0_f64.to_bits(),
            "FFT Nyquist mode",
        )?;
        let mut empty = [];
        match forward_inplace(&mut empty) {
            Err(SpectralError::PanelCount { .. }) => Ok(()),
            _ => Err("empty FFT was not refused".to_string()),
        }
    }

    #[test]
    fn cc_project_matches_direct_and_checks_shape() -> Result<(), String> {
        let m = 32;
        let modes = 16;
        let dimensions = 3;
        let beta = 2.5;
        let samples: Vec<f64> = (0..((m + 1) * dimensions))
            .map(|index| (index as f64 * 0.17).sin())
            .collect();
        let fft = cc_project(&samples, m, modes, dimensions, beta).map_err(|e| e.to_string())?;
        let direct =
            direct_project(&samples, m, modes, dimensions, beta).map_err(|e| e.to_string())?;
        for (left, right) in fft.iter().zip(&direct) {
            if (left - right).abs() >= 1.0e-12 {
                return Err("FFT and direct cosine projections differ".to_string());
            }
        }
        match cc_project(&samples[..3], m, modes, dimensions, beta) {
            Err(SpectralError::Shape { .. }) => Ok(()),
            _ => Err("short cosine samples were not refused".to_string()),
        }
    }

    #[test]
    fn cc_synthesize_matches_direct_and_checks_modes() -> Result<(), String> {
        let coefficients = [0.4, -0.2, 0.1, 0.03];
        let fft = cc_synthesize(&coefficients, 4, 1, 8).map_err(|e| e.to_string())?;
        let direct = direct_cosine_synthesize(&coefficients, 4, 1, 8).map_err(|e| e.to_string())?;
        for (left, right) in fft.iter().zip(&direct) {
            check((left - right).abs() < 1.0e-14, "cosine synthesis mismatch")?;
        }
        match cc_synthesize(&coefficients, 10, 1, 8) {
            Err(SpectralError::ModeCount { .. }) => Ok(()),
            _ => Err("excess cosine modes were not refused".to_string()),
        }
    }

    #[test]
    fn dst_project_matches_direct_reference_and_checks_modes() -> Result<(), String> {
        let m = 8;
        let samples: Vec<f64> = (0..=m)
            .map(|index| (PI * index as f64 / m as f64).sin())
            .collect();
        let projected = dst_project(&samples, m, 5, 1, 2.0).map_err(|e| e.to_string())?;
        check((projected[1] - 1.0).abs() < 1.0e-14, "sine projection")?;
        match dst_project(&samples, m, 10, 1, 2.0) {
            Err(SpectralError::ModeCount { .. }) => Ok(()),
            _ => Err("excess sine projection modes were not refused".to_string()),
        }
    }

    #[test]
    fn dst_synthesize_matches_direct_and_checks_shape() -> Result<(), String> {
        let coefficients = [0.0, 1.0, -0.2, 0.05];
        let fft = dst_synthesize(&coefficients, 4, 1, 8).map_err(|e| e.to_string())?;
        let direct = direct_sine_synthesize(&coefficients, 4, 1, 8).map_err(|e| e.to_string())?;
        for (left, right) in fft.iter().zip(&direct) {
            check((left - right).abs() < 1.0e-14, "sine synthesis mismatch")?;
        }
        match dst_synthesize(&coefficients[..2], 4, 1, 8) {
            Err(SpectralError::Shape { .. }) => Ok(()),
            _ => Err("short sine coefficients were not refused".to_string()),
        }
    }

    #[test]
    fn direct_project_checks_dimension() -> Result<(), String> {
        match direct_project(&[], 4, 2, 0, 1.0) {
            Err(SpectralError::Dimension { .. }) => Ok(()),
            _ => Err("zero-dimensional projection was not refused".to_string()),
        }
    }

    #[test]
    fn midpoint_integral_preserves_reference_and_refuses_zero_panels() -> Result<(), String> {
        let midpoint = midpoint_integral(1.0, 64, |point| Ok::<f64, SpectralError>(point * point))
            .map_err(|e| e.to_string())?;
        if (midpoint - 0.333_312_988_281_25).abs() >= 1.0e-15 {
            return Err("midpoint reference changed".to_string());
        }
        match midpoint_integral(1.0, 0, Ok::<f64, SpectralError>) {
            Err(error) => match input_error(error)? {
                SpectralError::PanelCount { .. } => Ok(()),
                _ => Err("wrong midpoint refusal".to_string()),
            },
            Ok(_) => Err("zero-panel midpoint rule was not refused".to_string()),
        }
    }

    #[test]
    fn monotone_exponential_enclosure_refuses_nonfinite_samples() -> Result<(), String> {
        let enclosure =
            monotone_exponential_enclosure(8.0, 32, |x| Ok::<f64, SpectralError>(1.0 - (-x).exp()))
                .map_err(|e| e.to_string())?;
        check(enclosure.monotone, "monotone factor rejected")?;
        match monotone_exponential_enclosure(8.0, 8, |_| Ok::<f64, SpectralError>(f64::NAN)) {
            Err(error) => match input_error(error)? {
                SpectralError::NonFinite { .. } => Ok(()),
                _ => Err("wrong non-finite enclosure refusal".to_string()),
            },
            Ok(_) => Err("NaN enclosure sample was accepted".to_string()),
        }
    }

    #[test]
    fn midpoint_trapezoid_enclosure_checks_segments() -> Result<(), String> {
        let (lower, upper, panels) =
            midpoint_trapezoid_enclosure(&[(0.0, 1.0, 4)], |x| x, |x| x, |x| x)
                .map_err(|e| e.to_string())?;
        check(
            lower.to_bits() == 0.5_f64.to_bits(),
            "midpoint enclosure lower",
        )?;
        check(
            upper.to_bits() == 0.5_f64.to_bits(),
            "trapezoid enclosure upper",
        )?;
        check(panels == 4, "enclosure panel count")?;
        match midpoint_trapezoid_enclosure(&[(0.0, 1.0, 0)], |x| x, |x| x, |x| x) {
            Err(SpectralError::PanelCount { .. }) => Ok(()),
            _ => Err("zero-count enclosure segment was not refused".to_string()),
        }
    }

    #[test]
    fn log_composite_trapezoid_checks_terms() -> Result<(), String> {
        let result = log_composite_trapezoid(&[(0.0, 1.0)], 8, |_| {
            Ok::<Option<f64>, SpectralError>(Some(0.0))
        })
        .map_err(|e| e.to_string())?
        .ok_or_else(|| "empty log composite result".to_string())?;
        check(result.abs() < 1.0e-15, "constant log composite result")?;
        match log_composite_trapezoid(&[(0.0, 1.0)], 4, |_| {
            Ok::<Option<f64>, SpectralError>(Some(f64::NAN))
        }) {
            Err(error) => match input_error(error)? {
                SpectralError::NonFinite { .. } => Ok(()),
                _ => Err("wrong log composite refusal".to_string()),
            },
            Ok(_) => Err("NaN log-composite term was accepted".to_string()),
        }
    }

    #[test]
    fn log_composite_by_exposes_weight_association() -> Result<(), String> {
        let result = log_composite_trapezoid_by(&[(0.0, 1.0)], 4, |_, spacing, weight| {
            Ok::<Option<f64>, SpectralError>(Some(spacing + weight))
        })
        .map_err(|e| e.to_string())?
        .ok_or_else(|| "empty associated log result".to_string())?;
        check(result.abs() < 1.0e-15, "associated log composite result")
    }

    #[test]
    fn log_stable_trapezoid_checks_nonfinite_values() -> Result<(), String> {
        let result = log_stable_trapezoid(0.0, 1.0, 8, |_| 0.0).map_err(|e| e.to_string())?;
        check(result.to_bits() == 1.0_f64.to_bits(), "log-stable constant")?;
        match log_stable_trapezoid(0.0, 1.0, 8, |_| f64::NAN) {
            Err(SpectralError::NonFinite { .. }) => Ok(()),
            _ => Err("NaN log-stable sample was not refused".to_string()),
        }
    }

    fn exponential_samples() -> [ExponentialSample; 2] {
        [
            ExponentialSample {
                coordinate: 0.0,
                exponent: 0.0,
                amplitude: 1.0,
            },
            ExponentialSample {
                coordinate: 1.0,
                exponent: 1.0,
                amplitude: 2.0,
            },
        ]
    }

    #[test]
    fn piecewise_exponential_integral_preserves_reference() -> Result<(), String> {
        let (trapezoid, linear) =
            piecewise_exponential_integral(&exponential_samples()).map_err(|e| e.to_string())?;
        check(
            trapezoid.to_bits() == (0.5 * (1.0 + 2.0 * (-1.0_f64).exp())).to_bits(),
            "piecewise trapezoid changed",
        )?;
        check(
            linear > 0.0 && linear.is_finite(),
            "piecewise linear integral",
        )
    }

    #[test]
    fn piecewise_exponential_by_refuses_nan_coordinate() -> Result<(), String> {
        let samples = [(0.0, 0.0, 1.0), (f64::NAN, 1.0, 2.0)];
        match piecewise_exponential_integral_by(&samples, |&(coordinate, exponent, amplitude)| {
            ExponentialSample {
                coordinate,
                exponent,
                amplitude,
            }
        }) {
            Err(SpectralError::NonFinite { .. }) => Ok(()),
            _ => Err("NaN piecewise coordinate was not refused".to_string()),
        }
    }

    #[test]
    fn piecewise_exponential_panel_checks_width() -> Result<(), String> {
        let value =
            piecewise_exponential_panel(1.0, 1.0, 2.0, 0.0, 1.0).map_err(|e| e.to_string())?;
        check(value > 0.0, "piecewise panel value")?;
        match piecewise_exponential_panel(0.0, 1.0, 2.0, 0.0, 1.0) {
            Err(SpectralError::Positive { .. }) => Ok(()),
            _ => Err("zero-width piecewise panel was not refused".to_string()),
        }
    }

    #[test]
    fn clenshaw_curtis_panel_integrates_and_refuses_nan() -> Result<(), String> {
        let panel = clenshaw_curtis_panel(0.0, 1.0, 8, |x| {
            Ok::<(f64, f64), SpectralError>((x * x, 0.0))
        })
        .map_err(|e| e.to_string())?;
        check(
            (panel.value - 1.0 / 3.0).abs() < 1.0e-14,
            "Clenshaw-Curtis value",
        )?;
        match clenshaw_curtis_panel(0.0, 1.0, 8, |_| {
            Ok::<(f64, f64), SpectralError>((f64::NAN, 0.0))
        }) {
            Err(error) => match input_error(error)? {
                SpectralError::NonFinite { .. } => Ok(()),
                _ => Err("wrong Clenshaw-Curtis refusal".to_string()),
            },
            Ok(_) => Err("NaN Clenshaw-Curtis sample was accepted".to_string()),
        }
    }

    #[test]
    fn clenshaw_and_bound_evaluate_a_linear_series() -> Result<(), String> {
        let coefficients = [1.0, 2.0];
        check(
            clenshaw(&coefficients, 0.25).to_bits() == 1.5_f64.to_bits(),
            "Clenshaw value",
        )?;
        let (value, error) = clenshaw_bound(&coefficients, 0.25);
        check(
            value.to_bits() == 1.5_f64.to_bits(),
            "bounded Clenshaw value",
        )?;
        check(error > 0.0 && error.is_finite(), "Clenshaw error bound")
    }

    #[test]
    fn chebyshev_antiderivative_differentiates_a_linear_series() -> Result<(), String> {
        let integrated = chebyshev_antiderivative(&[1.0, 2.0]);
        check(integrated.len() == 4, "antiderivative length")?;
        check(
            integrated[1].to_bits() == 1.0_f64.to_bits(),
            "constant antiderivative",
        )?;
        check(
            integrated[2].to_bits() == 0.5_f64.to_bits(),
            "linear antiderivative",
        )
    }

    #[test]
    fn direct_synthesis_public_items_preserve_endpoints() -> Result<(), String> {
        let cosine = direct_cosine_synthesize(&[1.0, 0.5], 2, 1, 4).map_err(|e| e.to_string())?;
        let sine = direct_sine_synthesize(&[0.0, 1.0], 2, 1, 4).map_err(|e| e.to_string())?;
        check(
            cosine[0].to_bits() == 1.5_f64.to_bits(),
            "direct cosine start",
        )?;
        check(sine[0].to_bits() == 0.0_f64.to_bits(), "direct sine start")
    }

    #[test]
    fn chirp_map_round_trip_and_typed_refusals() -> Result<(), String> {
        let map = ChirpMap::new(0.9, 0.08).map_err(|e| e.to_string())?;
        check(map.span().is_finite() && map.span() > 0.0, "chirp span")?;
        let point = map.inverse(0.7).map_err(|e| e.to_string())?;
        let phi = map.forward(point.theta).map_err(|e| e.to_string())?;
        check((phi - 0.7).abs() < 2.0e-15, "chirp round trip")?;
        match map.inverse(f64::NAN) {
            Err(SpectralError::MapCoordinate { .. }) => {}
            _ => return Err("NaN inverse coordinate was not refused".to_string()),
        }
        let folded_endpoint = PI + f64::EPSILON * PI;
        check(
            map.forward(folded_endpoint)
                .map_err(|error| error.to_string())?
                .to_bits()
                == PI.to_bits(),
            "rounded folded endpoint",
        )?;
        let outside = PI + 8.0 * f64::EPSILON * PI;
        match map.forward(outside) {
            Err(SpectralError::MapCoordinate { value, .. }) => check(
                value.to_bits() == outside.to_bits(),
                "forward refusal detail",
            )?,
            _ => return Err("out-of-range forward coordinate was not refused".to_string()),
        }
        match ChirpMap::new(0.9, 1.0e-320) {
            Err(SpectralError::ChirpWidth { value }) => check(
                value.to_bits() == 1.0e-320_f64.to_bits(),
                "chirp width detail",
            ),
            _ => Err("tiny chirp width was not refused".to_string()),
        }
    }

    #[test]
    fn chirp_grid_keeps_endpoints_and_refuses_huge_size() -> Result<(), String> {
        let map = ChirpMap::new(0.9, 0.08).map_err(|e| e.to_string())?;
        let grid = ChirpGrid::new(map, 7.0, 32).map_err(|e| e.to_string())?;
        check(grid.nodes() == 33, "chirp node count")?;
        check(
            grid.phi[0].to_bits() == 0.0_f64.to_bits(),
            "chirp phi start",
        )?;
        check(grid.phi[32].to_bits() == PI.to_bits(), "chirp phi end")?;
        check(
            grid.tau[0].to_bits() == 0.0_f64.to_bits(),
            "chirp tau start",
        )?;
        check(grid.tau[32].to_bits() == 3.5_f64.to_bits(), "chirp tau end")?;
        match ChirpGrid::new(map, 7.0, usize::MAX) {
            Err(SpectralError::SizeOverflow { .. }) => Ok(()),
            _ => Err("huge chirp grid was not refused".to_string()),
        }
    }

    #[test]
    fn evaluation_error_preserves_callback_failure() -> Result<(), String> {
        match midpoint_integral(1.0, 2, |_| Err::<f64, _>("callback")) {
            Err(EvaluationError::Evaluator("callback")) => Ok(()),
            _ => Err("callback error was not preserved".to_string()),
        }
    }

    #[test]
    fn nonincreasing_piecewise_coordinates_are_refused() -> Result<(), String> {
        let samples = [
            ExponentialSample {
                coordinate: 0.0,
                exponent: 0.0,
                amplitude: 1.0,
            },
            ExponentialSample {
                coordinate: 0.0,
                exponent: 1.0,
                amplitude: 2.0,
            },
        ];
        match piecewise_exponential_integral(&samples) {
            Err(SpectralError::Interval { .. }) => Ok(()),
            _ => Err("nonincreasing piecewise coordinates were not refused".to_string()),
        }
    }
}
