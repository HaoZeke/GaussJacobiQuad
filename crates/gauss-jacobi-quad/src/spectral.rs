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

/// Invalid input to a mapped spectral grid.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum SpectralError {
    #[error("the chirp centre must be finite and inside (0, pi)")]
    ChirpCentre,
    #[error("the chirp width must be finite and positive")]
    ChirpWidth,
    #[error("a mapped grid needs at least two panels")]
    PanelCount,
    #[error("the mapped-grid period must be finite and positive")]
    Period,
}

/// Endpoint-halved uniform cosine grid on a half period.
///
/// The returned vectors are `(theta, tau, weights)`, with
/// `theta[j] = pi*j/m`, `tau[j] = beta*theta[j]/(2*pi)`, and trapezoid
/// weights `beta/m`, halved at both endpoints. This operation order is kept
/// stable because downstream spectral calculations depend on bitwise
/// reproducibility.
pub fn half_period_grid(m: usize, beta: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut theta = vec![0.0; m + 1];
    let mut tau = vec![0.0; m + 1];
    let mut weights = vec![beta / m as f64; m + 1];
    for j in 0..=m {
        theta[j] = PI * j as f64 / m as f64;
        tau[j] = beta * theta[j] / (2.0 * PI);
    }
    weights[0] *= 0.5;
    weights[m] *= 0.5;
    (theta, tau, weights)
}

/// Two-point Gauss--Legendre rule on an interval described by centre and width.
pub fn gauss_legendre_two_point(centre: f64, width: f64) -> [(f64, f64); 2] {
    let offset = width / (2.0 * 3.0_f64.sqrt());
    let weight = 0.5 * width;
    [(centre - offset, weight), (centre + offset, weight)]
}

/// Endpoint-halved trapezoid for a fixed-size vector-valued integrand.
pub fn endpoint_trapezoid_array<const N: usize, F>(
    lower: f64,
    upper: f64,
    panels: usize,
    mut value: F,
) -> [f64; N]
where
    F: FnMut(f64) -> [f64; N],
{
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
            sum[component] += weight * sampled[component];
        }
    }
    sum
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

/// Unnormalized forward FFT with a plan cached on the calling thread.
pub fn forward_inplace(buf: &mut [Complex<f64>]) {
    let len = buf.len();
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        plan.process(buf);
    });
}

/// Weighted cosine projection on an endpoint-halved grid of `m` panels.
pub fn cc_project(samples: &[f64], m: usize, k_count: usize, d: usize, beta: f64) -> Vec<f64> {
    if m < 2 || k_count > m + 1 {
        return direct_project(samples, m, k_count, d, beta);
    }
    let mut out = vec![0.0; k_count * d];
    let scale = beta / (2.0 * m as f64);
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(2 * m);
        let mut wave = vec![Complex::new(0.0, 0.0); 2 * m];
        for a in 0..d {
            wave[0] = Complex::new(samples[a], 0.0);
            wave[m] = Complex::new(samples[m * d + a], 0.0);
            for t in 1..m {
                let value = Complex::new(samples[t * d + a], 0.0);
                wave[t] = value;
                wave[2 * m - t] = value;
            }
            plan.process(&mut wave);
            for k in 0..k_count {
                out[k * d + a] = wave[k].re * scale;
            }
        }
    });
    out
}

/// Cosine-series synthesis on an endpoint-halved grid.
pub fn cc_synthesize(coeff: &[f64], k_count: usize, d: usize, m: usize) -> Vec<f64> {
    if m < 2 || k_count > m + 1 {
        return direct_cosine_synthesize(coeff, k_count, d, m);
    }
    let mut out = vec![0.0; (m + 1) * d];
    let len = 2 * m;
    INVERSE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        let mut spectrum = vec![Complex::new(0.0, 0.0); len];
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
    });
    out
}

/// Weighted sine projection on an endpoint-halved grid.
pub fn dst_project(samples: &[f64], m: usize, k_count: usize, d: usize, beta: f64) -> Vec<f64> {
    let mut out = vec![0.0; k_count * d];
    let scale = -beta / (2.0 * m as f64);
    FORWARD.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(2 * m);
        let mut wave = vec![Complex::new(0.0, 0.0); 2 * m];
        for a in 0..d {
            wave.fill(Complex::new(0.0, 0.0));
            for t in 1..m {
                let value = samples[t * d + a];
                wave[t] = Complex::new(value, 0.0);
                wave[2 * m - t] = Complex::new(-value, 0.0);
            }
            plan.process(&mut wave);
            for k in 1..k_count {
                out[k * d + a] = wave[k].im * scale;
            }
        }
    });
    out
}

/// Sine-series synthesis on an endpoint-halved grid.
pub fn dst_synthesize(coeff: &[f64], k_count: usize, d: usize, m: usize) -> Vec<f64> {
    if m < 2 {
        return direct_sine_synthesize(coeff, k_count, d, m);
    }
    let mut out = vec![0.0; (m + 1) * d];
    let len = 2 * m;
    INVERSE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let plan = cache.plan(len);
        let mut spectrum = vec![Complex::new(0.0, 0.0); len];
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
    });
    out
}

/// Direct weighted cosine projection, provided as a small-grid reference.
pub fn direct_project(samples: &[f64], m: usize, k_count: usize, d: usize, beta: f64) -> Vec<f64> {
    let (theta, _, weights) = half_period_grid(m, beta);
    let mut out = vec![0.0; k_count * d];
    for k in 0..k_count {
        for t in 0..=m {
            let weighted_cosine = weights[t] * (theta[t] * k as f64).cos();
            for a in 0..d {
                out[k * d + a] += weighted_cosine * samples[t * d + a];
            }
        }
    }
    out
}

/// Kahan-summed midpoint rule on `[0, upper]`.
///
/// The callback receives each panel midpoint. Its error type is preserved.
pub fn midpoint_integral<E, F>(upper: f64, panels: usize, mut value: F) -> Result<f64, E>
where
    F: FnMut(f64) -> Result<f64, E>,
{
    let step = upper / panels as f64;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let point = (panel as f64 + 0.5) * step;
        let term = value(point)?;
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    Ok(step * sum)
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
) -> Result<MonotoneExponentialEnclosure, E>
where
    F: FnMut(f64) -> Result<f64, E>,
{
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
        let current = value(right)?;
        if current < previous {
            monotone = false;
            break;
        }
        let weight = (-left).exp() * panel_weight_factor;
        kahan_add(&mut lower, &mut lower_correction, previous * weight);
        kahan_add(&mut upper_sum, &mut upper_correction, current * weight);
        previous = current;
    }
    Ok(MonotoneExponentialEnclosure {
        lower,
        upper: upper_sum + (-upper).exp(),
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
) -> (T, T, usize)
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
    L: FnMut(T) -> T,
    U: FnMut(T) -> T,
    X: Fn(f64) -> T,
{
    let mut lower = exact(0.0);
    let mut upper = exact(0.0);
    let mut total_panels = 0;
    for &(start, end, count) in segments {
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
            total_panels += 1;
        }
    }
    (lower, upper, total_panels)
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
) -> Result<Option<f64>, E>
where
    F: FnMut(f64) -> Result<Option<f64>, E>,
{
    let mut terms = Vec::with_capacity(segments.len() * (panels + 1));
    for &(lower, upper) in segments {
        let spacing = (upper - lower) / panels as f64;
        for index in 0..=panels {
            let point = lower + spacing * index as f64;
            if let Some(log_integrand) = log_value(point)? {
                let weight: f64 = if index == 0 || index == panels {
                    0.5
                } else {
                    1.0
                };
                terms.push(spacing.ln() + weight.ln() + log_integrand);
            }
        }
    }
    let Some(largest) = terms.iter().copied().max_by(f64::total_cmp) else {
        return Ok(None);
    };
    let scaled_sum: f64 = terms.iter().map(|term| (*term - largest).exp()).sum();
    Ok(Some(largest + scaled_sum.ln()))
}

/// Composite trapezoid for a positive log integrand.
///
/// Each panel is scaled by its larger endpoint before exponentiation.
pub fn log_stable_trapezoid<F>(lower: f64, upper: f64, panels: usize, mut log_value: F) -> f64
where
    F: FnMut(f64) -> f64,
{
    let mut sum = 0.0;
    let mut previous = log_value(lower);
    for index in 0..panels {
        let left = lower + (upper - lower) * index as f64 / panels as f64;
        let right = lower + (upper - lower) * (index + 1) as f64 / panels as f64;
        let next = log_value(right);
        let width = right - left;
        let top = if previous >= next { previous } else { next };
        sum += width * top.exp() * 0.5 * ((previous - top).exp() + (next - top).exp());
        previous = next;
    }
    sum
}

/// One sample for piecewise-linear exponential quadrature.
#[derive(Clone, Copy, Debug)]
pub struct ExponentialSample {
    pub coordinate: f64,
    pub exponent: f64,
    pub amplitude: f64,
}

/// Trapezoid and exact integral of piecewise-linear amplitude and exponent.
pub fn piecewise_exponential_integral(samples: &[ExponentialSample]) -> (f64, f64) {
    piecewise_exponential_integral_by(samples, |sample| *sample)
}

/// Generic form of [`piecewise_exponential_integral`] without a copy.
pub fn piecewise_exponential_integral_by<T, F>(samples: &[T], sample: F) -> (f64, f64)
where
    F: Fn(&T) -> ExponentialSample,
{
    let mut trapezoid = 0.0;
    let mut linear = 0.0;
    for index in 0..samples.len().saturating_sub(1) {
        let left = sample(&samples[index]);
        let right = sample(&samples[index + 1]);
        let width = right.coordinate - left.coordinate;
        if width <= 0.0 {
            continue;
        }
        let amplitude0 = left.amplitude;
        let amplitude1 = right.amplitude;
        let exponent0 = left.exponent;
        let exponent1 = right.exponent;
        let value0 = (-exponent0).exp();
        let value1 = (-exponent1).exp();
        trapezoid += 0.5 * (amplitude0 * value0 + amplitude1 * value1) * width;
        linear += piecewise_exponential_panel(width, amplitude0, amplitude1, exponent0, exponent1);
    }
    (trapezoid, linear)
}

/// Exact panel integral for linear amplitude and exponent interpolants.
pub fn piecewise_exponential_panel(
    width: f64,
    amplitude0: f64,
    amplitude1: f64,
    exponent0: f64,
    exponent1: f64,
) -> f64 {
    let exponent_delta = exponent1 - exponent0;
    let amplitude_delta = amplitude1 - amplitude0;
    let value0 = (-exponent0).exp();
    if exponent_delta.abs() < 1.0e-10 {
        width * value0 * (amplitude0 + 0.5 * amplitude_delta)
    } else {
        let decay = (-exponent_delta).exp();
        let constant_integral = (1.0 - decay) / exponent_delta;
        let linear_integral =
            (1.0 - decay * (1.0 + exponent_delta)) / (exponent_delta * exponent_delta);
        width * value0 * (amplitude0 * constant_integral + amplitude_delta * linear_integral)
    }
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
) -> Result<ClenshawCurtisPanel, E>
where
    F: FnMut(f64) -> Result<(f64, f64), E>,
{
    let middle = 0.5 * (lower + upper);
    let half = 0.5 * (upper - lower);
    let mut samples = vec![0.0; degree + 1];
    let mut errors = vec![0.0; degree + 1];
    for index in 0..=degree {
        let coordinate = (PI * index as f64 / degree as f64).cos();
        let (value, error) = sample(middle + half * coordinate)?;
        samples[index] = value;
        errors[index] = error;
    }
    let mut coefficients = vec![0.0; degree + 1];
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
    let mut previous = 0.0_f64;
    let mut monotone = true;
    for index in 0..=degree {
        let coordinate = (PI * index as f64 / degree as f64).cos();
        let above = half * (at_plus - clenshaw(&antiderivative, coordinate));
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
pub fn direct_cosine_synthesize(coeff: &[f64], k_count: usize, d: usize, m: usize) -> Vec<f64> {
    direct_synthesize(coeff, k_count, d, m, false)
}

/// Direct sine-series synthesis, provided as a small-grid reference.
pub fn direct_sine_synthesize(coeff: &[f64], k_count: usize, d: usize, m: usize) -> Vec<f64> {
    direct_synthesize(coeff, k_count, d, m, true)
}

fn direct_synthesize(coeff: &[f64], k_count: usize, d: usize, m: usize, sine: bool) -> Vec<f64> {
    let (theta, _, _) = half_period_grid(m, 1.0);
    let mut out = vec![0.0; (m + 1) * d];
    for t in 0..=m {
        for k in 0..k_count {
            let phase = theta[t] * k as f64;
            let basis = if sine { phase.sin() } else { phase.cos() };
            for a in 0..d {
                out[t * d + a] += basis * coeff[k * d + a];
            }
        }
    }
    out
}

fn sinhc_pair(h: f64) -> (f64, f64) {
    let h2 = h * h;
    let mut value = 1.0;
    let mut slope = 0.0;
    let mut e = 1.0 / 6.0;
    let mut n = 1.0_f64;
    loop {
        let term = e * h2;
        let slope_term = 2.0 * n * e;
        let settled = term <= f64::EPSILON * value && slope_term <= f64::EPSILON * slope;
        value += term;
        slope += slope_term;
        if settled {
            break;
        }
        n += 1.0;
        e *= h2 / ((2.0 * n) * (2.0 * n + 1.0));
    }
    (value, slope)
}

fn asinh_gap(hi: f64, lo: f64, gap: f64) -> f64 {
    let roots = (1.0 + hi * hi).sqrt() + (1.0 + lo * lo).sqrt();
    let bracket = 1.0 + (hi + lo) / roots;
    (gap * bracket * (-lo.asinh()).exp()).ln_1p()
}

/// Monotone arcsinh coordinate map from a uniform phase to physical time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpMap {
    pub theta_star: f64,
    pub sigma: f64,
    s_plus: f64,
    s_minus: f64,
}

/// One mapped node and its first two coordinate derivatives.
#[derive(Debug, Clone, Copy)]
pub struct MapPoint {
    pub theta: f64,
    pub g: f64,
    pub g_phi: f64,
}

impl ChirpMap {
    pub fn new(theta_star: f64, sigma: f64) -> Result<Self, SpectralError> {
        if !(theta_star > 0.0 && theta_star < PI && theta_star.is_finite()) {
            return Err(SpectralError::ChirpCentre);
        }
        if !(sigma > 0.0 && sigma.is_finite()) {
            return Err(SpectralError::ChirpWidth);
        }
        let half = 0.5 * theta_star;
        let above = 2.0 * half.sin().powi(2) / sigma;
        let below = 2.0 * half.cos().powi(2) / sigma;
        Ok(Self {
            theta_star,
            sigma,
            s_plus: above.asinh(),
            s_minus: -below.asinh(),
        })
    }

    pub fn span(&self) -> f64 {
        self.s_plus - self.s_minus
    }

    pub fn inverse(&self, phi: f64) -> MapPoint {
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
        let (c_plus, e_plus) = sinhc_pair(h_plus);
        let (c_minus, e_minus) = sinhc_pair(h_minus);
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
        MapPoint {
            theta,
            g,
            g_phi: g * dlog,
        }
    }

    pub fn forward(&self, theta: f64) -> f64 {
        let above = 2.0 * (0.5 * theta).sin().powi(2) / self.sigma;
        let below = 2.0 * (0.5 * theta).cos().powi(2) / self.sigma;
        let a = self.s_plus.sinh();
        let c = -self.s_minus.sinh();
        let b = a - above;
        let upper = asinh_gap(a, b, above);
        let lower = asinh_gap(c, -b, below);
        2.0 * upper.sqrt().atan2(lower.sqrt())
    }
}

/// Uniform phase nodes with physical times and map rates.
#[derive(Debug, Clone)]
pub struct ChirpGrid {
    pub m: usize,
    pub beta: f64,
    pub map: ChirpMap,
    pub phi: Vec<f64>,
    pub tau: Vec<f64>,
    pub rate: Vec<f64>,
    pub weights: Vec<f64>,
}

impl ChirpGrid {
    pub fn new(map: ChirpMap, beta: f64, m: usize) -> Result<Self, SpectralError> {
        if m < 2 {
            return Err(SpectralError::PanelCount);
        }
        if !(beta > 0.0 && beta.is_finite()) {
            return Err(SpectralError::Period);
        }
        let omega = 2.0 * PI / beta;
        let mut phi = Vec::with_capacity(m + 1);
        let mut tau = Vec::with_capacity(m + 1);
        let mut rate = Vec::with_capacity(m + 1);
        let mut weights = vec![PI / m as f64; m + 1];
        weights[0] *= 0.5;
        weights[m] *= 0.5;
        for t in 0..=m {
            let phase = PI * t as f64 / m as f64;
            let point = map.inverse(phase);
            let mapped_rate = omega / point.g;
            phi.push(phase);
            tau.push(point.theta / omega);
            rate.push(mapped_rate);
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

    pub fn nodes(&self) -> usize {
        self.m + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fft_projection_matches_direct_rule() -> Result<(), &'static str> {
        let m = 32;
        let modes = 16;
        let dimensions = 3;
        let beta = 2.5;
        let samples: Vec<f64> = (0..((m + 1) * dimensions))
            .map(|index| (index as f64 * 0.17).sin())
            .collect();
        let fft = cc_project(&samples, m, modes, dimensions, beta);
        let direct = direct_project(&samples, m, modes, dimensions, beta);
        for (left, right) in fft.iter().zip(&direct) {
            if (left - right).abs() >= 1.0e-12 {
                return Err("FFT and direct cosine projections differ");
            }
        }
        Ok(())
    }

    #[test]
    fn generic_integrators_preserve_reference_values() -> Result<(), &'static str> {
        let midpoint = midpoint_integral(1.0, 64, |point| Ok::<f64, &'static str>(point * point))?;
        if (midpoint - 0.333_312_988_281_25).abs() >= 1.0e-15 {
            return Err("midpoint reference changed");
        }
        let samples = [
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
        ];
        let (trapezoid, linear) = piecewise_exponential_integral(&samples);
        if trapezoid.to_bits() != (0.5 * (1.0 + 2.0 * (-1.0_f64).exp())).to_bits() {
            return Err("piecewise trapezoid changed");
        }
        if !(linear > 0.0 && linear.is_finite()) {
            return Err("piecewise-linear integral is invalid");
        }
        Ok(())
    }

    #[test]
    fn chirp_round_trip_keeps_endpoints() -> Result<(), SpectralError> {
        let map = ChirpMap::new(0.9, 0.08)?;
        let grid = ChirpGrid::new(map, 7.0, 32)?;
        if grid.phi[0].to_bits() != 0.0_f64.to_bits()
            || grid.phi[32].to_bits() != PI.to_bits()
            || grid.tau[0].to_bits() != 0.0_f64.to_bits()
            || grid.tau[32].to_bits() != (3.5_f64).to_bits()
        {
            return Err(SpectralError::Period);
        }
        Ok(())
    }
}
