//! Utilities for ITM.
//!
//! These are miscellaneous functions that are functionally pure but also
//! implement well-known algorithms and formulae. They might be reusable
//! elsewhere and/or they might benefit from optimisation or being replaced by
//! calls to more efficient or correct versions.

/// Least-squares linear fit over evenly-spaced data between two points.
///
/// Returns _Z₀_ and _Zn_ for the line _Y = Z₀ + Zn × X_.
///
/// See ITM section `<53>`.
pub fn least_squares_linear_fit(interval: f64, data: &[f64], points: (f64, f64)) -> (f64, f64) {
    let xn = data.len() as f64;
    let mut xa = fortran_dim(points.0 / interval, 0.0) as isize as f64;
    let mut xb = xn - (fortran_dim(xn, points.1 / interval) as isize as f64);

    if xb <= xa {
        xa = fortran_dim(xa, 1.0);
        xb = xn - fortran_dim(xn, xb + 1.0);
    }

    let mut ja = xa as usize;
    let jb = xb as usize;

    let n = jb - ja;
    xa = xb - xa;

    let mut xx = -0.5 * xa;
    xb += xx;

    let mut a = 0.5 * (data[ja] + data[jb]);
    let mut b = 0.5 * (data[ja] - data[jb]) * xx;

    for _ in 2..=n {
        ja += 1;
        xx += 1.0;
        a += data[ja];
        b += data[ja] * xx;
    }

    a /= xa;
    b = b * 12.0 / ((xa * xa + 2.0) * xa);

    (a - b * xb, a + b * (xn - xb))
}

/// Fortran-style DIM operation.
pub fn fortran_dim(x: f64, y: f64) -> f64 {
    (x - y).max(0.0)
}
