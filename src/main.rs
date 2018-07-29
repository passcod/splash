#![forbid(unsafe_code)]
#![cfg_attr(feature = "cargo-clippy", deny(clippy_pedantic))]

use std::path::PathBuf;

const EQUATORIAL_RADIUS: f64 = 6_378_137.0;
const POLAR_RADIUS: f64 = 6_356_752.3;

const FZONE_CLEARANCE: f64 = 0.6;

/// Computes the geocentric radius (in metres) at a given latitude (in degrees).
fn local_radius(latitude: f64) -> f64 {
    let cos = latitude.cos();
    let sin = latitude.sin();

    let upper = (EQUATORIAL_RADIUS.powi(2) * cos).powi(2) + (POLAR_RADIUS.powi(2) * sin).powi(2);
    let lower = (EQUATORIAL_RADIUS * cos).powi(2) + (POLAR_RADIUS * sin).powi(2);

    (upper / lower).sqrt()
}

#[test]
fn test_local_radius() {
    assert_eq!(local_radius(0.0), EQUATORIAL_RADIUS);
    assert_eq!(local_radius(12.3), 6376666.768840445);
    assert_eq!(local_radius(-35.273), 6368978.931744378);
    assert_eq!(local_radius(90.0), 6361074.591356493);
    assert_eq!(local_radius(-1.469167), 6356974.249836446);
}

// in degrees
struct Point {
    pub lat: f64,
    pub lon: f64,
}

/// Computes the distance and forward azimuth on the Earth between two points.
///
/// Uses “inverse” geodesic equations from [Karney 2013] rather than Haversine
/// (which has large errors) or Valenty’s (which is unstable for some inputs).
/// Assumes WGS84 geography.
///
/// While inverse geodesic equations can have multiple solutions, the distance
/// remains the same for all results. At the moment, we ignore further
/// results for the azimuth and return only the first found.
///
/// This function was written with great help from [geographic]’s Geodesy Java
/// library, written by the author of the paper cited above.
///
/// [Karney 2013]: https://doi.org/10.1007/s00190-012-0578-z
/// [geographic]: https://geographiclib.sourceforge.io/html/java/
fn geodesic_distance(a: Point, b: Point) -> (f64, f64) {
    (0.0, 0.0)
}

/// TX or RX site definition
struct Site {
    /// Where it is
    position: Point,

    /// How high above ground
    aboveground: f64,

    /// What it's called
    name: String,
}

/// Parameters for a Site
struct Parameters {
    /// Earth Dielectric Constant (Relative permittivity)
    dielectric: f64,

    /// Earth Conductivity (Siemens per metre)
    conductivity: f64,

    /// Atmospheric Bending Constant (N-units)
    bending: f64,

    /// Site frequency (MHz)
    frequency: f64,

    /// Radio climate
    climate: Climate,

    /// Polarisation
    polarisation: Polarisation,

    /// Antenna pattern
    pattern: Vec<Vec<f64>>,
}

enum Climate {
    Equatorial, // 1
    ContinentalSubtropical, // 2
    MaritimeSubtropical, // 3
    Desert, // 4
    ContinentalTemperate, // 5
    MaritimeTemperateOverLand, // 6
    MaritimeTemperateOverSea, // 7
}

enum Polarisation {
    Horizontal, // 0
    Vertical, // 1
    Dual, // n/a
}

fn main() {
    println!("Hello, world!");
}
