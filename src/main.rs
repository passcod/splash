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

/// Computes a distance on the Earth between two points.
///
/// Uses “inverse” geodesic equations from [Karney 2013] rather than Haversine
/// (which has large errors) or Valenty’s (which is unstable for some inputs).
/// Assumes the WGS84 geodesic.
///
/// While inverse geodesic equations can have multiple solutions, the distance
/// remains the same for all results. Thus, this implementation ignores those.
///
/// This function was written with great help from [geographic]’s Geodesy Java
/// library, written by the author of the paper cited above.
///
/// [Karney 2013]: https://doi.org/10.1007/s00190-012-0578-z
/// [geographic]: https://geographiclib.sourceforge.io/html/java/
fn geodesic_distance(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    //
}

/// TX or RX site definition
struct Site {
    lat: f64,
    lon: f64,
    altitude: f64,
    name: String,
    filename: PathBuf,
}

struct RFPath {
    lats: Vec<f64>,
    lons: Vec<f64>,
    elevations: Vec<f64>,
    distances: Vec<f64>,
    length: f64,
}

/// Preprocessed elevation data block of regular gridding.
struct Dem {
    latitude_min: f64,
    latitude_max: f64,
    latitude_num: usize,
    longitude_min: f64,
    longitude_max: f64,
    longitude_num: usize,
    elevation_min: f64,
    elevation_max: f64,
    data: Vec<Vec<f64>>,
}

/// Parameters for a Site
struct Lrp {
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
