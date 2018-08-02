#![forbid(unsafe_code)]
#![cfg_attr(feature = "cargo-clippy", deny(clippy_pedantic))]

extern crate num_complex;
extern crate gdal;
extern crate geo;

use geo::Point;
use itm::*;

mod itm;

const EQUATORIAL_RADIUS: f64 = 6_378_137.0;
const POLAR_RADIUS: f64 = 6_356_752.3;

/// Computes the geocentric radius (in metres) at a given latitude (in degrees).
pub fn local_radius(latitude: f64) -> f64 {
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

/// Base site definition.
#[derive(Clone, Debug)]
pub struct Site {
    /// Where it is
    position: Point<f64>,

    /// How high above ground
    aboveground: f64,

    /// What it's called
    name: String,
}

/// An RF transmitter and its parameters.
#[derive(Clone, Debug)]
pub struct Transmitter {
    site: Site,

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

fn main() {
    println!("Hello, world!");
}
