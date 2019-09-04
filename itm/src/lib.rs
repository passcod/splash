//! The Longley-Rice Irregular Terrain Model for RF.
//!
//! This module was initially derived from a C++ translation of the original
//! FORTRAN implementation of the ITM, the [Irregular Terrain Model][ITM68],
//! also known as Longley-Rice, an empirical RF propagation model developed at
//! the U.S. National Telecommunications and Information Administration around
//! the 1960s by Anita Longley and Phil Rice.
//!
//! While the initial reconstruction was done referring to the above-mentioned
//! source and its adaptation in [SPLAT!], all variables and function names were
//! still unscrutable, as was the general operation of the model. Thus, the code
//! was then cross-referenced back to the [LaTeX documentation][ITM122] of the
//! ITM version 1.2.2 (sections within that document are referenced with `<N>`
//! in the code), and made more legible and understandable by deriving
//! meaningful names and incorporating documentation into this source by
//! referring back to George Hufford's 1999 memo describing
//! “[The Algorithm][GH1999]” (referenced as T.A. in the code).
//!
//! Note that the structure of the module is quite different from the above
//! implementations, making it more idiomatic to Rust and allowing more uses.
//! Remarkably, we are much more friendly to concurrent computes.
//!
//! This implementation is released in the Public Domain, although note that the
//! NTIA requests any use of the ITM is properly credited.
//!
//! [GH1999]: https://www.its.bldrdoc.gov/media/50676/itm_alg.pdf
//! [ITM122]: https://www.its.bldrdoc.gov/media/50674/itm.pdf
//! [ITM68]: https://www.its.bldrdoc.gov/resources/radio-propagation-software/itm/itm.aspx
//! [SPLAT!]: http://www.qsl.net/kd2bd/splat.html

#![forbid(unsafe_code)]
#![cfg_attr(feature = "cargo-clippy", deny(clippy_pedantic))]

use climate::{Climate, ClimateConstants};
use num_complex::Complex64;
use formulae::*;

pub mod climate;
pub mod formulae;

/// Propagation model instance.
///
/// Holds all state related to one instance of the Irregular Terrain Model, for
/// one particular transmitter, elevation profile from the transmitter outwards,
/// and set of options.
///
/// A `Model` instance can be used to query propagation at any distance from the
/// transmitter along the loaded elevation profile. All preliminary work is done
/// on initialisation: propagation queries only need immutable access and can
/// therefore be done concurrently.
#[derive(Clone, Debug, PartialEq)]
pub struct Model<'a> {
    length: f64, // <prop.dist>
    elevations: &'a Vec<f64>,
    heights: (f64, f64), // <hg>
    settings: &'a Settings,
    computed: Computed,
    cached: Cached,
}

/// Computed parameters of the model.
///
/// These are computed from the input settings and elevation profile.
///
/// In the source, `<letters>` indicate original variable names, such that one
/// may cross-reference these back to the memos and other implementations.
///
/// An exception is `<distance>`, here named `interval`. In the original, that
/// parameter was an input, and the total _length_ of the profile was computed.
/// That, along with the name, was confusing and possibly inflexible. Thus this
/// bit of logic is inverted.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Computed {
    /// [Wave number] of the carrier/central frequency (in radians per unit distance).
    ///
    /// [Wave number]: https://en.wikipedia.org/wiki/Wavenumber
    pub wave_number: f64, // <wn>

    /// General elevation: the mean elevation of the modeled system.
    pub general_elevation: f64, // <zsys>

    /// Earth's effective curvature at the system's elevation.
    pub effective_curvature: f64, // <gme>

    /// Effective surface refractivity at the system's elevation.
    pub effective_refractivity: f64, // <ens>

    /// Surface transfer impedance to the ground.
    pub transfer_impedance: Complex64, // <zgnd>

    /// Interval distance between each elevation.
    pub interval: f64, // <distance>

    /// Elevation angles of the horizons from each terminal at the heights of
    /// the antennas, in radians.
    pub elevation_angles: (f64, f64), // <the>

    /// Distances from each terminal to its radio horizon.
    pub horizon_distances: (f64, f64), // <dl>

    /// Terminal effective heights: adjusted against horizons (or obstructions).
    pub effective_heights: (f64, f64), // <he>

    pub terrain_irregularity: f64, // <dh>
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cached {
    pub lvar: isize,
    pub mdvar: isize,

    pub dmin: f64,
    pub xae: f64,
    pub wlos: bool,
    pub wscat: bool,

    pub ad: f64,
    pub rr: f64,
    pub etq: f64,
    pub h0s: f64,
}

impl Default for Cached {
    fn default() -> Self {
        Self {
            lvar: 5,
            mdvar: 12,

            dmin: 0.0,
            xae: 0.0,
            wlos: false,
            wscat: false,

            ad: 0.0,
            rr: 0.0,
            etq: 0.0,
            h0s: 0.0,
        }
    }
}

/// The polarisation of the radio wave.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Polarisation {
    Horizontal,
    Vertical,
}

impl Default for Polarisation {
    fn default() -> Self {
        Polarisation::Horizontal
    }
}

/// The carrier or central frequency (in MHz).
///
/// This is pre-computed (in the paper) assuming a speed of light **in air** of
/// 299.7 metres per microsecond. See Fig/1.1 in T.A.
const CARRIER_FREQUENCY: f64 = 47.7;

impl<'a> Model<'a> {
    /// Creates a new instance for a profile.
    ///
    ///  - `length` and `elevations` are the terrain input. We assume elevations
    ///    are along a line of length `length` metres, spaced equidistantly, and
    ///    that the transmitter is at distance zero.
    ///
    ///  - `heights` is a tuple of the _height above ground_ in metres of the
    ///    terminals (antennas), transmitter then receiver.
    ///
    ///  - `settings` is an instance of the `Settings` struct, containing
    ///    atmospheric, surface, radio, and statistical parameter values.
    ///
    /// Upon construction, we prepare all values needed for the individual
    /// calculations, perform bounds checking, and decide on the mode or modes
    /// of the model.
    ///
    /// The ITM has three different modes depending on the characteristics and
    /// distance modeled. We assume that all points within the given profile
    /// will be queried, so if the inputs span two or more modes, we prepare for
    /// each of those.
    ///
    ///  - If the _distance to site_ is less than line of sight, then path loss
    ///    is found from the Line-Of-Sight submodel, from optic horizon and the
    ///    two-ray model.
    ///
    ///  - If the _distance to site_ is larger than that, but still within the
    ///    radio horizon, then the Diffraction submodel is used, considering the
    ///    curvature of the Earth and knife-edge mechanisms.
    ///
    ///  - If the _distance to site_ is beyond the radio horizon, the Scatter
    ///    submodel is used, which computes constants for a linear relationship.
    ///
    /// After the attenuation is found, atmospheric and climate conditions are
    /// brought back to the global measurements made by the ITM team around the
    /// world and they are used to perform statistical quantile variability
    /// modelling given two input settings: % of time, and % of situations.
    pub fn new(
        length: f64,
        elevations: &'a Vec<f64>,
        heights: (f64, f64),
        settings: &'a Settings,
    ) -> Result<Self, String> {
        if elevations.len() < 3 {
            return Err("elevations should have more than 2 points".into());
        }

        if heights.0 < 0.0 || heights.1 < 0.0 {
            return Err("heights may not be negative".into());
        }

        let cached = Cached::default();
        let computed = Computed::default();
        let mut model = Self {
            length,
            elevations,
            heights,
            settings,
            computed,
            cached,
        };

        model.prepare_environment();
        model.find_horizons();
        // mdvar, lvar ???

        let xl = model.get_xl();
        model.computed.terrain_irregularity = d1thx(model.computed.interval, model.elevations, xl);

        model.adjust_horizons(xl);

        // lrprop(0.0, prop, propa);

        // mdvar = mode of variability

        Ok(model)
    }

    /// Various basic computed parameters.
    ///
    /// See ITM section `<41>`.
    fn prepare_environment(&mut self) {
        // Fig/1.1 in T.A.
        self.computed.wave_number = self.settings.frequency / CARRIER_FREQUENCY;

        let sum: f64 = self.elevations.iter().sum();
        self.computed.general_elevation = sum / (self.elevations.len() as f64);

        self.computed.effective_refractivity = self.settings.surface_refractivity;
        if self.computed.general_elevation != 0.0 {
            // Fig/1.2 in T.A.
            self.computed.effective_refractivity *=
                (-self.computed.general_elevation / 9460.0).exp();
        }

        // Fig/1.3 in T.A.
        self.computed.effective_curvature =
            157e-9 * (1.0 - 0.04665 * (self.computed.effective_refractivity / 179.3).exp());

        // Fig/1.5 in T.A.
        let complex_relative_permittivity = Complex64::new(
            self.settings.permittivity,
            376.62 * self.settings.conductivity / self.computed.wave_number,
        );

        // Fig/1.4 in T.A.
        self.computed.transfer_impedance = (complex_relative_permittivity - 1.0).sqrt();
        if self.settings.polarisation == Polarisation::Vertical {
            self.computed.transfer_impedance /= complex_relative_permittivity;
        }
    }

    /// Use the elevation profile to find the two horizons.
    ///
    /// See ITM section `<47>`.
    fn find_horizons(&mut self) {
        let len = self.elevations.len();
        let interval = self.length / (len - 1) as f64;
        self.computed.interval = interval; // see note in Computed

        // absolute heights of terminals
        let tx_z = self.elevations[0] + self.heights.0;
        let rx_z = self.elevations[len - 1] + self.heights.1;

        let half_curve = self.computed.effective_curvature / 2.0;
        let half_length = half_curve * self.length;

        let vertical_delta = (rx_z - tx_z) / self.length;
        let mut angle_tx = vertical_delta - half_length;
        let mut angle_rx = -vertical_delta - half_length;

        let mut dist_tx = self.length;
        let mut dist_rx = self.length;

        let mut along_tx = 0.0;
        let mut along_rx = self.length;

        let mut wq = true;

        // We advance along the elevation profile looking both from the TX and
        // from the RX at the same time. As we go we adjust the elevations to
        // the Earth's curvature.

        for i in 1..len {
            along_tx += interval;
            along_rx -= interval;

            let tx_adj = half_curve * along_tx + angle_tx;
            let tx_delta = self.elevations[i] - tx_adj * along_tx - tx_z;

            if tx_delta > 0.0 {
                angle_tx += half_length / along_tx;
                dist_tx = along_tx;
                wq = false;
            }

            if wq {
                continue;
            }

            let rx_adj = half_curve * along_rx + angle_rx;
            let rx_delta = self.elevations[i] - rx_adj * along_rx - rx_z;

            if rx_delta > 0.0 {
                angle_rx += half_length / along_rx;
                dist_rx = along_rx;
            }
        }

        self.computed.elevation_angles = (angle_tx, angle_rx);
        self.computed.horizon_distances = (dist_tx, dist_rx);
    }

    /// Intermediate values used in some prep calculations.
    ///
    /// The least of fifteen times the terminal height above ground, and 10% of
    /// the horizon distance.
    ///
    /// Unsure as to what it's about, really, but it's there.
    ///
    /// See ITM sections `<43>` and `<44>`.
    fn get_xl(&mut self) -> (f64, f64) {
        fn make_xl(h: f64, horiz_dist: f64) -> f64 {
            (15.0 * h).min(0.1 * horiz_dist)
        }

        (
            make_xl(self.heights.0, self.computed.horizon_distances.0),
            self.length - make_xl(self.heights.1, self.computed.horizon_distances.1),
        )
    }

    /// Given initial but naive horizon calculation, make some adjustments.
    ///
    /// This takes a hybrid approach of using the empirical formulae originally
    /// designed for the area mode to computed parameters for what is
    /// essentially the point-to-point mode used kinda like the area mode.
    ///
    /// This implementation does away completely with the distinction and
    /// eliminates all mode-switching logic. The original routine also had
    /// fallbacks for when a climate was not specified: we forbid that instead.
    ///
    /// See ITM sections `<43>`, `<45>`, `<46>`.
    fn adjust_horizons(&mut self, xl: (f64, f64)) {
        let mut q;
        let z;

        if self.computed.horizon_distances.0 + self.computed.horizon_distances.1 > 1.5 * self.length
        {
            // Redo light-of-sight horizons <45> if the path is line-of-sight

            let (_, nz) = z1sq1(self.computed.interval, self.elevations, xl);
            z = nz;
            self.computed.effective_heights = (
                // he = effective_heights
                self.heights.0 + fortran_dim(self.elevations[0], z.0),
                self.heights.1 + fortran_dim(self.elevations[self.elevations.len()], z.1),
            );

            fn make_dl(h: f64, curv: f64, terrain: f64) -> f64 {
                // curv = gme = effective_curvature
                (2.0 * h / curv).sqrt() * (-0.07 * (terrain / h.max(5.0)).sqrt()).exp()
            }

            self.computed.horizon_distances = (
                // dl = horizon_distances
                make_dl(
                    self.computed.effective_heights.0,
                    self.computed.effective_curvature,
                    self.computed.terrain_irregularity,
                ),
                make_dl(
                    self.computed.effective_heights.1,
                    self.computed.effective_curvature,
                    self.computed.terrain_irregularity,
                ),
            );

            q = self.computed.horizon_distances.0 + self.computed.horizon_distances.1;

            if q <= self.length {
                /* if there is a rounded horizon, or two obstructions, in the path */
                q = (self.length / q).powi(2);

                fn make_hedl(q: f64, he: f64, curv: f64, terrain: f64) -> (f64, f64) {
                    let he = he * q; /* tx effective height set to be path dist/self.computed.interval between obstacles */
                    (
                        he,
                        (2.0 * he / curv).sqrt() * (-0.07 * (terrain / he.max(5.0)).sqrt()).exp(),
                    )
                }

                let hedl = (
                    make_hedl(
                        q,
                        self.computed.effective_heights.0,
                        self.computed.effective_curvature,
                        self.computed.terrain_irregularity,
                    ),
                    make_hedl(
                        q,
                        self.computed.effective_heights.1,
                        self.computed.effective_curvature,
                        self.computed.terrain_irregularity,
                    ),
                );

                self.computed.effective_heights = ((hedl.0).0, (hedl.1).0); // he
                self.computed.horizon_distances = ((hedl.0).1, (hedl.1).1); // dl
            }

            /* original empirical adjustment?  uses delta-h to adjust grazing angles */
            fn make_qthe(he: f64, curv: f64, terrain: f64, horiz_dist: f64) -> f64 {
                let q = (2.0 * he / curv).sqrt();
                (0.65 * terrain * (q / horiz_dist - 1.0) - 2.0 * he) / q
            }

            self.computed.elevation_angles = (
                // the / theta
                make_qthe(
                    self.computed.effective_heights.0,
                    self.computed.effective_curvature,
                    self.computed.terrain_irregularity,
                    self.computed.horizon_distances.0,
                ),
                make_qthe(
                    self.computed.effective_heights.1,
                    self.computed.effective_curvature,
                    self.computed.terrain_irregularity,
                    self.computed.horizon_distances.1,
                ),
            );
        } else {
            // Get transhorizon effective heights <46>

            let (_, (z0, _)) = z1sq1(
                self.computed.interval,
                self.elevations,
                (xl.0, 0.9 * self.computed.horizon_distances.0),
            );
            let (_, (_, z1)) = z1sq1(
                self.computed.interval,
                self.elevations,
                (self.length - 0.9 * self.computed.horizon_distances.1, xl.1),
            );

            self.computed.effective_heights = (
                self.heights.0 + fortran_dim(self.elevations[0], z0),
                self.heights.1 + fortran_dim(self.elevations[self.elevations.len() - 1], z1),
            );
        }
    }

    pub fn attenuation_at(&self, distance_from_tx: f64) -> Result<f64, String> {
        if distance_from_tx < 0.0 {
            return Err("distance negative".into());
        }

        if distance_from_tx > self.length {
            return Err("distance above bounds".into());
        }

        Err("not implemented".into())
    }
}

/// Input settings for the model.
///
/// Refer to [ITU-R P.527] to derive ground permittivity and conductivity for
/// your region/terrain and frequency. Splash has related utilities.
/// Alternatively, you may use [SPLAT's simplified table][dielectrics].
///
/// [ITU-R P.527]: https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.527-4-201706-I!!PDF-E.pdf
/// [dielectrics]: http://www.qsl.net/n9zia/conduct.html
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Settings {
    /// Relative permittivity of the ground (aka "dielectric constant").
    pub permittivity: f64,

    /// Relative conductivity of the ground (in siemens per metre).
    pub conductivity: f64,

    /// Type of climate.
    pub climate: Climate,

    /// Surface refractivity reduced to sea level.
    pub surface_refractivity: f64,

    /// Frequency of modeled wave (MHz).
    pub frequency: f64,

    /// Polarisation of modeled wave.
    pub polarisation: Polarisation,

    // statistical fractions for the final quantisation
    pub conf: f64,
    pub rel: f64,
}

/// The inverse of the standard normal complementary probability function.
///
/// The standard normal complementary function is _Q(x) = 1 / √͞2͞π ∫ e^(-t²/2)_.
/// This inverse is the solution for _x_ to _q = Q(x)_, also noted _Q¯¹(q)_.
///
/// This function is used to scale the inputs (the desired fractions of time,
/// locations, situations to model) to later obtain normal quantiles.
///
/// The implementation is not the normal tables, but rather an approximation by
/// [Cecil Hastings][Hastings55], with a maximum error of 4.5 × 10¯⁴.
///
/// In the FORTRAN, this function was called `qerfi`
/// ("Q error function, inverted"), hence the constants. See <50>, <51>.
///
/// [Hastings55]: https://press.princeton.edu/titles/1133.html
fn inverse_normal_complementary(q: f64) -> f64 {
    let x = 0.5 - q;
    let mut t = (0.5 - x.abs()).max(0.000001);
    t = (-2.0 * t.ln()).sqrt();
    let v = t - ((QERFI_C.2 * t + QERFI_C.1) * t + QERFI_C.0)
        / (((QERFI_D.2 * t + QERFI_D.1) * t + QERFI_D.0) * t + 1.0);

    if x < 0.0 {
        -v
    } else {
        v
    }
}

/// C group of constants for the qerf/qerfi approximations.
const QERFI_C: (f64, f64, f64) = (2.515516698, 0.802853, 0.010328);

/// D group of constants for the qerf/qerfi approximations.
const QERFI_D: (f64, f64, f64) = (1.432788, 0.189269, 0.001308);


// below this point is "original" code that, once the legibility and rustification
// process is done, should completely disappear / not be used at all.

/// Point-to-Point propagation
pub fn point_to_point(
    distance: f64, // actually the interval
    elevations: &Vec<f64>,
    tx_height: f64,
    rx_height: f64,
    permittivity: f64,
    conductivity: f64,
    surfref: f64,
    freq: f64,
    climate: Climate,
    polarisation: Polarisation,
    conf: f64,
    rel: f64,
) -> Result<f64, ()> {
    let mut propa = PropA::default();

    let mut prop = Prop::default();
    prop.hg = (tx_height, rx_height);
    prop.ens = surfref;
    prop.wn = freq / 47.7;

    let mut propv = PropV::default();
    propv.klim = climate;
    propv.lvar = 5;
    propv.mdvar = 12;

    let subset = &elevations[1..(elevations.len() - 1)];
    let sum: f64 = subset.iter().sum();
    let zsys = sum / (subset.len() as f64);

    qlrps(zsys, polarisation, permittivity, conductivity, &mut prop);
    qlrpfl(
        distance,
        elevations,
        propv.mdvar,
        &mut prop,
        &mut propa,
        &mut propv,
    );

    let zc = inverse_normal_complementary(conf);
    let zr = inverse_normal_complementary(rel);
    let fs = 32.45 + 20.0 * (freq.log10() + (prop.dist / 1000.0).log10());
    let ret = avar(zr, 0.0, zc, &mut prop, &mut propv) + fs;

    if prop.kwx > 0 {
        Err(())
    } else {
        Ok(ret)
    }
}

// done (prepare_environment)
fn qlrps(zsys: f64, pol: Polarisation, dielect: f64, conduct: f64, prop: &mut Prop) {
    let gma = 157e-9;

    if zsys != 0.0 {
        prop.ens *= (-zsys / 9460.0).exp();
    }

    prop.gme = gma * (1.0 - 0.04665 * (prop.ens / 179.3).exp());

    let zq = Complex64::new(dielect, 376.62 * conduct / prop.wn);
    let zgnd = (zq - 1.0).sqrt();

    prop.zgnd = match pol {
        Polarisation::Horizontal => zgnd,
        Polarisation::Vertical => zgnd / zq,
    };
}

// mostly done (find_horizons, adjust_horizons)
fn qlrpfl(
    distance: f64,
    elevations: &Vec<f64>,
    mdvarx: isize,
    mut prop: &mut Prop,
    propa: &mut PropA,
    propv: &mut PropV,
) {
    prop.dist = elevations.len() as f64 * distance;
    let np = elevations.len();
    hzns(distance, &elevations, &mut prop);

    fn make_xl(hg: f64, dl: f64) -> f64 {
        (15.0 * hg).min(0.1 * dl)
    }

    let mut q;
    let z;
    let mut xl = (make_xl(prop.hg.0, prop.dl.0), make_xl(prop.hg.1, prop.dl.1));

    xl.1 = prop.dist - xl.1;
    prop.dh = d1thx(distance, elevations, xl);

    if prop.dl.0 + prop.dl.1 > 1.5 * prop.dist {
        let (_, nz) = z1sq1(distance, elevations, xl);
        z = nz;
        prop.he = (
            prop.hg.0 + fortran_dim(elevations[0], z.0),
            prop.hg.1 + fortran_dim(elevations[np], z.1),
        );

        fn make_dl(he: f64, gme: f64, dh: f64) -> f64 {
            (2.0 * he / gme).sqrt() * (-0.07 * (dh / he.max(5.0)).sqrt()).exp()
        }

        prop.dl = (
            make_dl(prop.he.0, prop.gme, prop.dh),
            make_dl(prop.he.1, prop.gme, prop.dh),
        );

        q = prop.dl.0 + prop.dl.1;

        if q <= prop.dist {
            /* if there is a rounded horizon, or two obstructions, in the path */
            let temp = prop.dist / q;
            q = temp * temp;

            fn make_hedl(q: f64, he: f64, gme: f64, dh: f64) -> (f64, f64) {
                let he = he * q; /* tx effective height set to be path dist/distance between obstacles */
                (
                    he,
                    (2.0 * he / gme).sqrt() * (-0.07 * (dh / he.max(5.0)).sqrt()).exp(),
                )
            }

            let hedl = (
                make_hedl(q, prop.he.0, prop.gme, prop.dh),
                make_hedl(q, prop.he.1, prop.gme, prop.dh),
            );

            prop.he = ((hedl.0).0, (hedl.1).0);
            prop.dl = ((hedl.0).1, (hedl.1).1);
        }

        /* original empirical adjustment?  uses delta-h to adjust grazing angles */
        fn make_qthe(he: f64, gme: f64, dh: f64, dl: f64) -> f64 {
            let q = (2.0 * he / gme).sqrt();
            (0.65 * dh * (q / dl - 1.0) - 2.0 * he) / q
        }

        prop.the = (
            make_qthe(prop.he.0, prop.gme, prop.dh, prop.dl.0),
            make_qthe(prop.he.1, prop.gme, prop.dh, prop.dl.1),
        );
    } else {
        let (_, (z0, _)) = z1sq1(distance, elevations, (xl.0, 0.9 * prop.dl.0));
        let (_, (_, z1)) = z1sq1(distance, elevations, (prop.dist - 0.9 * prop.dl.1, xl.1));

        prop.he = (
            prop.hg.0 + fortran_dim(elevations[0], z0),
            prop.hg.1 + fortran_dim(elevations[np - 1], z1),
        );
    }

    propv.lvar = propv.lvar.max(3);

    if mdvarx >= 0 {
        propv.mdvar = mdvarx;
        propv.lvar = propv.lvar.max(4);
    }

    propv.lvar = 5;

    // the big compute
    lrprop(0.0, prop, propa);
}

// done (find_horizons)
fn hzns(distance: f64, elevations: &Vec<f64>, prop: &mut Prop) {
    let np = elevations.len();
    let xi = distance;

    let za = elevations[0] + prop.hg.0;
    let zb = elevations[np - 1] + prop.hg.1;

    let qc = 0.5 * prop.gme;
    let mut q = qc * prop.dist;

    prop.the.1 = (zb - za) / prop.dist;
    prop.the.0 = prop.the.1 - q;
    prop.the.1 = -prop.the.1 - q;

    prop.dl = (prop.dist, prop.dist);

    if np >= 2 {
        let mut sa = 0.0;
        let mut sb = prop.dist;
        let mut wq = true;

        for i in 1..np {
            sa += xi;
            sb -= xi;
            q = elevations[i] - (qc * sa + prop.the.0) * sa - za;

            if q > 0.0 {
                prop.the.0 += q / sa;
                prop.dl.0 = sa;
                wq = false;
            }

            if !wq {
                q = elevations[i] - (qc * sb + prop.the.1) * sb - zb;

                if q > 0.0 {
                    prop.the.1 += q / sb;
                    prop.dl.1 = sb;
                }
            }
        }
    }
}

// <48> ("delta h over x", the interdecile range of elevations between x1 & x2)
fn d1thx(distance: f64, elevations: &Vec<f64>, xl: (f64, f64)) -> f64 {
    let np = elevations.len();
    let mut xa = xl.0 / distance;
    let xb = xl.1 / distance;

    if (xb - xa) < 2.0 {
        return 0.0;
    }

    let mut ka = (0.1 * (xb - xa + 8.0)) as usize;
    ka = 4.max(ka).min(25);

    let n = 10 * ka - 5;
    let kb = n - ka + 1;
    let sn = n - 1;
    let mut s = Vec::with_capacity(n);

    let xb = (xb - xa) / (sn as f64);
    let mut k = (xa + 1.0) as usize;
    xa -= k as f64;

    for j in 0..n {
        while xa > 0.0 && k < np {
            xa -= 1.0;
            k += 1;
        }

        s[j] = elevations[k] + (elevations[k] - elevations[k - 1]) * xa;
        xa = xa + xb;
    }

    let ((_, sn), (mut xa, mut xb)) = z1sq1(1.0, &s, (0.0, sn as f64));
    xb = (xb - xa) / (sn as f64);

    for j in 0..n {
        s[j] -= xa;
        xa = xa + xb;
    }

    let d1thxv = qtile(n - 1, &mut s, ka - 1) - qtile(n - 1, &mut s, kb - 1);
    d1thxv / (1.0 - 0.8 * (-(xl.1 - xl.0) / 50.0e3).exp())
}

// transliteration typo of zlsql
fn z1sq1(interval: f64, elevations: &Vec<f64>, x: (f64, f64)) -> ((f64, f64), (f64, f64)) {
    zlsql(interval, elevations, x)
}

// <53> done in least_squares_linear_fit
fn zlsql(interval: f64, elevations: &Vec<f64>, x: (f64, f64)) -> ((f64, f64), (f64, f64)) {
    (x, least_squares_linear_fit(interval, elevations, x))
}

// <52>. provides a quantile
fn qtile(nn: usize, elevations: &mut Vec<f64>, ir: usize) -> f64 {
    let mut m = 0;
    let mut n = nn;
    let k = 0.max(ir).min(n);

    let mut q = 0.0;
    let mut r: f64;
    let mut j: usize;
    let mut j1 = 0;
    let mut i0 = 0;

    let mut i;
    let mut goto10 = true;
    loop {
        if goto10 {
            q = elevations[k];
            i0 = m;
            j1 = n;
        }

        i = i0;

        while i <= n && elevations[i] >= q {
            i += 1;
        }

        if i > n {
            i = n;
        }

        j = j1;

        while j >= m && elevations[j] <= q {
            j -= 1;
        }

        if j < m {
            j = m;
        }

        if i < j {
            r = elevations[i];
            elevations[i] = elevations[j];
            elevations[j] = r;
            i0 = i + 1;
            j1 = j - 1;
            goto10 = false;
        } else if i < k {
            elevations[k] = elevations[i];
            elevations[i] = q;
            m = i + 1;
            goto10 = true;
        } else if j > k {
            elevations[k] = elevations[j];
            elevations[j] = q;
            n = j - 1;
            goto10 = true;
        } else {
            break;
        }
    }

    return q;
}

// <4>, <5>, the actual propagation compute
fn lrprop(d: f64, prop: &mut Prop, propa: &mut PropA) {
    // first part is mostly input checks
    // "kwx" is the error output... higher is worse, but really anything
    // besides zero is bad, so during rewrite just abort early anytime.

    if !prop.setup_done {
        // --- <6> --- prep secondary params

        propa.dls.0 = (2.0 * prop.he.0 / prop.gme).sqrt();
        propa.dls.1 = (2.0 * prop.he.1 / prop.gme).sqrt();

        propa.dlsa = propa.dls.0 + propa.dls.1;
        propa.dla = prop.dl.0 + prop.dl.1;
        propa.tha = (prop.the.0 + prop.the.1).max(-propa.dla * prop.gme);

        // --- <7> --- checks ranges

        if prop.wn < 0.838 || prop.wn > 210.0 {
            prop.kwx = prop.kwx.max(1);
        }

        fn make_kwx_hg(kwx: usize, hg: f64) -> usize {
            if hg < 1.0 || hg > 1000.0 {
                kwx.max(1)
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_hg(prop.kwx, prop.hg.0);
        prop.kwx = make_kwx_hg(prop.kwx, prop.hg.1);

        fn make_kwx_dl(kwx: usize, the: f64, dl: f64, dls: f64) -> usize {
            if (the.abs() > 200e-3) || (dl < 0.1 * dls) || (dl > 3.0 * dls) {
                kwx.max(3)
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_dl(prop.kwx, prop.the.0, prop.dl.0, propa.dls.0);
        prop.kwx = make_kwx_dl(prop.kwx, prop.the.1, prop.dl.1, propa.dls.1);

        if (prop.ens < 250.0)
            || (prop.ens > 400.0)
            || (prop.gme < 75e-9)
            || (prop.gme > 250e-9)
            || (prop.zgnd.re <= prop.zgnd.im.abs())
            || (prop.wn < 0.419)
            || (prop.wn > 420.0)
        {
            prop.kwx = 4; // fail here
        }

        fn make_kwx_hg_again(kwx: usize, hg: f64) -> usize {
            if hg < 0.5 || hg > 3000.0 {
                4 // fail here
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_hg_again(prop.kwx, prop.hg.0);
        prop.kwx = make_kwx_hg_again(prop.kwx, prop.hg.1);

        // --- <9> --- diffraction coefficients --- see T.A. 4.2 through 4.8
        prop.dmin = (prop.he.0 - prop.he.1).abs() / 200e-3;
        let q = adiff(0.0, prop, propa);
        prop.xae = (prop.wn * prop.gme.powi(2)).cbrt();
        let d3 = propa.dlsa.max(1.3787 * prop.xae + propa.dla);
        let d4 = d3 + 2.7574 * prop.xae;
        let a3 = adiff(d3, prop, propa);
        let a4 = adiff(d4, prop, propa);
        propa.emd = (a4 - a3) / (d4 - d3);
        propa.aed = a3 - propa.emd * d3;
    }

    prop.dist = d;

    // <8> distance bounds checks
    if prop.dist > 0.0 {
        if prop.dist > 1000e3 {
            prop.kwx = prop.kwx.max(1);
        }

        if prop.dist < prop.dmin {
            prop.kwx = prop.kwx.max(3);
        }

        if prop.dist < 1e3 || prop.dist > 2000e3 {
            prop.kwx = 4; // fail here
        }
    }

    // <15> line of sight calculations
    if prop.dist < propa.dlsa {
        if !prop.wlos { // <16> prep constants on first run
            alos(0.0, prop, propa);
            let d2 = propa.dlsa;
            let a2 = propa.aed + d2 * propa.emd;
            let mut d0 = 1.908 * prop.wn * prop.he.0 * prop.he.1; // T.A. 4.38

            let d1;
            if propa.aed >= 0.0 {
                d0 = d0.min(0.5 * propa.dla); // T.A. 4.28
                d1 = d0 + 0.25 * (propa.dla - d0); // T.A. 4.29
            } else {
                d1 = (-propa.aed / propa.emd).max(0.25 * propa.dla); // T.A. 4.30
            }

            let a1 = alos(d1, prop, propa); // T.A. 4.31

            if d0 < d1 {
                let a0 = alos(d0, prop, propa); // T.A. 4.30
                let q = (d2 / d0).ln();
                propa.ak2 = 0.0f64.max(
                    ((d2 - d0) * (a1 - a0) - (d1 - d0) * (a2 - a0))
                        / ((d2 - d0) * (d1 / d0).ln() - (d1 - d0) * q),
                ); // T.A. 4.32

                let wq = propa.aed >= 0.0 || propa.ak2 > 0.0;

                if wq {
                    propa.ak1 = (a2 - a0 - propa.ak2 * q) / (d2 - d0); // T.A. 4.33

                    if propa.ak1 < 0.0 {
                        propa.ak1 = 0.0; // T.A. 4.36
                        propa.ak2 = fortran_dim(a2, a0) / q; // T.A. 4.35

                        // T.A. 4.37
                        if propa.ak2 == 0.0 {
                            propa.ak1 = propa.emd;
                        }
                    }
                } else {
                    propa.ak1 = (a2 - a1) / (d2 - d1); // T.A. 4.40
                    propa.ak2 = 0.0; // T.A. 4.41

                    // T.A. 4.37
                    if propa.ak1 <= 0.0 {
                        propa.ak1 = propa.emd;
                    }
                }
            } else {
                // same as above
                propa.ak1 = (a2 - a1) / (d2 - d1);
                propa.ak2 = 0.0;

                if propa.ak1 <= 0.0 {
                    propa.ak1 = propa.emd;
                }
            }

            // T.A. 4.42
            propa.ael = a2 - propa.ak1 * d2 - propa.ak2 * d2.ln();

            prop.wlos = true;
        }

        // Do calculation when given real distance (dist = 0 is constant prep)
        if prop.dist > 0.0 {
            // T.A. 4.1
            prop.aref = propa.ael + propa.ak1 * prop.dist + propa.ak2 * prop.dist.ln();
        }
    }

    // <20> troposcatter calculations
    if prop.dist <= 0.0 || prop.dist >= propa.dlsa {
        if !prop.wscat { // <21> -- setup constants
            ascat(0.0, prop, propa);
            let d5 = propa.dla + 200e3; // T.A. 4.52
            let d6 = d5 + 200e3; // T.A. 4.53
            let a6 = ascat(d6, prop, propa); // T.A. 4.54
            let a5 = ascat(d5, prop, propa); // T.A. 4.55

            if a5 < 1000.0 {
                propa.ems = (a6 - a5) / 200e3; // T.A. 4.57
                propa.dx = propa.dlsa.max(
                    (propa.dla + 0.3 * prop.xae * (47.7 * prop.wn).ln())
                        .max((a5 - propa.aed - propa.ems * d5) / (propa.emd - propa.ems)),
                ); // T.A. 4.58
                propa.aes = (propa.emd - propa.ems) * propa.dx + propa.aed; // T.A. 4.59
            } else {
                propa.ems = propa.emd;
                propa.aes = propa.aed;
                propa.dx = 10.0e6; // T.A. 4.56
            }

            prop.wscat = true;
        }

        // T.A. 4.1
        if prop.dist > propa.dx {
            prop.aref = propa.aes + propa.ems * prop.dist;
        } else {
            prop.aref = propa.aed + propa.emd * prop.dist;
        }
    }

    prop.aref = prop.aref.max(0.0);
}

// <10> diffraction attenuation at distance d from site
fn adiff(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
    // first call with d == 0.0 is used to setup constants
    // should be extracted into two functions and the constants stored in a
    // struct or the cached structure or something.
    if d == 0.0 { // see <11>
        let mut q = prop.hg.1 * prop.hg.1;
        let qk = prop.he.0 * prop.he.1 - q;

        // if prop.mode == Mode::PointToPoint {
        //     q += 10.0;
        // }

        // "parts of Q" see T.A. 4.9
        let wd1 = (1.0 + qk / q).sqrt();
        let xd1 = propa.dla + propa.tha / prop.gme;

        // T.A. 4.10
        q = (1.0 - 0.8 * (-propa.dlsa / 50e3).exp()) * prop.dh;
        q *= 0.78 * (-(q / 16.0).powf(0.25)).exp();
        let afo = 15.0f64.min(2.171 * (1.0 + 4.77e-4 * prop.hg.0 * prop.hg.1 * prop.wn * q).ln());

        // T.A. 6.7
        let qk = 1.0 / prop.zgnd.norm_sqr().sqrt();
        let mut aht = 20.0;
        let mut xht = 0.0;

        if false {
            /// ??? probably used somehow, needs a reread
            fn make_axht(dl: f64, he: f64, wn: f64, qk: f64) -> (f64, f64) {
                let a = 0.5 * dl.powi(2) / he;
                let wa = (a * wn).cbrt();
                let pk = qk / wa;
                let q = (1.607 - pk) * 151.0 * wa * dl / a;
                (q, fht(q, pk))
            }

            let (x, a) = make_axht(prop.dl.0, prop.he.0, prop.wn, qk);
            xht += x;
            aht += a;

            let (x, a) = make_axht(prop.dl.1, prop.he.1, prop.wn, qk);
            xht += x;
            aht += a;
        }

        0.0 // returns 0 just because this is the dummy setup round
    } else { // see <12>

        // T.A. 4.12
        let th = propa.tha + d * prop.gme;
        let ds = d - propa.dla;
        let mut q = 0.0795775 * prop.wn * ds * th * th;

        // T.A. 4.14
        let adiffv =
            aknfe(q * prop.dl.0 / (ds + prop.dl.0)) + aknfe(q * prop.dl.1 / (ds + prop.dl.1));

        // Dummy values to get this compiling -- they're from the constants run
        let qk = 0.0;
        let xht = 0.0;
        let aht = 0.0;
        let wd1 = 0.0;
        let xd1 = 0.0;
        let afo = 0.0;

        // T.A. 4.16
        let a = ds / th;
        let wa = (a * prop.wn).cbrt();
        let pk = qk / wa; // T.A. 4.17
        q = (1.607 - pk) * 151.0 * wa * th + xht; // T.A. 4.18 and 6.2
        let ar = 0.05751 * q - 4.343 * q.ln() - aht; // T.A. 4.20
        q = (wd1 + xd1 / d) * 6283.2f64.min((1.0 - 0.8 * (-d / 50e3).exp()) * prop.dh * prop.wn);

        // T.A. 4.9
        let wd = 25.1 / (25.1 + q.sqrt());

        // T.A. 4.11
        ar * wd + (1.0 - wd) * adiffv + afo
    }
}

// <17> line of sight attenuation at distance d from site
fn alos(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
    let mut wls = 0.0;

    // see adiff comment on constant gen and splitting
    if d == 0.0 { // <18>
        // T.A. 4.43
        wls = 0.021 / (0.021 + prop.wn * prop.dh / 10e3f64.max(propa.dlsa));
        0.0
    } else { // <19>
        let mut q = (1.0 - 0.8 * (-d / 50e3).exp()) * prop.dh;
        let s = 0.78 * q * (-(q / 16.0).powf(0.25)).exp();
        q = prop.he.0 + prop.he.1;
        let sps = q / (d.powi(2) + q.powi(2)).sqrt();

        // T.A. 4.47
        let mut r = (sps - prop.zgnd) / (sps + prop.zgnd) * (-10.0f64.min(prop.wn * s * sps)).exp();
        q = r.norm_sqr();

        // T.A. 4.48
        if q < 0.25 || q < sps {
            r = r * (sps / q).sqrt();
        }

        let alosv = propa.emd * d + propa.aed; // T.A. 4.45
        q = prop.wn * prop.he.0 * prop.he.1 * 2.0 / d; // T.A. 4.49

        // T.A. 4.50
        if q > 1.57 {
            q = 3.14 - 2.4649 / q;
        }

        // T.A. 4.51 and 4.44
        let qq = Complex64::new(q.cos(), -q.sin());
        (-4.343 * (qq + r).norm_sqr().ln() - alosv) * wls + alosv
    }
}

// <22> scatter attenuation at distance d from site
// See TN101 for approximation method description
fn ascat(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
    // static double ad, rr, etq, h0s;
    // double h0, r1, r2, z0, ss, et, ett, th, q;
    // double ascatv, temp;

    // see adiff comment on constant gen and splitting
    if d == 0.0 { // <23>
        prop.ad = prop.dl.0 - prop.dl.1;
        prop.rr = prop.he.1 / prop.rch.0;

        if prop.ad < 0.0 {
            prop.ad = -prop.ad;
            prop.rr = 1.0 / prop.rr;
        }

        // T.A. 4.67 (partial)
        prop.etq = (5.67e-6 * prop.ens - 2.32e-3) * prop.ens + 0.031;
        prop.h0s = -15.0;
        0.0
    } else { // <24>
        let mut h0;
        if prop.h0s > 15.0 {
            h0 = prop.h0s;
        } else {
            let th = prop.the.0 + prop.the.1 + d * prop.gme; // T.A. 4.61

            // T.A. 4.62
            let mut r2 = 2.0 * prop.wn * th;
            let r1 = r2 * prop.he.0;
            r2 *= prop.he.1;

            // bounds check. 1001 is "error" value to exit out
            if r1 < 0.2 && r2 < 0.2 {
                return 1001.0;
            }

            let mut ss = (d - prop.ad) / (d + prop.ad); // T.A. 4.65

            // T.A. 4.66
            let mut q = prop.rr / ss;
            ss = ss.max(0.1);
            q = q.max(0.1).min(10.0);
            let z0 = (d - prop.ad) * (d + prop.ad) * th * 0.25 / d;

            // T.A. 4.67
            let temp = (z0 / 8.0e3).min(1.7).powi(6);
            let et = (prop.etq * (-temp).exp() + 1.0) * z0 / 1.7556e3;
            let ett = et.max(1.0);

            h0 = (h0f(r1, ett) + h0f(r2, ett)) * 0.5; // T.A. 6.12
            h0 += 1.38 - ett.ln().min(h0) * ss.ln() * q.ln() * 0.49; // T.A. 6.10 and 6.11
            h0 = fortran_dim(h0, 0.0);

            // T.A. 6.14
            if et < 1.0 {
                let temp = (1.0 + 1.4142 / r1) * (1.0 + 1.4142 / r2);
                h0 = et * h0
                    + (1.0 - et) * 4.343 * (temp.powi(2) * (r1 + r2) / (r1 + r2 + 2.8284)).ln();
            }

            // calc got out of bounds, revert back
            if h0 > 15.0 && prop.h0s >= 0.0 {
                h0 = prop.h0s;
            }
        }

        prop.h0s = h0;

        // T.A. 4.60
        let th = propa.tha + d * prop.gme;

        // T.A. 4.63 and 6.8
        ahd(th * d) + 4.343 * (47.7 * prop.wn * th.powi(4)).ln()
            - 0.1 * (prop.ens - 301.0) * (-th * d / 40e3).exp() + h0
    }
}

// <13> attenuation on a single knife edge
// this is an approximation of a Fresnel integral, see T.A. 6.1
fn aknfe(v2: f64) -> f64 {
    if v2 < 5.76 {
        6.02 + 9.11 * v2.sqrt() - 1.27 * v2
    } else {
        12.953 + 10.0 * v2.log10()
    }
}

// <14> Height gain over geodesic model -- here instead approximated to a
// smooth spherical earth. See T.A. 6.4
fn fht(x: f64, pk: f64) -> f64 {
    if x < 200.0 {
        let w = -pk.ln();

        if pk < 1.0e-5 || x * w.powi(3) > 5495.0 {
            if x > 1.0 {
                // this is changed from the original! to be investigated
                // <14> has 40.0 as 17.372 ref T.A. 6.5
                40.0 * x.log10() - 117.0
            } else {
                -117.0
            }
        } else {
            // T.A. 6.6
            2.5e-5 * x.powi(2) / pk - 8.686 * w - 15.0
        }
    } else {
        // T.A. 6.3
        let fhtv = 0.05751 * x - 10.0 * x.log10();

        if x < 2000.0 {
            let w = 0.0134 * x * (-0.005 * x).exp();
            (1.0 - w) * fhtv + w * (40.0 * x.log10() - 117.0) // T.A. 6.4
        } else {
            fhtv
        }
    }
}

// <25> H01 function for scatter fields, see T.A. §6
fn h0f(r: f64, et: f64) -> f64 {
    let a = [25.0, 80.0, 177.0, 395.0, 705.0];
    let b = [24.0, 45.0, 68.0, 80.0, 105.0];

    let (it, q) = if et <= 0.0 {
        (1, 0.0)
    } else if et >= 5.0 {
        (5, 0.0)
    } else {
        (et as usize, et.fract())
    };

    let x = (1.0 / r).powi(2);
    let h0fv = 4.343 * ((a[it - 1] * x + b[it - 1]) * x + 1.0).ln(); // T.A. 6.13

    if q == 0.0 {
        h0fv
    } else {
        (1.0 - q) * h0fv + q * 4.343 * ((a[it] * x + b[it]) * x + 1.0).ln()
    }
}

// <26> F(theta d) function for scatter fields
fn ahd(td: f64) -> f64 {
    let a = [133.4, 104.6, 71.8];
    let b = [0.332e-3, 0.212e-3, 0.157e-3];
    let c = [-4.343, -1.086, 2.171];

    // choice of constants
    let i = if td <= 10e3 {
        0
    } else if td <= 70e3 {
        1
    } else {
        2
    };

    a[i] + b[i] * td + c[i] * td.ln() // T.A. 6.9
}

const RT: f64 = 7.8;
const RL: f64 = 24.0;

// <27> "the statistics"
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct PropV { // used in avar
    pub sgc: f64, // stddev of confidence -- an output of avar

    pub lvar: isize, // control switch
    // this is an optimisation in the original program and will need to be
    // rewritten out, with the different sections it controls split out.
    // the idea is to run some preparation functions only when needed, but
    // it makes the whole thing incomprehensible. See <28> for more.

    pub mdvar: isize, // variability mode switch
    pub klim: Climate,
}

// <28>
fn avar(zzt: f64, zzl: f64, zzc: f64, prop: &mut Prop, propv: &mut PropV) -> f64 {
    // static int kdv;
    // static bool ws, w1;
    let mut kdv = 0;
    let mut ws = false;
    let mut w1 = false;

    // static double dexa, de, vs0, sgl, sgtm, sgtp, sgtd, tgtd
    let mut dexa = 0.0;
    let mut de = 0.0;
    let mut vmd = 0.0;
    let mut vs0 = 0.0;
    let mut sgl = 0.0;
    let mut sgtm = 0.0;
    let mut sgtp = 0.0;
    let mut sgtd = 0.0;
    let mut tgtd = 0.0;
    // ^ the "constants" to be set up (also see lvar) <27>

    // <29>, <30>, <32> select the set of constants to be used for climate adjustments
    let cc: ClimateConstants = propv.klim.into();

    if propv.lvar > 0 { // <31>
        match propv.lvar {
            4 => { // <33>
                kdv = propv.mdvar;
                ws = kdv >= 20;
                if ws {
                    kdv -= 20;
                }

                w1 = kdv >= 10;
                if w1 {
                    kdv -= 10;
                }

                if kdv < 0 || kdv > 3 {
                    kdv = 0;
                    prop.kwx = prop.kwx.max(2);
                }
            }
            2 => { // <35> system
                dexa = (18e6 * prop.he.0).sqrt()
                    + (18e6 * prop.he.1).sqrt()
                    + (575.7e12 / prop.wn).cbrt();
            }
            1 => { // <36> distance
                de = if prop.dist < dexa {
                    130e3 * prop.dist / dexa
                } else {
                    130e3 + prop.dist - dexa
                };
            }
            _ => {}
        }

        // <32> climate
        let refs = cc.reference_values(de, prop.wn);
        vmd = refs.0;
        sgtm = refs.1;
        sgtp = refs.2;
        sgtd = refs.3;
        tgtd = refs.4;

        // <36> distance again
        sgl = if w1 {
            0.0
        } else {
            let q = (1.0 - 0.8 * (-prop.dist / 50e3).exp()) * prop.dh * prop.wn;
            10.0 * q / (q + 13.0)
        };

        // <36> still distance
        vs0 = if ws {
            0.0
        } else {
            (5.0 + 3.0 * (-de / 100e3).exp()).powi(2)
        };

        propv.lvar = 0;
    }

    let mut zt = zzt;
    let mut zl = zzl;

    // <37> normal deviates
    match kdv {
        0 => {
            zt = zzc;
            zl = zzc;
        }
        1 => {
            zl = zzc;
        }
        2 => {
            zl = zt;
        }
        _ => {}
    };

    // <37>
    if zt.abs() > 3.1 || zl.abs() > 3.1 || zzc.abs() > 3.1 {
        prop.kwx = prop.kwx.max(1);
    }

    // <38> resolve standard deviations
    let sgt = if zt < 0.0 {
        sgtm
    } else if zt <= cc.zd {
        sgtp
    } else {
        sgtd + tgtd / zt
    };
    let vs = vs0 + (sgt * zt).powi(2) / (RT + zzc * zzc) + (sgl * zl).powi(2) / (RL + zzc * zzc);

    // <39> resolve deviations yr, yc
    let (yr, sgc) = match kdv {
        0 => (0.0, sgt.powi(2) + sgl.powi(2) + vs),
        1 => (sgt * zt, sgl.powi(2) + vs),
        2 => ((sgt.powi(2) + sgl.powi(2)).sqrt() * zt, vs),
        _ => (sgt * zt + sgl * zl, vs),
    };
    propv.sgc = sgc.sqrt();

    // T.A. 5.1
    let avarv = prop.aref - vmd - yr - propv.sgc * zzc;
    if avarv < 0.0 {
        avarv * (29.0 - avarv) / (29.0 - 10.0 * avarv) // T.A. 5.2
    } else {
        avarv
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Prop {
    /// Reference attenuation
    pub aref: f64,

    /// Distance from tx to rx
    pub dist: f64,

    /// Antenna structural heights (tx, rx)
    pub hg: (f64, f64),

    pub rch: (f64, f64),

    /// Wave number (radio frequency)
    pub wn: f64,

    /// Terrain irregularity parameter
    pub dh: f64,

    pub dhd: f64,

    /// Surface refractivity
    pub ens: f64,

    pub encc: f64,
    pub cch: f64,
    pub cd: f64,

    /// Earth's effective curvature
    pub gme: f64,

    /// Surface transfer impedance to the ground
    pub zgnd: Complex64,

    /// Antenna effective heights (tx, rx)
    pub he: (f64, f64),

    /// Horizon distances (tx, rx)
    pub dl: (f64, f64),

    /// Horizon elevation angles (tx, rx)
    pub the: (f64, f64),

    pub tiw: f64,
    pub ght: f64,
    pub ghr: f64,
    pub rph: f64,
    pub hht: f64,
    pub hhr: f64,
    pub tgh: f64,
    pub tsgh: f64,
    pub thera: f64,
    pub thenr: f64,
    pub rpl: isize,

    /// Error indicator
    pub kwx: usize,

    /// Whether initial setup has been done.
    ///
    /// Usually ITM will be called to get a sequence of results for varying
    /// distances, and this switch gets set after the first iteration so common
    /// parameters are only computed once.
    pub setup_done: bool,

    pub ptx: isize,
    pub los: isize,

    // statics below
    pub dmin: f64,
    pub xae: f64,
    pub wlos: bool,
    pub wscat: bool,

    pub ad: f64,
    pub rr: f64,
    pub etq: f64,
    pub h0s: f64,
}

/// Secondary parameters computed in LRProp.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct PropA {
    /// Line of sight distance
    pub dlsa: f64,

    /// Scatter distance
    pub dx: f64,

    /// Line of sight coefficient L
    pub ael: f64,

    /// Line of sight coefficient 1
    pub ak1: f64,

    /// Line of sight coefficient 2
    pub ak2: f64,

    /// Diffraction coefficient 1
    pub aed: f64,

    /// Diffraction coefficient 2
    pub emd: f64,

    /// Scatter coefficient 1
    pub aes: f64,

    /// Scatter coefficient 2
    pub ems: f64,

    /// Smooth earth horizon distances (tx, rx)
    pub dls: (f64, f64),

    /// Total horizon distance
    pub dla: f64,

    /// Total bending angle
    pub tha: f64,
}
