// #![deny(missing_docs)]

//! Climate enum, constants, and Long Term Fading calculations.
//!
//! TODO: LTF explainer.

/// Radio-climate regions.
///
/// TODO: explainer.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Climate {
    Equatorial,
    ContinentalSubtropical,
    MaritimeSubtropical,
    Desert,
    ContinentalTemperate,
    MaritimeTemperateOverLand,
    MaritimeTemperateOverSea,
}

impl Default for Climate {
    fn default() -> Self {
        Climate::ContinentalTemperate
    }
}

/// Climate constants for Long Term Fading calculations.
///
/// Instead of just having a bunch of arrays to hold these, they are grouped
/// in meaningful bunches (although the names remain inscrutable, their usage
/// clearly indicates this arrangement holds some kind of meaning), and one
/// struct holds all the constants for one climate.
///
/// To get the constants for a climate, use:
///
/// ```
/// # use splash::itm::climate::{Climate, ClimateConstants};
/// # let climate = Climate::default();
/// let cc: ClimateConstants = climate.into();
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct ClimateConstants {
    pub name: Climate,

    pub bv: (f64, f64),
    pub xv: (f64, f64, f64),
    pub bsm: (f64, f64),
    pub xsm: (f64, f64, f64),
    pub bsp: (f64, f64),
    pub xsp: (f64, f64, f64),

    /// CD (Table 5.1 of T.A.)
    pub cd: f64,
    /// ZD (Table 5.1 of T.A.)
    pub zd: f64,

    pub bfm: (f64, f64, f64),
    pub bfp: (f64, f64, f64),
}

impl ClimateConstants {
    /// Creates from contants.
    fn new(
        name: Climate,
        bv: (f64, f64),
        xv: (f64, f64, f64),
        bsm: (f64, f64),
        xsm: (f64, f64, f64),
        bsp: (f64, f64),
        xsp: (f64, f64, f64),
        cd: f64,
        zd: f64,
        bfm: (f64, f64, f64),
        bfp: (f64, f64, f64),
    ) -> Self {
        Self {
            name,
            bv,
            xv,
            bsm,
            xsm,
            bsp,
            xsp,
            cd,
            zd,
            bfm,
            bfp,
        }
    }

    /// Computes long-term fading reference values from climate constants.
    pub fn reference_values(self, de: f64, wn: f64) -> (f64, f64, f64, f64, f64) {
        let q = (0.133 * wn).ln();
        let gm = self.bfm.0 + self.bfm.1 / ((self.bfm.2 * q).powi(2) + 1.0);
        let gp = self.bfp.0 + self.bfp.1 / ((self.bfp.2 * q).powi(2) + 1.0);

        let vmd = Self::curve(self.bv, self.xv, de);
        let sgtm = Self::curve(self.bsm, self.xsm, de) * gm;
        let sgtp = Self::curve(self.bsp, self.xsp, de) * gp;

        let sgtd = sgtp * self.cd;
        let tgtd = (sgtp - sgtd) * self.zd;

        (vmd, sgtm, sgtp, sgtd, tgtd)
    }

    /// Long-term fading climate curves.
    ///
    /// This function's only reference is in the FORTRAN source. No comment is given
    /// as to how it was derived, and whether the figure in the research are from
    /// this function, or whether the function is fit from the figure.
    ///
    /// The figure is available in [Technical Note 101 Volume I][TN101-I] and
    /// [Volume II][TN101-II], sections 10 and III respectively.
    ///
    /// For more background details see the `avar` function documentation.
    ///
    /// [TN101-I]: https://www.its.bldrdoc.gov/publications/2726.aspx
    /// [TN101-II]: https://www.its.bldrdoc.gov/publications/2727.aspx
    pub fn curve(b: (f64, f64), x: (f64, f64, f64), de: f64) -> f64 {
        (b.0 + b.1 / (1.0 + ((de - x.1) / x.2).powi(2))) * (de / x.0).powi(2)
            / (1.0 + (de / x.0).powi(2))
    }
}

impl From<Climate> for ClimateConstants {
    fn from(climate: Climate) -> Self {
        match climate {
            name @ Climate::Equatorial => Self::new(
                name,
                (-9.67, 12.7),
                (144.9e3, 190.3e3, 133.8e3),
                (2.13, 159.5),
                (762.2e3, 123.6e3, 94.5e3),
                (2.11, 102.3),
                (636.9e3, 134.8e3, 95.6e3),
                1.224,
                1.282,
                (1.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ),
            name @ Climate::ContinentalSubtropical => Self::new(
                name,
                (-0.62, 9.19),
                (228.9e3, 205.2e3, 143.6e3),
                (2.66, 7.67),
                (100.4e3, 172.5e3, 136.4e3),
                (6.87, 15.53),
                (138.7e3, 143.7e3, 98.6e3),
                0.801,
                2.161,
                (1.0, 0.0, 0.0),
                (0.93, 0.31, 2.00),
            ),
            name @ Climate::MaritimeSubtropical => Self::new(
                name,
                (1.26, 15.5),
                (262.6e3, 185.2e3, 99.8e3),
                (6.11, 6.65),
                (138.2e3, 242.2e3, 178.6e3),
                (10.08, 9.60),
                (165.3e3, 225.7e3, 129.7e3),
                1.380,
                1.282,
                (1.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ),
            name @ Climate::Desert => Self::new(
                name,
                (-9.21, 9.05),
                (84.1e3, 101.1e3, 98.6e3),
                (1.98, 13.11),
                (139.1e3, 132.7e3, 193.5e3),
                (3.68, 159.3),
                (464.4e3, 93.1e3, 94.2e3),
                1.000,
                20.0,
                (1.0, 0.0, 0.0),
                (0.93, 0.19, 1.79),
            ),
            name @ Climate::ContinentalTemperate => Self::new(
                name,
                (-0.62, 9.19),
                (228.9e3, 205.2e3, 143.6e3),
                (2.68, 7.16),
                (93.7e3, 186.8e3, 133.5e3),
                (4.75, 8.12),
                (93.2e3, 135.9e3, 113.4e3),
                1.224,
                1.282,
                (0.92, 0.25, 1.77),
                (0.93, 0.31, 2.00),
            ),
            name @ Climate::MaritimeTemperateOverLand => Self::new(
                name,
                (-0.39, 2.86),
                (141.7e3, 315.9e3, 167.4e3),
                (6.86, 10.38),
                (187.8e3, 169.6e3, 108.9e3),
                (8.58, 13.97),
                (216.0e3, 152.0e3, 122.7e3),
                1.518,
                1.282,
                (1.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ),
            name @ Climate::MaritimeTemperateOverSea => Self::new(
                name,
                (3.15, 857.9),
                (2222.0e3, 164.8e3, 116.3e3),
                (8.51, 169.8),
                (609.8e3, 119.9e3, 106.6e3),
                (8.43, 8.19),
                (136.2e3, 188.5e3, 122.9e3),
                1.518,
                1.282,
                (1.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ),
        }
    }
}
