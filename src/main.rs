#![forbid(unsafe_code)]
#![cfg_attr(feature = "cargo-clippy", deny(clippy_pedantic))]

extern crate num_complex;
extern crate gdal;
extern crate geo;

use geo::Point;
use num_complex::Complex64;
use std::path::PathBuf;

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

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Climate {
    Equatorial, // 1
    ContinentalSubtropical, // 2
    MaritimeSubtropical, // 3
    Desert, // 4
    ContinentalTemperate, // 5
    MaritimeTemperateOverLand, // 6
    MaritimeTemperateOverSea, // 7
}

impl Default for Climate {
    fn default() -> Self {
        Climate::ContinentalTemperate
    }
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Polarisation {
    Horizontal, // 0
    Vertical, // 1
    Dual, // n/a
}

impl Default for Polarisation {
    fn default() -> Self {
        Polarisation::Horizontal
    }
}

//
// void point_to_point_ITM(double elev[], double tht_m, double rht_m, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double frq_mhz, int radio_climate, int pol, double conf, double rel, double &dbloss, char *strmode, int &errnum)
//
// /******************************************************************************
//
// Note that point_to_point has become point_to_point_ITM for use as the old ITM
//
// 	pol:
// 		0-Horizontal, 1-Vertical
//
// 	radio_climate:
// 		1-Equatorial, 2-Continental Subtropical,
// 		3-Maritime Tropical, 4-Desert, 5-Continental Temperate,
// 		6-Maritime Temperate, Over Land, 7-Maritime Temperate,
// 		Over Sea
//
// 	conf, rel: .01 to .99
//
// 	elev[]: [num points - 1], [delta dist(meters)],
// 	        [height(meters) point 1], ..., [height(meters) point n]
//
// 	errnum: 0- No Error.
// 		1- Warning: Some parameters are nearly out of range.
// 		            Results should be used with caution.
// 		2- Note: Default parameters have been substituted for
// 		         impossible ones.
// 		3- Warning: A combination of parameters is out of range.
// 			    Results are probably invalid.
// 		Other-  Warning: Some parameters are out of range.
// 			Results are probably invalid.
//
// *****************************************************************************/
// {
//  // var decls
// 	prop_type   prop;
// 	propv_type  propv;
// 	propa_type  propa;
// 	double zsys=0;
// 	double zc, zr;
// 	double eno, enso, q;
// 	long ja, jb, i, np;
// 	double fs;
//
//  // prop things
// 	prop.hg[0]=tht_m;
// 	prop.hg[1]=rht_m;
// 	propv.klim=radio_climate;
// 	prop.kwx=0;
// 	propv.lvar=5;
// 	prop.mdp=-1;
//
// 	// qerfis
// 	zc=qerfi(conf);
// 	zr=qerfi(rel);
//
// 	np=(long)elev[0];
// 	eno=eno_ns_surfref;
// 	enso=0.0;
// 	q=enso;
//
//  // q
// 	if (q<=0.0)
// 	{
// 		ja=(long)(3.0+0.1*elev[0]);  /* added (long) to correct */
// 		jb=np-ja+6;
//
// 		for (i=ja-1; i<jb; ++i)
// 			zsys+=elev[i];
//
// 		zsys/=(jb-ja+1);
// 		q=eno;
// 	}
//
// 	propv.mdvar=12;
//
// 	qlrps(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,prop);
// 	qlrpfl(elev,propv.klim,propv.mdvar,prop,propa,propv);
//
// 	// final
// 	fs=32.45+20.0*log10(frq_mhz)+20.0*log10(prop.dist/1000.0);
// 	dbloss=avar(zr,0.0,zc,prop,propv)+fs;
// 	errnum=prop.kwx;
// }

/// Computes the propagation at a point along a point-to-point line.
///
/// This function is black magic to me. I initially ported it straight from the
/// SPLAT source, then adjusted it to fit a Rusty aesthetic. The original source
/// is unscrutable. Even deriving the meaning of variables and subroutine names
/// was an Herculean effort.
///
/// The original returned, in addition to a single number, a modestring
/// containing what I guessed was diagnostic or analysis results for two or
/// three possibilities along two axes, and a warning code which indicated
/// when results would be non-sensical because parameters were out of range.
///
/// I have completely discarded the modestring, and simplified the code through
/// that. I have also taken a harsher route on the out of range thing: if params
/// are out of range, or results are non-sensical, an error is returned. That
/// might make Splash somewhat less useful for extreme cases: fine by me.
///
/// ## Errors
///
/// Note that this function errors early, so one error may shadow another.
///
///  - If any parameter is out of range, an `OutOfRange(ParamName)` is
///    returned. `ParamName` is an enum describing which parameter was found to
///    be out of range.
///
///  - If any parameter is **impossible**, i.e. not out of range but, in
///    combination with the rest of the parameters, of a non-sensical value,
///    an `Impossible(ParamName)` is returned.
pub fn propagation(
    distance: f64,
    elevations: &Vec<f64>,
    tx_height: f64,
    rx_height: f64,
    dielectric: f64,
    conductivity: f64,
    surfref: f64,
    freq: f64,
    climate: Climate,
    polarisation: Polarisation,
    conf: f64,
    rel: f64
) -> Result<f64, PropagationError> {
    // var decls
    let mut zsys = 0.0;

    // prop things
    let mut propa = PropA::default();

    let mut prop = Prop::default();
    prop.hg = (tx_height, rx_height);
    prop.mdp = -1;
	prop.ens = surfref;
	prop.wn = freq / 47.7;

    let mut propv = PropV::default();
    propv.klim = climate;
    propv.lvar = 5;
    propv.mdvar = 12;

    // qerfis
    let zc = qerfi(conf);
    let zr = qerfi(rel);

    // As far as I can tell, this takes the sum of elevations
    // from the second element to the penultimate one, and then
    // takes the average, but it does so through fucked up math.
    //
    if false { // original
        // Gotta keep in mind that this math is for the original
        // elevations array, which was:
        //
        // [ num points, distance, ...elevations 1 through N... ]
        //
        // and `np = elev[0] - 1`. That explains some of the fuckedup,
        // but not why it was trying to do array indexing with floats,
        // nor the arbitrary numbers all over the place.
        let np = elevations.len() as f64 - 1.0;

        let ja = 1.0 + 0.1 * np;
        let jb = np - ja + 6.0;

        for i in ((ja - 1.0) as usize)..(jb as usize) {
            zsys += elevations[i];
        }

        zsys /= jb - ja + 1.0;
    } else { // rusty
        let subset = &elevations[1..(elevations.len() - 1)];
        let sum: f64 = subset.iter().sum();
        zsys = sum / (subset.len() as f64);
    }

    qlrps(zsys, polarisation, dielectric, conductivity, &mut prop);
    qlrpfl(distance, elevations, climate, propv.mdvar, &mut prop, &mut propa, &mut propv);

    let fs = 32.45 + 20.0 * (freq.log10() + (prop.dist / 1000.0).log10());
    println!(
        "avar(zr: {:?}, 0.0, zc: {:?}, prop: {:?}, propv: {:?}) + fs: {:?}",
        zr, zc, prop, propv, fs
    );
    // Ok(avar(zr, 0.0, zc, &mut prop, &mut propv) + fs);

    Ok(0.0)
}

/// Some kind of normalising function? Who knows though.
fn qerfi(q: f64) -> f64 {
	let c0 = 2.515516698;
	let c1 = 0.802853;
	let c2 = 0.010328;
	let d1 = 1.432788;
	let d2 = 0.189269;
	let d3 = 0.001308;

	let x = 0.5 - q;
	let mut t=(0.5 - x.abs()).max(0.000001);
	t = (-2.0 * t.ln()).sqrt();
	let v = t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0);

	if x < 0.0 { -v } else { v }
}

fn qlrps(
    zsys: f64,
    pol: Polarisation,
    dielect: f64,
    conduct: f64,
    prop: &mut Prop
) {
	let gma = 157e-9;

	if zsys != 0.0 {
		prop.ens *= (-zsys/9460.0).exp();
    }

	prop.gme = gma * (1.0 - 0.04665 * (prop.ens/179.3).exp());

    let zq = Complex64::new(dielect, 376.62 * conduct / prop.wn);
    let mut prop_zgnd = (zq - 1.0).sqrt();

	if pol == Polarisation::Vertical {
		prop_zgnd = prop_zgnd / zq;
    }

	prop.zgndreal = prop_zgnd.re;
	prop.zgndimag = prop_zgnd.im;
}

fn qlrpfl(
    distance: f64,
    elevations: &Vec<f64>,
    klimx: Climate,
    mdvarx: isize,
    mut prop: &mut Prop,
    mut propa: &mut PropA,
    mut propv: &mut PropV
) {
	prop.dist = elevations.len() as f64 * distance;
	let np = elevations.len();
	hzns(distance, &elevations, &mut prop);

    fn make_xl(hg: f64, dl: f64) -> f64 {
        (15.0 * hg).min(0.1 * dl)
    }

    let mut q;
    let z;
    let mut xl = (
        make_xl(prop.hg.0, prop.dl.0),
        make_xl(prop.hg.1, prop.dl.1)
    );

	xl.1 = prop.dist - xl.1;
	prop.dh = d1thx(distance, elevations, xl);

	if prop.dl.0 + prop.dl.1 > 1.5 * prop.dist {
        let (_, nz) = z1sq1(distance, elevations, xl);
        z = nz;
		prop.he = (
            prop.hg.0 + fortran_dim(elevations[0], z.0),
		    prop.hg.1 + fortran_dim(elevations[np], z.1)
        );

        fn make_dl(he: f64, gme: f64, dh: f64) -> f64 {
	        (2.0 * he / gme).sqrt() * (-0.07 * (dh / he.max(5.0)).sqrt()).exp()
        }

        prop.dl = (
            make_dl(prop.he.0, prop.gme, prop.dh),
            make_dl(prop.he.1, prop.gme, prop.dh)
        );

		q = prop.dl.0 + prop.dl.1;

		if q <= prop.dist { /* if there is a rounded horizon, or two obstructions, in the path */
			let temp = prop.dist / q;
			q = temp * temp;

            fn make_hedl(q: f64, he: f64, gme: f64, dh: f64) -> (f64, f64) {
                let he = he * q; /* tx effective height set to be path dist/distance between obstacles */
                (he, (2.0 * he / gme).sqrt() * (-0.07 * (dh / he.max(5.0)).sqrt()).exp())
            }

            let hedl = (
                make_hedl(q, prop.he.0, prop.gme, prop.dh),
                make_hedl(q, prop.he.1, prop.gme, prop.dh)
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
            make_qthe(prop.he.1, prop.gme, prop.dh, prop.dl.1)
        );
	} else {
        let (_, (z0, _)) = z1sq1(distance, elevations, (xl.0, 0.9 * prop.dl.0));
        let (_, (_, z1)) = z1sq1(distance, elevations, (prop.dist - 0.9 * prop.dl.1, xl.1));

		prop.he = (
            prop.hg.0 + fortran_dim(elevations[0], z0),
            prop.hg.1 + fortran_dim(elevations[np - 1], z1),
        );
	}

	prop.mdp = -1;
	propv.lvar = propv.lvar.max(3);

	if mdvarx >= 0 {
		propv.mdvar = mdvarx;
		propv.lvar = propv.lvar.max(4);
	}

    propv.klim = klimx;
    propv.lvar = 5;

    println!("lrprop(0.0, prop: {:?}, propa: {:?})", prop, propa);
	//lrprop(0.0, prop, propa);
}

fn hzns(distance: f64, elevations: &Vec<f64>, prop: &mut Prop) {
	let np = elevations.len();
    let xi = distance;

	let za = elevations[0] + prop.hg.0;
	let zb = elevations[np-1] + prop.hg.1;

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
					prop.the.1 += q/sb;
					prop.dl.1 =sb;
				}
			}
		}
	}
}

fn d1thx(distance: f64, elevations: &Vec<f64>, xl: (f64, f64)) -> f64 {
	let np = elevations.len();
	let mut xa = xl.0 / distance;
	let xb = xl.1 / distance;

	if (xb - xa) < 2.0 { return 0.0; }

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

fn z1sq1 (
    distance: f64,
    elevations: &Vec<f64>,
    x: (f64, f64)
) -> ((f64, f64), (f64, f64)) {
	let xn = elevations.len() as f64;
	let mut xa = fortran_dim(x.0 / distance, 0.0) as isize as f64;
	let mut xb = xn - (fortran_dim(xn, x.1 / distance) as isize as f64);

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

	let mut a = 0.5 * (elevations[ja] + elevations[jb]);
	let mut b = 0.5 * (elevations[ja] - elevations[jb]) * xx;

	for _ in 2..=n {
		ja += 1;
		xx += 1.0;
		a += elevations[ja];
		b += elevations[ja] * xx;
	}

	a /= xa;
	b = b * 12.0 / ((xa * xa + 2.0) * xa);

    (x, (a - b * xb, a + b * (xn - xb)))
}

fn fortran_dim(x: f64, y: f64) -> f64 {
    (x - y).max(0.0)
}

fn qtile (nn: usize, elevations: &mut Vec<f64>, ir: usize) -> f64 {
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

		while i <= n && elevations[i] >= q { i += 1; }

		if i > n { i = n; }

		j = j1;

		while j >= m && elevations[j] <= q { j -= 1; }

		if j < m { j = m; }

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

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Prop {
    pub aref: f64,
	pub dist: f64,

    /// Heights (tx, rx)
	pub hg: (f64, f64),
	pub rch: (f64, f64),
	pub wn: f64,
	pub dh: f64,
	pub dhd: f64,
	pub ens: f64,
	pub encc: f64,
	pub cch: f64,
	pub cd: f64,
	pub gme: f64,
	pub zgndreal: f64,
	pub zgndimag: f64,
	pub he: (f64, f64),
	pub dl: (f64, f64),
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
	pub kwx: isize,
	pub mdp: isize,
	pub ptx: isize,
	pub los: isize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct PropV {
    pub sgc: f64,
	pub lvar: isize,
	pub mdvar: isize,
	pub klim: Climate,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct PropA {
    pub dlsa: f64,
    pub dx:  f64,
    pub ael: f64,
    pub ak1: f64,
    pub ak2: f64,
    pub aed: f64,
    pub emd: f64,
    pub aes: f64,
    pub ems: f64,
    pub dls: (f64, f64),
    pub dla: f64,
    pub tha: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum PropagationError {
    OutOfBounds(ParamName),
    Impossible(ParamName)
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum ParamName {
    Elevation,
    Dielectric,
    Conductivity,
    SurfRef,
    Frequency,
    Climate,
    Polarisation,
    Conf,
    Rel,
}

fn main() {
    println!("Hello, world!");
}
