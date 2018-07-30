#![forbid(unsafe_code)]
#![cfg_attr(feature = "cargo-clippy", deny(clippy_pedantic))]

extern crate gdal;
extern crate geo;

use geo::Point;
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

/// Base site definition.
struct Site {
    /// Where it is
    position: Point<f64>,

    /// How high above ground
    aboveground: f64,

    /// What it's called
    name: String,
}

/// An RF transmitter and its parameters.
struct Transmitter {
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
fn propagation(
    distance: f64,
    elevations: Vec<f64>,
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
    /*
    let prop = Prop {
        hg: (tx_height, rx_height),
        kwx: 0,
        mdp: -1,
    };

    let propv = PropV {
        klim: climate,
        lvar: 5,
    };
    */

    // qerfis

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

        zsys /= (jb - ja + 1.0);
    } else { // rusty
        let subset = &elevations[1..(elevations.len() - 2)];
        let sum: f64 = subset.iter().sum();
        let zsys = sum / (subset.len() as f64);
    }

    // q = surfref

    // qlrps
    // qlrpfl

    // final

    Ok(0.0)
}

/// Some kind of normalising function? Who knows though.
fn qerfi(q: f64) -> f64 {
	let c0=2.515516698;
	let c1=0.802853;
	let c2=0.010328;
	let d1=1.432788;
	let d2=0.189269;
	let d3=0.001308;

	let x=0.5-q;
	let mut t=(0.5-x.abs()).max(0.000001);
	t=(-2.0*t.ln()).sqrt();
	let mut v=t-((c2*t+c1)*t+c0)/(((d3*t+d2)*t+d1)*t+1.0);

	if x < 0.0 { -v } else { v }
}

struct Prop {
    pub aref: f64,
	pub dist: f64,

    /// Heights (tx, rx)
	pub hg: (f64, f64),
	pub rch: [f64; 2],
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
	pub he: [f64; 2],
	pub dl: [f64; 2],
	pub the: [f64; 2],
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

struct PropV {
    pub sgc: f64,
	pub lvar: isize,
	pub mdvar: isize,
	pub klim: Climate,
}

enum PropagationError {
    OutOfBounds(ParamName),
    Impossible(ParamName)
}

enum ParamName {
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
