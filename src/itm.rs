//! The Longley-Rice Irregular Terrain Model for RF.
//!
//! This module was initially derived from a C++ translation of the original
//! FORTRAN based of the ITM, the [Irregular Terrain Model][ITM68], also known
//! as Longley-Rice, an empirical RF propagation model developed at the National
//! Telecommunications and Information Administration of the United States in
//! the 1960s by Anita Longley and Phil Rice.
//!
//! While the initial reconstruction was done referring to the above-mentioned
//! source and its adaptation in [SPLAT!], all variables and function names were
//! still unscrutable, as was the general operation of the model. Thus, the code
//! was then cross-referenced back to the [LaTeX documentation][ITM122] of the
//! ITM version 1.2.2, and made more legible and understandable by deriving
//! meaningful names and incorporating documentation into this source, as well
//! as doing some restructuring by referring back to George Hufford's 1999 memo
//! describing “[The Algorithm][GH1999]” in exceedingly pleasant prose.
//!
//! The goal is to be able to read this documentation or source, and
//! _understand_ the algorithm, without referring back to other documents listed
//! above. Please file bugs if something is unclear.
//!
//! This implementation is released in the Public Domain, although note that the
//! NTIA requests any use of the ITM is properly credited.
//!
//! [GH1999]: https://www.its.bldrdoc.gov/media/50676/itm_alg.pdf
//! [ITM122]: https://www.its.bldrdoc.gov/media/50674/itm.pdf
//! [ITM68]: https://www.its.bldrdoc.gov/resources/radio-propagation-software/itm/itm.aspx
//! [SPLAT!]: http://www.qsl.net/kd2bd/splat.html

use num_complex::Complex64;

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

/// Point-to-Point propagation
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

	lrprop(0.0, prop, propa);
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

fn lrprop(d: f64, prop: &mut Prop, propa: &mut PropA) {
    // these 4 are cpp statics too
    let mut dmin: f64 = 0.0;
    let mut xae: f64 = 0.0;
    let mut wlos: bool = false;
	let mut wscat: bool = false;

    let prop_zgnd = Complex64::new(prop.zgndreal, prop.zgndimag);

	if prop.mdp != 0 {
        propa.dls.0 = (2.0 * prop.he.0 / prop.gme).sqrt();
        propa.dls.1 = (2.0 * prop.he.1 / prop.gme).sqrt();

		propa.dlsa = propa.dls.0 + propa.dls.1;
		propa.dla = prop.dl.0 + prop.dl.1;
		propa.tha = (prop.the.0 + prop.the.1).max(-propa.dla * prop.gme);

        wlos = false;
        wscat = false;


		if prop.wn < 0.838 || prop.wn > 210.0 {
			prop.kwx = prop.kwx.max(1);
        }

        fn make_kwx_hg(kwx: isize, hg: f64) -> isize {
            if hg < 1.0 || hg > 1000.0 {
                kwx.max(1)
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_hg(prop.kwx, prop.hg.0);
        prop.kwx = make_kwx_hg(prop.kwx, prop.hg.1);

        fn make_kwx_dl(kwx: isize, the: f64, dl: f64, dls: f64) -> isize {
            if (the.abs() > 200e-3) || (dl < 0.1 * dls) || (dl > 3.0 * dls) {
                kwx.max(3)
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_dl(prop.kwx, prop.the.0, prop.dl.0, propa.dls.0);
        prop.kwx = make_kwx_dl(prop.kwx, prop.the.1, prop.dl.1, propa.dls.1);

		if (prop.ens < 250.0) || (prop.ens > 400.0) || (prop.gme < 75e-9) || (prop.gme > 250e-9) || (prop_zgnd.re <= prop_zgnd.im.abs()) || (prop.wn < 0.419) || (prop.wn > 420.0) {
			prop.kwx = 4;
        }

        fn make_kwx_hg_again(kwx: isize, hg: f64) -> isize {
            if hg < 0.5 || hg > 3000.0 {
                4
            } else {
                kwx
            }
        }

        prop.kwx = make_kwx_hg_again(prop.kwx, prop.hg.0);
        prop.kwx = make_kwx_hg_again(prop.kwx, prop.hg.1);

		dmin = (prop.he.0 - prop.he.1).abs() / 200e-3;
		adiff(0.0, prop, propa);
		xae = (prop.wn * prop.gme.powi(2)).cbrt();
		let d3 = propa.dlsa.max(1.3787 * xae + propa.dla);
		let d4 = d3 + 2.7574 * xae;
		let a3 = adiff(d3, prop, propa);
		let a4 = adiff(d4, prop, propa);
		propa.emd = (a4 - a3) / (d4 - d3);
		propa.aed = a3 - propa.emd * d3;
	}

	if prop.mdp >= 0 {
		prop.mdp = 0;
		prop.dist = d;
	}

	if prop.dist > 0.0 {
		if prop.dist > 1000e3 {
			prop.kwx = prop.kwx.max(1);
        }

		if prop.dist < dmin {
			prop.kwx = prop.kwx.max(3);
        }

		if prop.dist < 1e3 || prop.dist > 2000e3 {
			prop.kwx = 4;
        }
	}

	if prop.dist < propa.dlsa {
		if !wlos {
			alos(0.0, prop, propa);
			let d2 = propa.dlsa;
			let a2 = propa.aed + d2 * propa.emd;
			let mut d0 = 1.908 * prop.wn * prop.he.0 * prop.he.1;

            let d1;
			if propa.aed >= 0.0 {
				d0 = d0.min(0.5 * propa.dla);
				d1 = d0 + 0.25 * (propa.dla - d0);
			} else {
				d1 = (-propa.aed / propa.emd).max(0.25 * propa.dla);
            }

			let a1 = alos(d1, prop, propa);

			if d0 < d1 {
				let a0 = alos(d0, prop, propa);
				let q = (d2 / d0).ln();
				propa.ak2 = 0.0f64.max(
                    ((d2 - d0) * (a1 - a0) - (d1 - d0) * (a2 - a0))
                    /
                    ((d2 - d0) * (d1 / d0).ln() - (d1 - d0) * q)
                );

                let wq = propa.aed >= 0.0 || propa.ak2 > 0.0;

				if wq {
					propa.ak1 = (a2 - a0 - propa.ak2 * q) / (d2 - d0);

					if propa.ak1 < 0.0 {
						propa.ak1 = 0.0;
						propa.ak2 = fortran_dim(a2, a0) / q;

						if propa.ak2 == 0.0 {
							propa.ak1 = propa.emd;
                        }
					}
				} else {
					propa.ak2 = 0.0;
					propa.ak1 = (a2 - a1) / (d2 - d1);

					if propa.ak1 <= 0.0 {
						propa.ak1 = propa.emd;
                    }
				}
			} else {
				propa.ak1 = (a2 - a1) / (d2 - d1);
				propa.ak2 = 0.0;

				if propa.ak1 <= 0.0 {
					propa.ak1 = propa.emd;
                }
			}

			propa.ael = a2 - propa.ak1 * d2 - propa.ak2 * d2.ln();
			wlos = true;
		}

		if  prop.dist > 0.0 {
			prop.aref = propa.ael + propa.ak1 * prop.dist + propa.ak2 * prop.dist.ln();
        }
	}

	if prop.dist <= 0.0 || prop.dist >= propa.dlsa {
		if !wscat {
			ascat(0.0, prop, propa);
			let d5 = propa.dla + 200e3;
			let d6 = d5 + 200e3;
			let a6 = ascat(d6, prop, propa);
			let a5 = ascat(d5, prop, propa);

			if a5 < 1000.0 {
				propa.ems = (a6 - a5) / 200e3;
				propa.dx = propa.dlsa.max(
                    (propa.dla + 0.3 * xae * (47.7 * prop.wn).ln()).max(
                        (a5 - propa.aed - propa.ems * d5) / (propa.emd - propa.ems)
                    )
                );
				propa.aes = (propa.emd - propa.ems) * propa.dx + propa.aed;
			} else {
				propa.ems = propa.emd;
				propa.aes = propa.aed;
				propa.dx = 10.0e6;
			}

			wscat = true;
		}

		if prop.dist > propa.dx {
			prop.aref = propa.aes + propa.ems * prop.dist;
        } else {
			prop.aref = propa.aed + propa.emd * prop.dist;
        }
	}

	prop.aref = prop.aref.max(0.0);
}

fn adiff(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
    let prop_zgnd = Complex64::new(prop.zgndreal, prop.zgndimag);

    // this is setting up data for iterations beyond the first (d=0 == first)
    // gotta figure out how to do this in rust...
    // ...without sharing data between *different* runs!
	if d == 0.0 {
		let mut q = prop.hg.1 * prop.hg.1;
		let qk = prop.he.0 * prop.he.1 - q;

		if prop.mdp < 0 {
			q += 10.0;
        }

		let wd1 = (1.0 + qk / q).sqrt();
		let xd1 = propa.dla + propa.tha / prop.gme;

		q = (1.0 - 0.8 * (-propa.dlsa / 50e3 ).exp()) * prop.dh;
		q *= 0.78 * (-(q / 16.0).powf(0.25)).exp();

		let afo = 15.0f64.min(2.171 * (1.0 + 4.77e-4 * prop.hg.0 * prop.hg.1 * prop.wn * q).ln());
		let qk = 1.0 / prop_zgnd.norm_sqr().sqrt();
		let mut aht = 20.0;
		let mut xht = 0.0;

        if false {
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

		0.0 // ???
	} else {
		let th = propa.tha + d * prop.gme;
		let ds = d - propa.dla;
		let mut q = 0.0795775 * prop.wn * ds * th * th;
		let adiffv = aknfe(q * prop.dl.0 / (ds + prop.dl.0))
            + aknfe(q * prop.dl.1 / (ds + prop.dl.1));

        // These were added in as they were not set before in this branch wtf???
        let qk = 0.0;
        let xht = 0.0;
        let aht = 0.0;
        let wd1 = 0.0;
        let xd1 = 0.0;
        let afo = 0.0;

        let a = ds / th;
		let wa = (a * prop.wn).cbrt();
		let pk = qk / wa;
		q = (1.607 - pk) * 151.0 * wa * th + xht;
		let ar = 0.05751 * q - 4.343 * q.ln() - aht;
		q = (wd1 + xd1 / d) * 6283.2f64.min(
            (1.0 - 0.8 * (-d / 50e3).exp()) * prop.dh * prop.wn
        );

        let wd = 25.1 / (25.1 + q.sqrt());
		ar * wd + (1.0 - wd) * adiffv + afo
	}
}

fn alos(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
	let prop_zgnd = Complex64::new(prop.zgndreal, prop.zgndimag);
    let mut wls = 0.0;

    // see adiff comment
	if d == 0.0 {
		wls = 0.021 / (0.021 + prop.wn * prop.dh / 10e3f64.max(propa.dlsa));
		0.0
	} else {
		let mut q = (1.0 - 0.8 * (-d / 50e3).exp()) * prop.dh;
		let s = 0.78 * q * (-(q / 16.0).powf(0.25)).exp();
		q = prop.he.0 + prop.he.1;
		let sps = q / (d.powi(2) + q.powi(2)).sqrt();
		let mut r = (sps - prop_zgnd) / (sps + prop_zgnd)
            * (-10.0f64.min(prop.wn * s * sps)).exp();
		q = r.norm_sqr();

		if q < 0.25 || q < sps {
			r = r * (sps / q).sqrt();
        }

		let alosv = propa.emd * d + propa.aed;
		q = prop.wn * prop.he.0 * prop.he.1 * 2.0 / d;

		if q > 1.57 {
			q = 3.14 - 2.4649 / q;
        }

        let qq = Complex64::new(q.cos(), -q.sin());
		(-4.343 * (qq + r).norm_sqr().ln() - alosv) * wls + alosv
	}
}

fn ascat(d: f64, prop: &mut Prop, propa: &mut PropA) -> f64 {
	// static double ad, rr, etq, h0s;
	// double h0, r1, r2, z0, ss, et, ett, th, q;
	// double ascatv, temp;

    let ad = 0.0;
    let rr = 0.0;
    let mut etq = 0.0;
    let mut h0s = 0.0;

    // see adiff comment
	if d == 0.0 {
		let mut ad = prop.dl.0 - prop.dl.1;
		let mut rr = prop.he.1 / prop.rch.0;

		if ad < 0.0 {
			ad = -ad;
			rr = 1.0/rr;
		}

		etq = (5.67e-6 * prop.ens - 2.32e-3) * prop.ens + 0.031;
		h0s =- 15.0;
		0.0
	} else {
        let mut h0;
		if h0s > 15.0 {
			h0 = h0s;
		} else {
			let th = prop.the.0 + prop.the.1 + d * prop.gme;
			let mut r2 = 2.0 * prop.wn * th;
			let r1 = r2 * prop.he.0;
			r2 *= prop.he.1;

			if r1 < 0.2 && r2 < 0.2 { return 1001.0; }

			let mut ss = (d - ad) / (d + ad);
			let mut q = rr / ss;
			ss = ss.max(0.1);
			q = q.max(0.1).min(10.0);
			let z0 = (d - ad) * (d + ad) * th * 0.25 / d;

			let temp = (z0/8.0e3).min(1.7).powi(6);
			let et = (etq * (-temp).exp() + 1.0) * z0 / 1.7556e3;

			let ett = et.max(1.0);
			h0 = (h0f(r1, ett) + h0f(r2, ett)) * 0.5;
			h0 += 1.38 - ett.ln().min(h0) * ss.ln() * q.ln() * 0.49;
			h0 = fortran_dim(h0, 0.0);

			if et < 1.0 {
				let temp = (1.0 + 1.4142 / r1) * (1.0 + 1.4142 / r2);
				h0 = et * h0 + (1.0 - et) * 4.343
                    * (temp.powi(2) * (r1 + r2) / (r1 + r2 + 2.8284)).ln();
			}

			if h0 > 15.0 && h0s >= 0.0 {
                h0 = h0s;
            }
		}

		h0s = h0;
		let th = propa.tha + d * prop.gme;
		ahd(th * d) + 4.343 * (47.7 * prop.wn * th.powi(4)).ln() - 0.1
            * (prop.ens - 301.0) * (-th * d / 40e3).exp() + h0
	}
}

fn aknfe(v2: f64) -> f64
{
	if v2 < 5.76 {
		6.02 + 9.11 * v2.sqrt() - 1.27 * v2
	} else {
		12.953 + 10.0 * v2.log10()
    }
}

fn fht(x: f64, pk: f64) -> f64 {
    let mut fhtv;

    if x < 200.0 {
		let w = -pk.ln();

		if pk < 1.0e-5 || x * w.powi(3) > 5495.0 {
			fhtv =- 117.0;

			if x > 1.0 {
				fhtv = 40.0 * x.log10() + fhtv;
            }
		} else {
			fhtv = 2.5e-5 * x.powi(2) / pk - 8.686 * w - 15.0;
        }
	} else {
		fhtv = 0.05751 * x - 10.0 * x.log10();

		if x < 2000.0 {
			let w = 0.0134 * x * (-0.005 * x).exp();
			fhtv = (1.0 - w) * fhtv + w * (40.0 * x.log10() - 117.0);
		}
	}

	fhtv
}

fn h0f(r: f64, et: f64) -> f64 {
	let a = [25.0, 80.0, 177.0, 395.0, 705.0];
	let b = [24.0, 45.0,  68.0,  80.0, 105.0];

	let (it, q) = if et <= 0.0 {
		(1, 0.0)
	} else if et >= 5.0 {
		(5, 0.0)
	} else {
        (et as usize, et.fract())
    };

    let x = (1.0 / r).powi(2);
	let h0fv = 4.343 * ((a[it - 1] * x + b[it - 1]) * x + 1.0).ln();

	if q == 0.0 {
        h0fv
    } else {
		(1.0 - q) * h0fv + q * 4.343 * ((a[it] * x + b[it]) * x + 1.0).ln()
    }
}

fn ahd(td: f64) -> f64 {
	let a = [   133.4,    104.6,     71.8];
	let b = [0.332e-3, 0.212e-3, 0.157e-3];
	let c = [  -4.343,   -1.086,    2.171];

	let i = if td <= 10e3 {
		0
    } else if td <= 70e3 {
		1
    } else {
		2
    };

    a[i] + b[i] * td + c[i] * td.ln()
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
