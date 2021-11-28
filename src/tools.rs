use num_complex::Complex;

pub fn create_range(start: f64, end: f64, samples: usize) -> Vec<(usize, f64)> {
    let mut vec: Vec<(usize, f64)> = Vec::with_capacity(samples);

    let diff = match end - start {
        n if n > 0.0 => n,
        _ => panic!("Start greater than end"),
    };

    let mut increment: f64 = start;
    for counter in 0..samples {
        vec.push((counter, increment));
        increment += diff / ((samples - 1) as f64);
    }

    vec
}

pub fn create_colortable(colors: (f64, f64, f64), samples: usize) -> Vec<(u8, u8, u8)> {
    let pi = std::f64::consts::PI;
    // let rgb = (0.11, 0.02, 0.92);
    let mut vec: Vec<(u8, u8, u8)> = Vec::new();
    for i in 0..(samples + 1) {
        let value = 2.0 * pi * (i as f64) / samples as f64;
        vec.push((
            ((0.5 + 0.5 * (value * colors.0).sin()) * 255.0) as u8,
            ((0.5 + 0.5 * (value * colors.1).sin()) * 255.0) as u8,
            ((0.5 + 0.5 * (value * colors.2).sin()) * 255.0) as u8,
        ));
    }

    vec
}

pub struct MandelbrotPixel {
    pub x: u32,
    pub y: u32,
    pub xp: f64,
    pub yp: f64,
    pub x2: f64,
    pub y2: f64,
    pub n: u32,
    pub max_iterations: u32,
    pub orbitCount: f64,
    pub lambert_shading: f64,
    pub stripes: f64,
    pub surface_normal: (f64, f64),
    pub distance: f64,
    pub stripe_sig: f64,
    pub stripe_a: f64,
    pub diag: f64,
    pub anglex: f64,
    pub angley: f64,
    pub iteration: f64,
}

impl MandelbrotPixel {
    pub fn normalized_color(&self) -> u8 {
        // if self.n > 0 {
        //     println!("{}", self.n);
        // }
        (255.0 * ((self.n as f64) / (self.max_iterations as f64))) as u8
    }
    pub fn stripe_color(&self) -> u8 {
        // if self.n > 0 {
        //     println!("{}", self.n);
        // }
        // Calculate accorfing to asim
        (self.lambert_shading * 255.0 * self.stripe_a) as u8
    }
    // TODO: Write in python so its one to one
    pub fn smooth_iter_count_complex(&mut self) {
        let mut stripe_t: f64 = 0.0;
        let R: f64 = 10f64.powf(10.0);
        let mut nn = 0;

        let log2: f64 = 2.0f64.ln();
        let c = Complex::new(self.xp, self.yp);
        let mut z = Complex::new(0., 0.);
        let mut dz = Complex::new(1.0, 0.);
        for n in 0..self.max_iterations {
            nn = n;
            dz = dz * 2.0 * z + 1.;
            z = z * z + c;

            stripe_t = ((self.stripes * z.im.atan2(z.re)).sin() + 1.0) / 2.0;
            if z.re * z.re + z.im * z.im > (R) {
                break;
            }
            // stripe_t = (math.sin(stripe_s*math.atan2(pzi, pzr)) + 1) / 2

            self.stripe_a = self.stripe_a * self.stripe_sig + stripe_t * (1.0 - self.stripe_sig);
        }
        if nn < self.max_iterations - 1 {
            let modz = z.norm();
            let modz_ln = modz.ln();
            // let dii = dx * dx + dy * dy;

            let log_ratio = 2.0 * modz_ln / (R).ln();

            let mut smooth_i = 1.0 - log_ratio.ln() / log2;
            // smooth_i = 1.0 + ((R).ln() / zii.ln()).log2();

            self.stripe_a = self.stripe_a * (1.0 + smooth_i * (self.stripe_sig - 1.0))
                + stripe_t * smooth_i * (1.0 - self.stripe_sig);
            self.stripe_a = self.stripe_a
                / (1.0
                    - self.stripe_sig.powf(nn as f64) * (1.0 + smooth_i * (self.stripe_sig - 1.0)));

            if self.stripe_a > 1.0 {
                self.stripe_a = 1.0
            }
            if self.stripe_a < 0.0 {
                self.stripe_a = 0.0
            }
            self.iteration = smooth_i + nn as f64;
        } else {
            self.lambert_shading = 0.0;
            self.distance = 0.0;
            self.stripe_a = 0.0;
        }
    }

    pub fn smooth_iter_count(&mut self) {
        let mut zrr: f64 = 0.0;
        let mut zii: f64 = 0.0;
        let mut zr: f64 = self.xp; // x point
        let mut zi: f64 = self.yp; // y point
        let mut dx: f64 = 1.0;
        let mut dy: f64 = 0.0;
        let mut pzr: f64 = zr;
        let mut pzi: f64 = zi;
        let mut stripe_t: f64 = 0.0;
        let R: f64 = 10f64.powf(10.0);
        let h = 1.5;
        let mut nn = 0;
        let mut twori: f64 = 0.0;
        let mut new_dx: f64 = 0.0;
        let mut new_dy: f64 = 0.0;
        let Rii = R * R;
        let log2: f64 = 2.0f64.ln();

        for n in 0..self.max_iterations {
            nn = n;
            new_dx = 2.0 * (dx * pzr - dy * pzi + 0.5);
            new_dy = 2.0 * (dx * pzi + dy * pzr);
            zrr = zr * zr;
            zii = zi * zi;
            twori = 2.0 * zr * zi;
            zr = zrr - zii + self.xp;
            zi = twori + self.yp;
            if zrr + zii > (R) {
                break;
            }
            // stripe_t = (math.sin(stripe_s*math.atan2(pzi, pzr)) + 1) / 2
            stripe_t = ((self.stripes * zi.atan2(zr)).sin() + 1.0) / 2.0;
            dx = new_dx;
            dy = new_dy;
            pzr = zr;
            pzi = zi;

            self.stripe_a = self.stripe_a * self.stripe_sig + stripe_t * (1.0 - self.stripe_sig);
        }
        if nn < self.max_iterations - 1 {
            let zii = zi * zi + zr * zr;
            let modz = zii.sqrt();
            let modz_ln = modz.ln();
            let dii = dx * dx + dy * dy;

            let log_ratio = 2.0 * modz_ln / (R).ln();

            let smooth_i = 1.0 - log_ratio.ln() / log2;
            // smooth_i = 1.0 + ((R).ln() / zii.ln()).log2();
            self.stripe_a = self.stripe_a * (1.0 + smooth_i * (self.stripe_sig - 1.0))
                + stripe_t * smooth_i * (1.0 - self.stripe_sig);
            self.stripe_a = self.stripe_a
                / (1.0
                    - self.stripe_sig.powf(nn as f64) * (1.0 + smooth_i * (self.stripe_sig - 1.0)));
            self.distance = (modz * modz_ln / (dii).sqrt()) / 2.0;
            self.distance = self.distance / self.diag; //diag
            self.distance = -self.distance.ln() / 12.0;
            self.distance = 1.0 / (1.0 + (-10.0 * ((2.0 * self.distance - 1.0) / 2.0)).exp());

            let mut ux = (zr * dx + (zi * dy)) / (dx * dx + dy * dy);
            let mut uy = (-zr * dy + (zi * dx)) / (dx * dx + dy * dy);

            let uxx = ux / (ux * ux + uy * uy).sqrt();
            uy = uy / (ux * ux + uy * uy).sqrt();
            ux = uxx;
            let mut t = ux * self.anglex + uy * self.angley + h;
            t = t / (1.0 + h);
            if t > 1.0 {
                t = 1.0;
            }
            if t < 0.0 {
                t = 0.0;
            }
            self.lambert_shading = t;
            self.surface_normal = (ux, uy);
            self.iteration = smooth_i + nn as f64;
        } else {
            self.lambert_shading = 0.0;
            self.distance = 0.0;
            self.stripe_a = 0.0;
        }
    }

    pub fn calculate_according_to_asim(&mut self) -> () {
        // WORKING
        // let mut orbitCount += 0.5+0.5*sin(stripes*atan2(zi, zr));
        // let stripes = 15.; Nice
        // let stripes = 3.;

        // let mut zr: f64 = self.xp;
        // let mut zi: f64 = self.yp;
        let mut orbitCount = 0.0;
        // let R = 10f64.powf(10.0);
        // let mut zrr = 0.0;
        // let mut zii = 0.0;
        // let mut twori = 0.0;
        // let mut dx = 1.0;
        // let mut dy = 0.0;
        let h = 1.5;
        // let mut new_dx = 0.0;
        // let mut new_dy = 0.0;
        let mut zrr: f64 = 0.0;
        let mut zii: f64 = 0.0;
        let mut zr: f64 = self.xp; // x point
        let mut zi: f64 = self.yp; // y point
        let mut dx: f64 = 1.0;
        let mut dy: f64 = 0.0;
        let mut pzr: f64 = zr;
        let mut pzi: f64 = zi;
        let mut stripe_t: f64 = 0.0;
        let R: f64 = 10f64.powf(10.0);
        let h = 1.5;
        let mut nn = 0;
        let mut twori: f64 = 0.0;
        let mut new_dx: f64 = 0.0;
        let mut new_dy: f64 = 0.0;
        let Rii = R * R;
        let log2: f64 = 2.0f64.ln();

        while self.n < self.max_iterations {
            dx += 2.0 * zr + 1.0;
            dy += 2.0 * zi;
            zrr = zr * zr;
            zii = zi * zi;
            twori = 2.0 * zr * zi;
            zr = zrr - zii + self.xp;
            zi = twori + self.yp;
            new_dx = 2.0 * (dx * pzr - dy * pzi + 0.5);
            new_dy = 2.0 * (dx * pzi + dy * pzr);
            dx = new_dx;
            dy = new_dy;
            // pzr = zr;
            if zrr + zii > R {
                break;
            };
            // pzi = zi;

            // orbitCount += 0.5 + 0.5 * (stripes * zi.atan2(zr)).sin();
            orbitCount += ((self.stripes * zi.atan2(zr)).sin() + 1.0) / 2.0;
            // orbitCount = orbitCount * self.stripe_sig + stripe_t * (1.0 - self.stripe_sig);
            pzr = zr;
            pzi = zi;
            self.n += 1;
        }
        if self.n == self.max_iterations {
            self.lambert_shading = 0.0;
            self.stripe_a = 0.0;
        } else {
            // println!("{} {}", dx, dy);
            let lastOrbit = 0.5 + 0.5 * (self.stripes * pzi.atan2(pzr)).sin();
            let mut smallCount = orbitCount - lastOrbit;
            orbitCount = orbitCount / (self.n as f64);
            smallCount = smallCount / ((self.n as f64) - 1.0);
            // let frac = -1.0 + (2.0 * R.ln()).ln() / 2f64.ln()
            //     - (0.5 * (pzr * pzr + pzi * pzi).ln()).ln()
            //         / 2f64.ln();
            let sq = (pzr * pzr + pzi * pzi);
            let frac = 1.0 + ((R).ln() / sq.ln()).log2();
            let mix = frac * orbitCount + (1.0 - frac) * smallCount;

            let modz = (zi * zi + zr * zr).sqrt();
            //  modulus(z)*log(modulus(z))/modulus(dz)
            self.distance = modz * modz.ln() / (dx * dx + dy * dy).sqrt() / 2.0;
            self.distance = self.distance / self.diag; //diag
            self.distance = 1.0 / (1.0 + (-10.0 * ((2.0 * self.distance - 1.0) / 2.0)).exp());
            self.distance = -self.distance.ln() / 12.0;
            let mut ux = (zr * dx + (zi * dy)) / (dx * dx + dy * dy);
            let mut uy = (-zr * dy + (zi * dx)) / (dx * dx + dy * dy);

            let uxx = ux / (ux * ux + uy * uy).sqrt();
            uy = uy / (ux * ux + uy * uy).sqrt();
            ux = uxx;
            let mut t = ux * self.anglex + uy * self.angley + h;
            t = t / (1.0 + h);
            if t > 1.0 {
                t = 1.0;
            }
            if t < 0.0 {
                t = 0.0;
            }

            self.lambert_shading = t;
            self.surface_normal = (ux, uy);
            self.stripe_a = mix;
            self.iteration = (self.n as f64) + frac;
        }
    }
}
