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

pub struct MandelbrotPixel {
    pub x: u32,
    pub y: u32,
    pub xp: f64,
    pub yp: f64,
    pub x2: f64,
    pub y2: f64,
    pub n: u32,
    pub max_iterations: u32,
    pub orbit_count: f64,
    pub lambert_shading: f64,
    pub stripes: f64,
    pub surface_normal: (f64, f64),
    pub distance: f64,
    pub stripe_sig: f64,
    pub stripe_avg: f64,
    pub diag: f64,
    pub anglex: f64,
    pub angley: f64,
    pub iteration: f64,
}

impl MandelbrotPixel {
    pub fn normalized_color(&self) -> u8 {
        (255.0 * ((self.n as f64) / (self.max_iterations as f64))) as u8
    }
    pub fn stripe_color(&self) -> u8 {
        (self.lambert_shading * 255.0 * self.stripe_avg) as u8
    }
    pub fn smooth_iter_count_complex(&mut self) {
        // CREDIT https://github.com/jlesuffleur/gpu_mandelbrot
        let mut stripe_t: f64 = 0.0;
        let r: f64 = 10f64.powf(10.0);
        let mut nn = 0;

        let log2: f64 = 2.0f64.ln();
        let c = Complex::new(self.xp, self.yp);
        let mut z = Complex::new(0., 0.);
        let mut dz = Complex::new(1.0, 0.);
        let h = 1.5;
        for n in 0..self.max_iterations {
            nn = n;
            dz = dz * 2.0 * z + 1.;
            z = z * z + c;

            stripe_t = ((self.stripes * z.im.atan2(z.re)).sin() + 1.0) / 2.0;
            if z.re * z.re + z.im * z.im > (r) {
                break;
            }

            self.stripe_avg =
                self.stripe_avg * self.stripe_sig + stripe_t * (1.0 - self.stripe_sig);
        }
        if nn < self.max_iterations - 1 {
            let modz = z.norm();
            let modz_ln = modz.ln();

            let log_ratio = 2.0 * modz_ln / (r).ln();

            let smooth_i = 1.0 - log_ratio.ln() / log2;

            self.stripe_avg = self.stripe_avg * (1.0 + smooth_i * (self.stripe_sig - 1.0))
                + stripe_t * smooth_i * (1.0 - self.stripe_sig);
            self.stripe_avg = self.stripe_avg
                / (1.0
                    - self.stripe_sig.powf(nn as f64) * (1.0 + smooth_i * (self.stripe_sig - 1.0)));

            if self.stripe_avg > 1.0 {
                self.stripe_avg = 1.0;
            }
            if self.stripe_avg < 0.0 {
                self.stripe_avg = 0.0;
            }

            self.distance = (modz * modz_ln) / dz.norm() / 2.0;
            self.distance = self.distance / self.diag; //diag
            self.distance = -self.distance.ln() / 12.0;
            self.distance = 1.0 / (1.0 + (-10.0 * ((2.0 * self.distance - 1.0) / 2.0)).exp());

            let mut ux = (z.re * dz.re + (z.im * dz.im)) / (dz.re * dz.re + dz.im * dz.im);
            let mut uy = (-z.re * dz.im + (z.im * dz.re)) / (dz.re * dz.re + dz.im * dz.im);

            let uxx = ux / (ux * ux + uy * uy).sqrt();
            uy = uy / (ux * ux + uy * uy).sqrt();
            ux = uxx;
            let mut t = ux * self.anglex + uy * self.angley + h;
            t = t / (1.0 + h);

            self.lambert_shading = t;
            self.surface_normal = (ux, uy);
            self.iteration = smooth_i + nn as f64;
        } else {
            self.lambert_shading = 0.0;
            self.distance = 0.0;
            self.stripe_avg = 0.0;
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
        let r: f64 = 10f64.powf(10.0);
        let h = 1.5;
        let mut nn = 0;
        let mut twori: f64;
        let mut new_dx: f64;
        let mut new_dy: f64;
        let log2: f64 = 2.0f64.ln();

        for n in 0..self.max_iterations {
            nn = n;
            new_dx = 2.0 * (dx * pzr - dy * pzi + 0.5);
            new_dy = 2.0 * (dx * pzi + dy * pzr);
            if zrr + zii > (r) {
                break;
            }
            zrr = zr * zr;
            zii = zi * zi;
            twori = 2.0 * zr * zi;
            zr = zrr - zii + self.xp;
            zi = twori + self.yp;
            stripe_t = ((self.stripes * zi.atan2(zr)).sin() + 1.0) / 2.0;
            dx = new_dx;
            dy = new_dy;
            pzr = zr;
            pzi = zi;

            self.stripe_avg =
                self.stripe_avg * self.stripe_sig + stripe_t * (1.0 - self.stripe_sig);
        }
        if nn < self.max_iterations - 1 {
            let zii = zi * zi + zr * zr;
            let modz = zii.sqrt();
            let modz_ln = modz.ln();
            let dii = dx * dx + dy * dy;

            let log_ratio = 2.0 * modz_ln / (r).ln();

            let smooth_i = 1.0 - log_ratio.ln() / log2;
            self.stripe_avg = self.stripe_avg * (1.0 + smooth_i * (self.stripe_sig - 1.0))
                + stripe_t * smooth_i * (1.0 - self.stripe_sig);
            self.stripe_avg = self.stripe_avg
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

            self.lambert_shading = t;
            self.surface_normal = (ux, uy);
            self.iteration = smooth_i + nn as f64;
        } else {
            self.lambert_shading = 0.0;
            self.distance = 0.0;
            self.stripe_avg = 0.0;
        }
    }

    pub fn smooth_iteration_count_2(&mut self) -> () {
        // CREDIT user asimes: https://www.fractalforums.com/general-discussion/stripe-average-coloring/

        let mut orbit_count = 0.0;
        let mut zrr: f64;
        let mut zii: f64;
        let mut zr: f64 = self.xp; // x point
        let mut zi: f64 = self.yp; // y point
        let mut dx: f64 = 1.0;
        let mut dy: f64 = 0.0;
        let mut pzr: f64 = zr;
        let mut pzi: f64 = zi;
        let r: f64 = 10f64.powf(10.0);
        let h = 1.5;
        let mut twori: f64;
        let mut new_dx: f64;
        let mut new_dy: f64;

        let mut nn = 0;

        for n in 0..self.max_iterations {
            nn = n;
            new_dx = 2.0 * (dx * zr - dy * zi + 0.5);
            new_dy = 2.0 * (dx * zi + dy * zr);
            zrr = zr * zr;
            zii = zi * zi;
            twori = 2.0 * zr * zi;
            let _zr = zrr - zii + self.xp;
            let _zi = twori + self.yp;
            if zrr + zii > (r) {
                break;
            }
            zr = _zr;
            zi = _zi;
            orbit_count += ((self.stripes * zi.atan2(zr)).sin() + 1.0) / 2.0;
            dx = new_dx;
            dy = new_dy;
            pzr = zr;
            pzi = zi;
        }

        if nn < self.max_iterations - 1 {
            let last_orbit = 0.5 + 0.5 * (self.stripes * pzi.atan2(pzr)).sin();
            let mut small_count = orbit_count - last_orbit;
            orbit_count = orbit_count / (nn as f64);
            small_count = small_count / ((nn as f64) - 1.0);
            let sq = pzr * pzr + pzi * pzi;
            let frac = 1.0 + ((r).ln() / sq.ln()).log2();
            let mix = frac * orbit_count + (1.0 - frac) * small_count;

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
            // println!("{}", mix);
            self.lambert_shading = t;
            self.surface_normal = (ux, uy);
            self.stripe_avg = mix;
            self.iteration = (nn as f64) + frac;
        } else {
            self.lambert_shading = 0.0;
            self.stripe_avg = 0.0;
        }
    }
}
