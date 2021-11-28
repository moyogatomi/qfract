use image;
use image::{Rgb, RgbImage};
use ndarray::{Array3, Array4, ArrayD, ArrayViewD, ArrayViewMutD};
use num_complex::Complex;
use numpy::{IntoPyArray, PyArray3, PyArray4, PyArrayDyn, PyReadonlyArrayDyn};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};
use rayon::prelude::*;
use std::time::Instant;

mod tools;
pub use crate::tools::{create_colortable, create_range, MandelbrotPixel};

pub fn calculate_mandelbrot(
    // center: (f64, f64),
    xstart: f64,
    xend: f64,
    ystart: f64,
    yend: f64,
    max_iterations: u32,
    xsize: f64,
    anglex: f64,
    angley: f64,
    stripes: f64,
    method: u16,
) -> Array3<f64> {
    println!("{},{},{},{}", xstart, xend, ystart, yend);

    let ysize = &xsize * { (yend - ystart) / (xend - xstart) };

    let x_range = create_range(xstart, xend, xsize as usize);
    let y_range = create_range(ystart, yend, ysize as usize);
    let diag = ((xend - xstart).powf(2.0) + (yend - ystart).powf(2.0)).sqrt();

    println!("image size: {}x{}", xsize, ysize as usize);
    println!("iterations: {}", max_iterations);

    let before = Instant::now();
    let colortable = create_colortable((0.11, 0.02, 0.92), 255);
    let mut a: Array3<f64> = Array3::zeros((xsize as usize, ysize as usize, 6));
    println!(
        "Elapsed time: {:.2?}, megapixels: {}",
        before.elapsed(),
        (&xsize * &ysize / 1000000.),
    );

    let results: Vec<Vec<MandelbrotPixel>> = x_range
        .par_iter()
        .map(|(xx, x)| {
            let new_iterator = y_range
                .iter()
                .map(|(yy, y)| {
                    let mut pixel = MandelbrotPixel {
                        x: *xx as u32,
                        y: *yy as u32,
                        xp: *x,
                        yp: *y,
                        x2: 0.,
                        y2: 0.,
                        n: 0,
                        max_iterations: max_iterations,
                        orbitCount: 0.0,
                        lambert_shading: 0.0,
                        stripes: stripes,
                        surface_normal: (0.0, 0.0),
                        distance: 0.0,
                        stripe_sig: 0.9,
                        stripe_a: 0.0,
                        diag: diag,
                        anglex: anglex,
                        angley: angley,
                        iteration: 0.0,
                    };
                    if method == 0 {
                        pixel.smooth_iter_count_complex();
                    } else if method == 1 {
                        pixel.smooth_iter_count();
                    } else if method == 2 {
                        pixel.calculate_according_to_asim();
                    }
                    pixel
                })
                .collect::<Vec<MandelbrotPixel>>();
            new_iterator
        })
        .collect::<Vec<Vec<MandelbrotPixel>>>();
    println!("Elapsed time: {:.2?}, after calc", before.elapsed(),);
    for i in results {
        for pixel in i {
            // let color = pixel.normalized_color();
            // let color = pixel.stripe_color();
            // let shade = pixel.shade;
            // let n = (255.0 * (pixel.n as f64 / pixel.max_iterations as f64));
            // // println!("{}", n);
            // let (colorx, colory, colorz) = colortable[n as usize];
            // let image_pixel = Rgb([
            //     (colorx as f64 * shade) as u8,
            //     (colory as f64 * shade) as u8,
            //     (colorz as f64 * shade) as u8,
            // ]);
            // let image_pixel = Rgb([color, color, color]);
            // img.put_pixel(pixel.x, pixel.y, image_pixel);
            a[[pixel.x as usize, pixel.y as usize, 0]] = pixel.lambert_shading;
            a[[pixel.x as usize, pixel.y as usize, 1]] = pixel.stripe_a;
            a[[pixel.x as usize, pixel.y as usize, 2]] = pixel.distance;
            a[[pixel.x as usize, pixel.y as usize, 3]] = pixel.surface_normal.0;
            a[[pixel.x as usize, pixel.y as usize, 4]] = pixel.surface_normal.1;
            a[[pixel.x as usize, pixel.y as usize, 5]] = pixel.iteration;
        }
    }
    println!(
        "Elapsed time: {:.2?}, megapixels: {}",
        before.elapsed(),
        (&xsize * &ysize / 1000000.),
    );
    // img.save("testx.png").expect("Not possible");
    // true
    a
}
#[pymodule]
fn rust_ext(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // immutable example
    fn axpy(a: f64, x: ArrayViewD<'_, f64>, y: ArrayViewD<'_, f64>) -> ArrayD<f64> {
        a * &x + &y
    }

    // mutable example (no return)
    fn mult(a: f64, mut x: ArrayViewMutD<'_, f64>) {
        x *= a;
    }

    // wrapper of `calculate_mandelbort`
    #[pyfn(m, "calculate")]
    fn calculate_py<'py>(
        py: Python<'py>,
        xstart: f64,
        xend: f64,
        ystart: f64,
        yend: f64,
        max_iterations: u32,
        xsize: f64,
        anglex: f64,
        angley: f64,
        stripes: f64,
        method: u16,
    ) -> &'py PyArray3<f64> {
        calculate_mandelbrot(
            xstart,
            xend,
            ystart,
            yend,
            max_iterations,
            xsize,
            anglex,
            angley,
            stripes,
            method,
        )
        .into_pyarray(py)
        // axpy(a, x, y).into_pyarray(py)
    }
    // wrapper of `axpy`
    #[pyfn(m, "axpy")]
    fn axpy_py<'py>(
        py: Python<'py>,
        a: f64,
        x: PyReadonlyArrayDyn<f64>,
        y: PyReadonlyArrayDyn<f64>,
    ) -> &'py PyArrayDyn<f64> {
        let x = x.as_array();
        let y = y.as_array();
        axpy(a, x, y).into_pyarray(py)
    }

    // wrapper of `mult`
    #[pyfn(m, "mult")]
    fn mult_py(_py: Python<'_>, a: f64, x: &PyArrayDyn<f64>) -> PyResult<()> {
        let x = unsafe { x.as_array_mut() };
        mult(a, x);
        Ok(())
    }

    Ok(())
}
