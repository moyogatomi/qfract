use ndarray::Array3;
use numpy::{IntoPyArray, PyArray3};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};
use rayon::prelude::*;
// use std::time::Instant;

mod tools;
pub use crate::tools::{create_range, MandelbrotPixel};

pub fn calculate_mandelbrot(
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
    // println!("{},{},{},{}", xstart, xend, ystart, yend);

    let ysize = &xsize * { (yend - ystart) / (xend - xstart) };

    let x_range = create_range(xstart, xend, xsize as usize);
    let y_range = create_range(ystart, yend, ysize as usize);
    let diag = ((xend - xstart).powf(2.0) + (yend - ystart).powf(2.0)).sqrt();

    // println!("image size: {}x{}", xsize, ysize as usize);
    // println!("iterations: {}", max_iterations);

    // let before = Instant::now();
    let mut a: Array3<f64> = Array3::zeros((xsize as usize, ysize as usize, 6));

    // println!(
    //     "Elapsed time: {:.2?}, megapixels: {}",
    //     before.elapsed(),
    //     (&xsize * &ysize / 1000000.),
    // );

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
                        orbit_count: 0.0,
                        lambert_shading: 0.0,
                        stripes: stripes,
                        surface_normal: (0.0, 0.0),
                        distance: 0.0,
                        stripe_sig: 0.9,
                        stripe_avg: 0.0,
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
                        pixel.smooth_iteration_count_2();
                    }
                    pixel
                })
                .collect::<Vec<MandelbrotPixel>>();
            new_iterator
        })
        .collect::<Vec<Vec<MandelbrotPixel>>>();
    // println!("Elapsed time: {:.2?}, after calc", before.elapsed(),);
    for i in results {
        for pixel in i {
            a[[pixel.x as usize, pixel.y as usize, 0]] = pixel.lambert_shading;
            a[[pixel.x as usize, pixel.y as usize, 1]] = pixel.stripe_avg;
            a[[pixel.x as usize, pixel.y as usize, 2]] = pixel.distance;
            a[[pixel.x as usize, pixel.y as usize, 3]] = pixel.surface_normal.0;
            a[[pixel.x as usize, pixel.y as usize, 4]] = pixel.surface_normal.1;
            a[[pixel.x as usize, pixel.y as usize, 5]] = pixel.iteration;
        }
    }
    // println!(
    //     "Elapsed time: {:.2?}, megapixels: {}",
    //     before.elapsed(),
    //     (&xsize * &ysize / 1000000.),
    // );
    a
}
#[pymodule]
fn qfract(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // wrapper of `calculate_mandelbort`
    #[pyfn(m)]
    #[pyo3(name = "calculate")]
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
    }

    Ok(())
}
