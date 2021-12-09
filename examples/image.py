import numpy as np
from qfract.tools import Mandelbrot
from PIL import Image


x = Mandelbrot(
      -0.57027936,
      -0.5303794,
      -0.6373346,
      -0.6145346,
    max_iter=5000,
    size=1500,
).calculate(3,method=2)


Image.fromarray((255 * x.lambert_shading*x.stripe_avg).astype(np.uint8)).show()
