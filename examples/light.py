import cv2

import numpy as np
from qfract.tools import Mandelbrot

angle = 0
frame = np.ones(shape=[1000,500]).astype(np.uint8)
while True:
    mset = Mandelbrot.center(0,0,resolution=(200,100))
    angles = mset.light_vector(angle)
    angle+=1
    mset.calculate(stripes=3,method=1,anglex=angles[0],angley=angles[1])
    # cv2.imshow('Frame',mset.lambert_shading*mset.stripe_avg)
    frame = mset.lambert_shading
    cv2.imshow('Frame',frame)

    if cv2.waitKey(25) & 0xFF == ord('q'):
        break
