import cv2

import numpy as np
from qfract.tools import Mandelbrot


x,y = -0.1956945679119091,1.1000283477436004

print(x,y)
angle = 0
zoom_increment = 0
while True:
    mset = Mandelbrot.center(x,y,resolution=(300,300),zoom=1/(1.05**zoom_increment),max_iter=50)
    angles = mset.light_vector(angle)
    zoom_increment+=1
    angle+=3
    mset.calculate(stripes=2,method=0,anglex=angles[0],angley=angles[1])
    frame = mset.stripe_avg*mset.lambert_shading
    cv2.imshow('Frame',frame)

    if cv2.waitKey(25) & 0xFF == ord('q'):
        break
