import numpy as np
from .qfract import calculate


class Mandelbrot:
    result = None
    def __init__(
        self,
        x1=-2.25,
        x2=0.75,
        y1=-1.5,
        y2=1.5,
        size=1000,
        max_iter=5000,
    ):
        self.x = [x1, x2]
        self.y = [y1, y2]
        self.max_iter = max_iter
        self.size = size

    def light_vector(self,angle):
        v = np.exp(1j*angle*2*np.pi/360)  
        return (v.real,v.imag)

    @classmethod
    def center(cls, x, y, zoom=1, resolution=(1366, 768), reference=3, max_iter=150):
        x_ = 0.5 * zoom * reference
        y_ = 0.5 * zoom * reference / (resolution[0] / resolution[1])
        kwargs = dict(
            x1=x - x_,
            x2=x + x_,
            y1=y - y_,
            y2=y + y_,
            size=resolution[0],
            max_iter=max_iter,
        )
        return cls(**kwargs)

    def calculate(self, stripes=3,anglex= 0.7071067811865476,angley= 0.7071067811865476, method=2):
        result = np.transpose(
            calculate(
                *self.x,
                *self.y,
                self.max_iter,
                self.size,
                anglex,
                angley,
                stripes,
                method
            ),
            axes=(1, 0, 2),
        )
        self.result = result
        return self

    @property
    def lambert_shading(self):
        if self.result is not None:
            return self.result[:, :, 0]
        else:
            raise Exception("No result found")
    @property        
    def stripe_avg(self):
        if self.result is not None:
            return self.result[:, :, 1]
        else:
            raise Exception("No result found")
    @property        
    def distance(self):
        if self.result is not None:
            return self.result[:, :, 2]
        else:
            raise Exception("No result found")
    @property        
    def nx(self):
        if self.result is not None:
            return self.result[:, :, 3]
        else:
            raise Exception("No result found")
    @property        
    def ny(self):
        if self.result is not None:
            return self.result[:, :, 4]
        else:
            raise Exception("No result found")
    @property        
    def iteration(self):
        if self.result is not None:
            return self.result[:, :, 5]
        else:
            raise Exception("No result found")
