import cv2

from qfract.tools import Mandelbrot

resolution = (1000,700)
stripes = 3
#Center and zoom start position

# Spiral
start_center,zoom = [-0.20456008226012767,-0.6539432845897316], 163
#Common start
start_center,zoom = [-0.7,0], -2



# If none, it will be adapted based on zoom, otherwise use max_iter 100-5000
max_iter= None

# zoom on click
zoom_increment=15

instructions = """
Instructions:

Q                 : Exit 
Left mouse click  : Zoom in
Middle mouse click: Zoom out

For SPEED decrese max_iter (lowers quality with increasing zoom)
"""

print(instructions)

def adapt_max_iter(max_iter,zoom):
    if max_iter is None:
        if zoom <=30:
            max_iter = 100
        if zoom > 30:
            max_iter = 200
        if zoom > 50:
            max_iter = 500
        if zoom > 80:
            max_iter = 800
        if zoom > 110:
            max_iter = 1200
        if zoom > 140:
            max_iter = 2500
        if zoom > 180:
            max_iter = 5000
        return max_iter
    else:
        return max_iter

def get_center(box,mouse_pos,max_width):
    ratio = mouse_pos/max_width
    return box[0]+ratio*(box[1]-box[0])

def capture_event(event,x,y,flags,params):
    global angle,zoom,resolution, mset, max_iter
    xx = mset.x
    yy = mset.y
    modified_max_iter = adapt_max_iter(max_iter,zoom)
    x_center = get_center(xx,x,resolution[0])
    y_center = get_center(yy,y,resolution[1])

    if event==cv2.EVENT_LBUTTONDOWN:
        zoom+=zoom_increment
    if event==cv2.EVENT_MBUTTONDOWN:
        zoom-=zoom_increment
    if event==cv2.EVENT_LBUTTONDOWN or event==cv2.EVENT_MBUTTONDOWN:
        mset = Mandelbrot.center(x_center,y_center,resolution=resolution,zoom=1/(1.1**zoom),max_iter=modified_max_iter)
        mset.calculate(stripes=stripes,method=0)
        print(f"Location [x0,x1,y0,y1]: {[*mset.x,*mset.y]}")
        print(f"center [x0,y0] - zoom: [{x_center},{y_center}], {zoom}")
        print(f"max_iter: {modified_max_iter}")

        frame = mset.stripe_avg*mset.lambert_shading
        cv2.imshow(window,frame)




# set the name to the window
window = 'Simple Explorer'


mset = Mandelbrot.center(*start_center,resolution=resolution,zoom=1/(1.1**zoom),max_iter=adapt_max_iter(max_iter,zoom))
mset.calculate(stripes=stripes,method=0)
frame = mset.stripe_avg*mset.lambert_shading
cv2.namedWindow(window)
cv2.imshow(window,frame)
cv2.setMouseCallback(window,capture_event)

while True:
    if cv2.waitKey(25) & 0xFF == ord('q'):
        break
cv2.destroyAllWindows()