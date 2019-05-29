import math as m

Ldomain = 640e+3
dxCrse = 8.0e+3

def basal_flux(x,y,*etc):
    #maintain calving front at the bottom and right edges,
    #aligned with a cell face on the coarsesst grid
    huge = 1.0e+4
    f = 0.0
    
    if (x > Ldomain - dxCrse):
        f = -huge
    if (y < dxCrse ):
        f = -huge

    return f


def thickness(x,y,*etc):

    return 1.0


def topography(x,y,*etc):

    return 1.0

def topography_ripple(x,y,*etc):

    xx = 16.0 * m.pi * x / 640.0e+3
    yy = 16.0 * m.pi * y / 640.0e+3
    
    return 0.0 + m.cos(xx)*m.cos(yy)
