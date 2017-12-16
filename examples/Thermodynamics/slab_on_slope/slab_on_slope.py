#constant thickness slab on slope
L = 64.0e+3

def thickness(x,y):
    return 1000.0 if x < 0.75*L else 0.0

def topography(x,y):
    return 1000.0 * (1.0 - x/L)

def velocity(x,y,*etc):
    return x / L * 1.0e+3

def stemp(x,y,t,thck,topg,*etc):
    s = thck + topg
    return 273.15 - 0.01 * s

def acab(x,y,t,*etc):
    return 0.5
