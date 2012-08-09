import math

def slope(x):
    return 50.0 - 4.0 * x / 1.0e+3

def thickness(x,y):
    thickness = 0.0
  
    if (x < 70.0e+3):
        surface = 300.0 + slope(x)
        thickness = surface - topography(x,y)
        thickness = max(0.0, thickness)

    return thickness

def topography(x,y):
    w = 4.0e+3 + 16.0e+3 * math.exp(-x/16.0e+3)
    y = y - 32.0e+3
    topography = slope(x) + 1.0e+3 - 1.0e+3 * math.exp(-(y/w * y/w))
    
    return topography
