#Geometry from Gudmundson et al, TC vol 6, p 1467

import math as m

Bo = 300.0
fc = 5.0e+3
dc = 1.0e+3
wc = 40.0e+3

def Bx(x):
    xx = x/750.0e3
    xx2 = xx*xx
    xx4 = xx2*xx2
    xx6 = xx4*xx2
    return Bo - 2184.8 * xx2 + 1031.72 * xx4 - 151.72 * xx6 

def By(y):
    return dc/(1.0+m.exp(-2.0*(y-wc)/fc)) + dc/(1.0+m.exp(2.0*(y+wc)/fc))

def slope(x):
    return 50.0 - 4.0 * x / 1.0e+3

def width(x):
    return 4.0e+3 + 16.0e+3 * math.exp(-x/16.0e+3)

def tag(x,y,dx,thck,topg):
    tag = -1.0 #don't tag unless...
    return tag

def thickness(x,y):
    thickness = 0.0
  
    if (x < 1800.0e+3):
        thickness = 100.0
        
    return thickness

def topography(x,y):
    y = y - 120e+3
    topography = Bx(x) + By(y)
    
    return topography


def accumulation(x,y,t,thck,topg):
    acc = 0.3
    if (x > 1800e+3):
        acc = -1000.0
    return acc

def friction(x,y,t,thck,topg):
    friction = 24133.38
    return friction


def linear100(x,y,t,thck,topg):
    zb = -thck * 0.9
    a = m.sin(t*2.0/100.0*m.pi)
    a = a**2 * 0.1
    arg = max(1.0,1.0+zb-topg)
    
    return a*zb*m.log(arg)


def linlogB(x,y,t,thck,topg):
    zb = -thck * 0.9
    a = m.sin(t*2.0/100.0*m.pi)
    a = a**2 * 0.025
    arg = max(1.0,1.0+zb-topg)
    
    return a*zb*m.log(arg)
