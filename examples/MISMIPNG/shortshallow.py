import math as m

B0 = -150.0 #vs 0 from gael

B2 = -728.8 #from gael
B4 = 343.91 #from gael
B6 = - 50.57 #from gael

fc = 1.0e+3
dc = 5.0e+2
wc = 16.0e+3 # 10 from Gael

def Bx(x):
    xx = x/300.0e3
    xx2 = xx*xx
    xx4 = xx2*xx2
    xx6 = xx4*xx2
    return B0 + B2*xx2 + B4*xx4 + B6*xx6 

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
  
    if (x < 700.0e+3):
        thickness =  100.0
        
    return thickness

def topography(x,y):
    y = y - 40e+3
    topography = Bx(x) + By(y)
    
    return topography


def accumulation(x,y,t,thck,topg):
    acc = 0.2
    if (x > 700e+3):
        acc = -1000.0
    return acc


def friction(x,y,t,thck,topg):
    friction = 1.0e4
    return friction





