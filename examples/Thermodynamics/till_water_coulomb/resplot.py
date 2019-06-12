from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

amrio.freeAll()

def readplot(file, nlev):
    amrID = amrio.load(file)
    level = nlev
    #first, work out the domain corners
    lo,hi = amrio.queryDomainCorners(amrID, level)
    order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
    x,y,thk = amrio.readBox2D(amrID, level, lo, hi, 'thickness', order)
    x,y,usrf = amrio.readBox2D(amrID, level, lo, hi, 'Z_surface', order)
    x,y,ux = amrio.readBox2D(amrID, level, lo, hi, 'xVel', order)
    x,y,uy = amrio.readBox2D(amrID, level, lo, hi, 'yVel', order)
    u = np.sqrt(ux**2 + uy**2) + 1.0e-10
    
    amrio.free(amrID)
    
    return x*1e-3,y*1e-3,thk,usrf,u
    
      

plt.figure(figsize=(6,8))
n = 3
m = 2

for lev in [0,1,2,3,4]:

    plt.subplot(n,m,lev+1,aspect='equal')
    file = sorted(glob.glob('v1_480/plot.twc.{}lev.*.2d.hdf5'.format(lev)))[-1]
    print file
    x,y,h,s,u = readplot(file,lev) 
    plt.pcolormesh(x,y,np.log10(u),vmin = 0, vmax = 3, cmap = 'RdYlBu_r')
    plt.title(r'$ \Delta x_L =  $ {} km'.format(x[1] - x[0])) 
    plt.xticks([])
plt.savefig('mesh.png',dpi=300)


