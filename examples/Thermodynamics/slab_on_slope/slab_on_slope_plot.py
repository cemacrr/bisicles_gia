from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt


amrio.freeAll()

def readplot(file, nlayer):
    amrID = amrio.load(file)
    level = 0
    #first, work out the domain corners
    lo,hi = amrio.queryDomainCorners(amrID, level)
    hi[1] = 1
    order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
    x,y,thk = amrio.readBox2D(amrID, level, lo, hi, 'thickness', order)

    sigma = np.zeros(nlayer + 2) # center of each layer, surface, base
    dsigma = 1.0 / float(nlayer) # even layers
    sigma[1:nlayer + 1] = np.arange(0.5*dsigma, 1.0, dsigma)
    sigma[nlayer + 1] = 1.0
    T = np.zeros( (len(y),len(x), len(sigma)))
    namebase = 'internalEnergy'
    for l,s in enumerate(sigma):
        
        name = namebase + '{:04d}'.format(l-1)
        if (l == 0):
            name = namebase + 'Surface'
        elif (l == nlayer + 1):
            name = namebase + 'Base'
            
        x,y,E = amrio.readBox2D(amrID, level, lo, hi, name, order)
        E = E / 2009.0
        Tpmp = 273.15 - 9.7456e-8 * 918.0 * 9.81 * thk * s
        T[:,:,l] = np.where(E > Tpmp, Tpmp, E)
    
    amrio.free(amrID)
    
    return x,y,sigma,thk,T
    
      

plt.figure(figsize=(8,4))

color = ['red','brown','orange','blue','magenta','black','cyan']
nlayer = [4,8,16,32,64,128,256]

for n,col in zip(nlayer,color):
    file = 'plot.slab_on_slope_{}.004000.2d.hdf5'.format(n)
    x,y,s,h,T = readplot(file,n)  
    plt.subplot(1,2,1)
    plt.plot(T[0,47,:]-273.15,1-s,'.-',color=col,label=n)
    plt.subplot(1,2,2)
    plt.plot(T[0,47,:]-273.15,1-s,'.-',color=col,label=n)
    
plt.subplot(1,2,1)    
plt.legend()
plt.xlabel(r'$T$($^\circ$C)')
plt.ylabel(r'$1- \sigma$')


plt.subplot(1,2,2)    
plt.ylim((0,0.1))
plt.xlabel(r'$T$($^\circ$C)')

plt.savefig('T_sigma.png')


