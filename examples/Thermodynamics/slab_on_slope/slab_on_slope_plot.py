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
    
      

plt.figure(figsize=(12,12))

color = ['red','brown','purple','orange','blue','magenta','cyan','black']
nlayer = np.array([4,8,16,32,64,128,256,512])
Tmin = np.zeros(len(nlayer))
nl = 0 
for n,col in zip(nlayer,color):
    file = 'plot.slab_on_slope_{}.004000.2d.hdf5'.format(n)
    x,y,s,h,T = readplot(file,n)  
    TT = T[0,47,:]-273.15
    plt.subplot(2,2,1)
    plt.plot(TT,1-s,'.-',color=col,label=n)
    
    plt.subplot(2,2,2)
    plt.plot(TT,1-s,'.-',color=col,label=n)
    
    Tmin[nl] = np.min(TT)
    nl = nl + 1

plt.subplot(2,2,3)
Tslice = np.flipud( (T[0,:,:]).transpose() ) - 273.15
plt.pcolormesh(Tslice)
plt.colorbar()
plt.xlim((0,48))
plt.ylim((0,np.max(nlayer)))

plt.subplot(2,2,4)
N = len(Tmin)
plt.loglog((nlayer[0:N-1]),np.abs(Tmin[0:N-1]-Tmin[1:N]))
plt.xlabel('n layers')
plt.ylabel('|T_min_n - T_min_n+1|' )
plt.loglog(nlayer,1/nlayer)

   
plt.subplot(2,2,1)    
plt.legend()
plt.xlabel(r'$T$($^\circ$C)')
plt.ylabel(r'$1- \sigma$')


plt.subplot(2,2,2)    
plt.ylim((-0.01,0.1))
plt.xlabel(r'$T$($^\circ$C)')

plt.savefig('T_sigma.png',dpi=300)


