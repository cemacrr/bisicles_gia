#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:37:14 2019

produce initial data Antarctica models
from bedmachine and measures netcdf files.
Output is a netcdf file conatining the 
minimal BISICLES initial data: 

   x,y cartesian coordinates (EPSG 3413)
   topg (bedrock), 
   thk (ice thickness),
   umod (observed speed),
   umodc ( weight w(x,y) in the misfit f_m(x,y) =  w (|u_model| - |u_obs|)^2)
   btrc (initial guess for C in BISICLES inverse problem)
@author: stephen
"""

from netCDF4 import Dataset
import numpy as np
from  scipy.interpolate import RectBivariateSpline
from scipy import ndimage
from osgeo import ogr, osr


def plot(x,y,z,zmin,zmax,name):
    import matplotlib.pyplot as plt
    print ('plotting ' + name)
    plt.figure()
    plt.subplot(1,1,1,aspect='equal')
    plt.pcolormesh(x[::4],y[::4],z[::4,::4])
    plt.colorbar()
    plt.title(name)
    plt.savefig('{}.png'.format(name),dpi=300)

def zeros_2D(x,y):
    return np.zeros((len(y),len(x)))


def read_umod_mouginot(x,y,measures_nc):
    # compute |u| from the Mouginot data (450 m res)
    ncu = Dataset(measures_nc, 'r')
    xu = ncu.variables["x"][:]
    yu = np.flipud(ncu.variables["y"][:])
    vx = np.flipud( ncu.variables["VX"][:,:] ) # image data - top left to bottom right
    vy = np.flipud( ncu.variables["VY"][:,:] )
    umod = np.sqrt(vx**2 + vy**2)
    print ('max |u| = {} m/a '.format(np.max(umod)))
    print ('umod.shape = {} xu.shape = {} yu.shape = {}'.format(umod.shape,xu.shape,yu.shape))
    print ('interpolating ...')
    #interpolation of velocity data onto bedmachine grid
    spl = RectBivariateSpline(yu,xu,umod,kx=1,ky=1)    
    return spl(y,x) 


def preprocess(output_nc, bedmachine_nc, measures_nc):


    C_MAX = 1.0e+4 # maximum value for C
    C_MIN = 1.0e+2 # minimum value for C

    #desired dimensions
    NX = 6144*2
    NY = 6144*2
    
    ncbm = Dataset(bedmachine_nc, 'r')
    xbm = ncbm.variables["x"][:]
    ybm = np.flipud(ncbm.variables["y"][:])

    print ('xbm[0] = {}, ybm[0] = {}'.format(xbm[0],ybm[0]))

    #bed machine data dimensions
    dx = xbm[1] - xbm[0]

    #bedmachine data
    topg =  np.flipud(ncbm.variables["bed"][:,:])
    thk =  np.flipud(ncbm.variables["thickness"][:,:])
    usrf_bm =  np.flipud(ncbm.variables["surface"][:,:])
    mask = np.flipud(ncbm.variables["mask"][:,:])

    #raise lake vostok
    topg = np.where(mask == 4, usrf_bm - thk, topg)

    #speed data
    umod = read_umod_mouginot(xbm,ybm,measures_nc)
                    
    #trim to desired size
    iy0 = int( (len(xbm) - NX)/2 ) + 1
    ix0 = int( (len(ybm) - NY)/2 ) + 1
    s = (  slice( iy0 , iy0 + NY ), slice ( ix0, ix0 + NX) )
    x = xbm[s[1]]
    y = ybm[s[0]]
    thk = thk[s]
    topg = topg[s]
    umod = umod[s]
    usrf_bm = usrf_bm[s]
    nx = NX
    ny = NY

    #dependents
    rhoi = 917.0
    rhoo = 1027.0
    sg = topg + thk
    sf = (1.0 - rhoi/rhoo)*thk
    eps = 1.0e-6
    grounded = np.logical_and( thk > eps, sg + eps > sf)
    usrf = np.where( grounded, sg, sf )
    plot(x,y,np.log(np.abs(usrf-usrf_bm)),-3,3,'surface_discrepancy')


    print ('umod c ...')
    #umodc is the weight w(x,y) in the misfit f_m(x,y) =  w (|u_model| - |u_obs|)^2
    umodc = np.where(umod > 1.0, 1.0, 0.0)
    umodc = np.where(thk > 10.0, umodc, 0.0)

    #surface gradient
    print ('grad s ...')
    usrf = ndimage.gaussian_filter(usrf, 4) # smooth
    grads = zeros_2D(x,y)
    grads[1:ny-1,1:nx-1] = 0.5 / dx *  np.sqrt(
        (usrf[1:ny-1,0:nx-2] - usrf[1:ny-1,2:nx])**2 + 
        (usrf[0:ny-2,1:nx-1] - usrf[2:ny,1:nx-1])**2 )
 
    #initial guess for C
    print ('btrc...')
    btrc = rhoi * 9.81 * grads * thk / (umod + 1.0)
    btrc = np.where(umod > 1, btrc, C_MAX)
    btrc = np.where(btrc < C_MAX, btrc, C_MAX)
    btrc = np.where(btrc > C_MIN, btrc, C_MIN)
    #smooth with slippy bias
    print ('    ...filtering')
    btrcs = ndimage.minimum_filter(btrc, 8)
    btrcs = ndimage.gaussian_filter(btrcs, 32)
    btrc = np.where(btrc < btrcs, btrc, btrcs) # retain slippy spots
    
    #no ice value for C
    btrc = np.where(thk > 0, btrc, 100.0)
    
    #ouput netcdf
    print ('writing ...')
    ncout = Dataset(output_nc,'w')
    #dimensions
    xdim = ncout.createDimension('x',size=NX)
    ydim = ncout.createDimension('y',size=NY)
    #var defs


    crsv =  ncout.createVariable('crs','int')
    #projection information. too lazy to read from files :)
    #        char mapping ;
    #                mapping:grid_mapping_name = "polar_stereographic" ;
    #                mapping:latitude_of_projection_origin = -90. ;
    #                mapping:standard_parallel = -71. ;
    #               mapping:straight_vertical_longitude_from_pole = 0. ;
    #                mapping:semi_major_axis = 6378273. ;
    #                mapping:inverse_flattening = 298.27940504282 ;
    #                mapping:false_easting = 0. ;
    EPSG = 3031
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(EPSG)
    crs_wkt = crs.ExportToWkt()
    ncout.setncattr('spatial_ref',crs_wkt)
    ncout.setncattr('Conventions','CF-1.7') 
    crsv.setncattr('EPSG',int(EPSG))
    crsv.setncattr('crs_wkt',crs_wkt)
    crsv.setncattr('grid_mapping_name','polar_stereographic')
    crsv.setncattr('latitude_of_projection_origin', -90.0)
    crsv.setncattr('straight_vertical_longitude_from_pole', 0.0)
    crsv.setncattr('scale_factor',1.0)
    crsv.setncattr('standard_parallel',-71.0)
    crsv.setncattr('false_easting',0.0)
    crsv.setncattr('false_northing',0.0)
    
    
    xv = ncout.createVariable('x','f8',('x'))
    xv.setncattr('standard_name','projection_x_coordinate')
    xv.setncattr('units','meter')
    
    yv = ncout.createVariable('y','f8',('y'))
    yv.setncattr('standard_name','projection_y_coordinate')
    yv.setncattr('units','meter')
    
    def create2D(name):
        v = ncout.createVariable(name,'f8',('y','x'))
        v.setncattr('grid_mapping','crs')
        return v

    topgv = create2D('topg')
    thkv = create2D('thk')
    umodv = create2D('umod')
    umodcv = create2D('umodc')
    btrcv = create2D('btrc')
    
    #data
    xv[:] = x
    yv[:] = y
    topgv[:,:] = topg
    thkv[:,:] = thk
    umodv[:,:] = umod
    umodcv[:,:] = umodc
    btrcv[:,:] = btrc

    ncout.close()

    dx = x[1] - x[0]
    print( ' {} < x < {} '.format(np.min(x) - 0.5 * dx, np.max(x) + 0.5*dx))
    dy = y[1] - y[0]
    print( ' {} < y < {} '.format(np.min(y) - 0.5 * dy, np.max(y) + 0.5*dy))

