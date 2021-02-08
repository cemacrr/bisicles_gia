#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 10:49:03 2020

@author: scornford
"""
from amrfile import io as amrio
amrio.freeAll()
import numpy as np
from osgeo import gdal, osr
from netCDF4 import Dataset
from osgeo.gdalconst import *
import os

def _system(cmd):
    print(cmd)
    os.system(cmd)
   

def nctoamr(nc_file, hdf5_file, var_string):
    _system('nctoamr {} {} {}'.format(nc_file, hdf5_file, var_string))
    
class BisiclesData:
    """
    Store data from a BISICLES output and provide
    derived data (e.g speed from velocity)
    """
    
    def __init__(self,file_name, level=0, origin=(0,0),
                 croplo = (0,0), crophi = (1,1), plot_file = True):
        """

        """
    
        amrid = amrio.load(file_name)
        lo,hi = amrio.queryDomainCorners(amrid, level)

        lo_0 = lo
    
        for dir in [0,1]:
            L = hi[dir] - lo[dir]
            hi[dir] = int( lo_0[dir] + crophi[dir]*L )
            lo[dir] = int( lo_0[dir] + croplo[dir]*L )
            
        iord = 1
        def read(name):
            return amrio.readBox2D(amrid, level, lo, hi, name, iord)
        
        self.x,self.y,self.thk = read('thickness')

        #adjust x,y origin
        self.x += origin[0]
        self.y += origin[1]
        
        x,y,self.usrf = read('Z_surface')
        x,y,self.topg = read('Z_base')
        
        if plot_file:
            #plot.* file data
            x,y,self.xvel = read('xVel')
            x,y,self.yvel = read('yVel')
            x,y,self.beta = read('dragCoef')
            x,y,self.hmu = read('viscosityCoef')
            x,y,self.acab = read('surfaceThicknessSource')
        else:
            #ctrl.* file data
            x,y,self.xvel = read('xVelb')
            x,y,self.yvel = read('yVelb')
            x,y,self.xvel_o = read('xVelo')
            x,y,self.yvel_o = read('yVelo')
            x,y,self.velc = read('velc')
            x,y,self.mucoef = read('muCoef')
            x,y,self.beta = read('Cwshelf')

        amrio.free(amrid)
           
        hf = np.where(self.topg < 0.0, -self.topg*1027.0/917.0, 0.0)
        self.hab = self.thk - hf
        self.speed = np.sqrt(self.xvel**2 + self.yvel**2)
        self.Tb = self.beta * self.speed
        
        if not (plot_file):
            tol = 0.1
            self.speed_o = np.ma.masked_array( np.sqrt(self.xvel_o**2 + self.yvel_o**2), self.velc < tol)
            self.misfit = np.ma.masked_array( (self.speed - self.speed_o)*self.velc, self.velc < tol)
    


        
def write_raster(data,path, EPSG=3031):
    """
    write BisiclesData speeed,Tb to a geotiff path
    """
    N_BANDS = 2
    x = data.x
    y = data.y
    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    dx = x[1]-x[0]
    dataset = driver.Create(path, len(x), len(y), N_BANDS, gdal.GDT_Float64)
    dataset.SetGeoTransform( ( np.min(x), dx, 0, np.max(y), 0, -dx) ) 
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    wkt_projection = srs.ExportToWkt()
    dataset.SetProjection(wkt_projection)
    np.ma.set_fill_value(data,-9999.0)
    
    #speed
    d = np.flipud((data.speed)/365.0)
    dataset.GetRasterBand(1).WriteArray(d)
    dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
    
    
    #basal traction
    d = np.flipud(data.Tb/1.0e3)
    dataset.GetRasterBand(2).WriteArray(d)
    dataset.GetRasterBand(2).SetNoDataValue(-9999.0) 

    dataset.FlushCache()  # Write to disk. 


def read_raster(raster_path):
    """
    Opens a tiff as specified by the user

    Returns an array of the raster with co-oordinates
    """
    from osgeo import gdal,gdalconst  
    import numpy as np
    
    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    src = gdal.Open(raster_path, gdalconst.GA_ReadOnly)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    
    
    data=src.ReadAsArray()
    print("Opened %s" %(raster_path))
    #print(src.GetMetadata())
    
    tol = 1.0e-6
    x = np.arange(ulx,lrx-tol, +xres)
    y = np.arange(lry,uly-tol, -yres)
    
    return x,y,np.flipud(data[:,:])*365 # m/day -> m/a
