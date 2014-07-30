#AMRFile hdf5 i/o functions via libamrfile (a C++ shared lib with C interface)

from ctypes import *;
import numpy as np
libamrfile = CDLL("libamrfile.so")

def __error__(status):
    if (status.value != 0):
        print 'error code ',status.value
    

def load(filename):
    status = c_int(-1)
    amrID = c_int(-1)
    libamrfile.amr_read_file(pointer(status), pointer(amrID), filename)
    __error__(status)
    return amrID

def free(amrID):
    status = c_int(-1)
    libamrfile.amr_free(pointer(status), pointer(amrID)) 
    __error__(status)

def freeAll(amrID):
    status = c_int(-1)
    libamrfile.amr_free_all(pointer(c_int(status)))
    __error__(status)

def readBox2D(amrID, level, lo, hi, component, interpolationOrder = 0):
    status = c_int(-1)
    nx = hi[0] - lo[0]+1
    ny = hi[1] - lo[1]+1
    x = np.zeros(nx)
    y = np.zeros(ny)
    v = np.asfortranarray(np.zeros((nx,ny)))
    nplo = np.intc(lo)
    nphi = np.intc(hi)

    libamrfile.amr_read_box_data_2d(pointer(status), 
                                    v.ctypes.data_as(POINTER(c_double)), 
                                    x.ctypes.data_as(POINTER(c_double)),
                                    y.ctypes.data_as(POINTER(c_double)), 
                                    pointer(amrID), 
                                    pointer(c_int(level)), 
                                    nplo.ctypes.data_as(POINTER(c_int)), 
                                    nphi.ctypes.data_as(POINTER(c_int)), 
                                    pointer(c_int(component)), 
                                    pointer(c_int(interpolationOrder)))

    __error__(status)
    return x,y,v.transpose()
