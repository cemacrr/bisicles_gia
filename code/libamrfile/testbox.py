from ctypes import *;
import numpy as np
import matplotlib.pyplot as plt

libamrfile = CDLL("libamrfile.so")


filename="plot.amundsen.2d.hdf5"
status=0
amrID=0
libamrfile.amr_read_file(pointer(c_int(status)), pointer(c_int(amrID)), filename) 
ox = 0
oy = 0
nx = 127
ny = 127
x = np.zeros(nx)
y = np.zeros(ny)
thk = np.asfortranarray(np.zeros((nx,ny)))
lo = np.intc([ox+0,oy+0])
hi = np.intc([ox+nx-1,oy+ny-1])
interp = 0
lev=0
comp=0

libamrfile.amr_read_box_data_2d(pointer(c_int(status)), 
                                thk.ctypes.data_as(POINTER(c_double)), 
                                x.ctypes.data_as(POINTER(c_double)),
                                y.ctypes.data_as(POINTER(c_double)), 
                                pointer(c_int(amrID)), 
                                pointer(c_int(lev)), 
                                lo.ctypes.data_as(POINTER(c_int)), 
                                hi.ctypes.data_as(POINTER(c_int)), 
                                pointer(c_int(comp)), 
                                pointer(c_int(interp)))

libamrfile.amr_free(pointer(c_int(status)), pointer(c_int(amrID))) 
xx, yy = np.meshgrid(x, y)
plt.pcolormesh(thk.transpose())
plt.savefig("test.png")
