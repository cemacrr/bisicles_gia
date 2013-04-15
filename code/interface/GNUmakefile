# This makefile is used primarily to make bike_driver, then build a BISICLES library
# for use by the Fortran code, Glimmer.
# Doug Ranken 5/14/10
  
#CHOMBO_HOME = /home/ranken/util/Chombo/lib
CHOMBO_HOME = $(PWD)/../../../Chombo/lib

CISM_HOME = $(PWD)/../../../seacism

#LIBDYCORE = libdycore
LIBDYCORE = $(CISM_HOME)/libdycore

include $(CHOMBO_HOME)/mk/Make.defs.local
include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config
# include $(CHOMBO_HOME)/../example/AMRINS/util/CellToEdge.H

#HDF5 = /home/ranken/util/hdf/5-1.6.10-linux-x86_64-static
#HDF5 = /home/ranken/util/hdf/hdf5-1.8.4/hdf5

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic  AMRTimeDependent AMRTools BoxTools BaseTools

SRC_DIR = ../src
UTIL_DIR = ../util 

src_dirs = $(SRC_DIR) $(UTIL_DIR)  

XTRACPPFLAGS+=-I$(LIBDYCORE)

#FC = gfortran
#FC = f95
#CC = CC

# link arguments used in Chombo makefile
# LINK_ARGS = -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-25 -m64
 LINK_ARGS = -m64 -lz 
 
# need?: -lstdc++


driver: driver.o libChomboLibs.a libBisicles.a
	$(FC) -o driver driver.o \
        $(LINK_ARGS) \
	libBisicles.a \
	libChomboLibs.a \
	$(HDF5)/lib/libhdf5.a \
	/home/ranken/util/hdf/szip-2.1/src/.libs/libsz.a \
	-lstdc++

#libs: libChomboLibs libBisicles libdycore

$(CHOMBO_HOME)/libamrelliptic$(config).a :
	cd $(CHOMBO_HOME); \
	$(MAKE) amrelliptic

$(CHOMBO_HOME)/libamrtimedependent$(config).a :
	cd $(CHOMBO_HOME); \
	$(MAKE) amrtimedependent

$(CHOMBO_HOME)/libamrtools$(config).a :
	cd $(CHOMBO_HOME); \
	$(MAKE) amrtools

$(CHOMBO_HOME)/libboxtools$(config).a :
	cd $(CHOMBO_HOME); \
	$(MAKE) boxtools

$(CHOMBO_HOME)/libbasetools$(config).a :
	cd $(CHOMBO_HOME); \
	$(MAKE) basetools

libChomboLibs.a: $(CHOMBO_HOME)/libamrelliptic$(config).a $(CHOMBO_HOME)/libamrtimedependent$(config).a $(CHOMBO_HOME)/libamrtools$(config).a $(CHOMBO_HOME)/libboxtools$(config).a $(CHOMBO_HOME)/libbasetools$(config).a
	mkdir tmp_lib_objects; \
	cd tmp_lib_objects; \
		pwd; \
		ar -x $(CHOMBO_HOME)/libamrelliptic$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtimedependent$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libboxtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libbasetools$(config).a; 
	ar -cvq libChomboLibs.a  tmp_lib_objects/*.o
	rm tmp_lib_objects/*.o
	rmdir tmp_lib_objects

../exec2D/driver.o:
../exec2D/driver$(config).ex:

bisiclesObjects:
	cd ../exec2D; $(MAKE) nolink

libs: bisiclesObjects libBisicles.a libChomboLibs.a

#dependency on ../exec2D/o/$(config) is designed to catch updated .o files
# so far, it seems to do the trick (DFM -- 7/9/12)
libBisicles.a: o/$(config)/bike_driver.o ../exec2D/o/$(config) libChomboLibs.a
	if test -f libBisicles.a; then rm libBisicles.a; fi
	mkdir -p o/$(config) 
	cp -rp ../exec2D/o/$(config)  o/
	if test -f o/$(config)/driver.o; then \
           mv o/$(config)/driver.o o/$(config)/driver.o.tmp; \
        fi 
	if test -f o/$(config)/faces.o; then \
           mv o/$(config)/faces.o o/$(config)/faces.o.tmp; \
        fi 
	if test -f o/$(config)/glfaces.o; then \
           mv o/$(config)/glfaces.o o/$(config)/glfaces.o.tmp; \
        fi 
	if test -f o/$(config)/stats.o; then \
           mv o/$(config)/stats.o o/$(config)/stats.o.tmp; \
        fi 
	ar -cvq libBisicles.a o/$(config)/*.o 
	if test -f o/$(config)/driver.o.tmp; then \
	   mv o/$(config)/driver.o.tmp o/$(config)/driver.o; \
	fi 
	if test -f o/$(config)/faces.o.tmp; then \
	   mv o/$(config)/faces.o.tmp o/$(config)/faces.o; \
	fi 
	if test -f o/$(config)/glfaces.o.tmp; then \
	   mv o/$(config)/glfaces.o.tmp o/$(config)/glfaces.o; \
	fi 
	if test -f o/$(config)/stats.o.tmp; then \
	   mv o/$(config)/stats.o.tmp o/$(config)/stats.o; \
	fi 



#driver.o: driver.cpp
#	pwd && $(CC) -c -g driver.cpp -I $(CHOMBO_HOME)/include -I ../src \
#	-DCH_SPACEDIM=2 -DCH_USE_DOUBLE -DCH_Linux -DCH_LANG_CC \
#	-DCH_USE_SETVAL -DCH_USE_COMPLEX -DUSE_MEMORY_TRACKING \
#	-DCH_USE_64 -DCH_MPI


bike_driver_call: bike_driver_call.o
	$(FC) -o bike_driver_call bike_driver_call.o \
	bike_driver.o \
        $(LINK_ARGS) \
	libBisicles.a \
	libChomboLibs.a \
	$(HDF5)/lib/libhdf5.a \
	-lz \
	-lstdc++	

bike_driver_call.o: bike_driver_call.F90
	$(FC) -c bike_driver_call.F90 -o bike_driver_call.o


# o/$(config)/bike_driver.o: bike_driver.cpp
# 	$(CC) -c -g bike_driver.cpp -I $(CHOMBO_HOME)/include -I ./ -I ../src \
# 	-I $(LIBDYCORE) \
# 	-DCH_SPACEDIM=2 -DCH_USE_DOUBLE -DCH_Linux -DCH_LANG_CC \
# 	-DCH_USE_SETVAL -DCH_USE_COMPLEX -DUSE_MEMORY_TRACKING \
# 	-DCH_USE_64 -DCH_MPI

#libBisicles.a: bike_driver.o
#	if test -f libBisicles.a; then rm libBisicles.a; fi
#	if test -f o/$(config)/driver.o; then \
#           mv o/$(config)/driver.o o/$(config)/driver.o.tmp; \
#        fi \

#	ar -cvq libBisicles.a o/$(config)/*.o 

#	if test -f driver.o; then \
#	   cp driver.o o/$(config); \
#	fi
 
testGradientCall.o: testGradientCall.F90
	$(FC) -c testGradientCall.F90 -o testGradientCall.o

#	@echo "$(config)"

#include $(CHOMBO_HOME)/mk/Make.example
BISICLES_MK = ../mk/
include $(BISICLES_MK)/Make.example
