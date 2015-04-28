# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
include ../mk/Make.defs

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)
DIM=2

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 
src_dirs = ../src ../util

include $(CHOMBO_HOME)/mk/Make.example

VERBOSE=-v
BISICLES_MK = ../mk/
#include ../mk/Make.example

ifneq ($(PYTHON_INC),)
CPPFLAGS+=-I$(PYTHON_INC) -DHAVE_PYTHON	
LIBFLAGS+=$(PYTHON_LIBS)
endif

ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif

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

libChomboLibs$(config).a: $(CHOMBO_HOME)/libamrelliptic$(config).a $(CHOMBO_HOME)/libamrtimedependent$(config).a $(CHOMBO_HOME)/libamrtools$(config).a $(CHOMBO_HOME)/libboxtools$(config).a $(CHOMBO_HOME)/libbasetools$(config).a
	mkdir tmp_lib_objects; \
	cd tmp_lib_objects; \
		pwd; \
		ar -x $(CHOMBO_HOME)/libamrelliptic$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtimedependent$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libboxtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libbasetools$(config).a; 
	ar -cvq libChomboLibs$(config).a  tmp_lib_objects/*.o
	rm tmp_lib_objects/*.o
	rmdir tmp_lib_objects

bisiclesObjects$(config):
	cd ../exec2D; $(MAKE) nolink


o/$(config)/cwrapper.o: cwrapper.cpp 
	$(cxxc) $(CXXFLAGS) $(CPPFLAGS) -DCH_LANG_CC -c -o  o/$(config)/cwrapper.o cwrapper.cpp 

libBisicles$(config).a: ../exec2D/o/$(config) o/$(config)/cwrapper.o libChomboLibs$(config).a
	if test -f libBisicles$(config).a; then rm libBisicles$(config).a; fi
	mkdir -p o/$(config) 
	cp -rp ../exec2D/o/$(config)  o/
	if test -f o/$(config)/driver.o; then \
           mv o/$(config)/driver.o o/$(config)/driver.o.tmp; \
        fi 
	ar -cvq libBisicles$(config).a o/$(config)/*.o 


libs: bisiclesObjects$(config) libBisicles$(config).a libChomboLibs$(config).a




fcc := $(subst FALSE,$(FC),$(subst TRUE,mpif90,$(MPI)))

fcflag := $(subst FALSE,,$(subst TRUE,-DCH_MPI,$(MPI)))
#not sure if this only works with open mpi : we need the C++ mpi lib as well as the fortran lib
lmpicxx := $(subst FALSE,,$(subst TRUE,-lmpi_cxx,$(MPI)))
ftestwrapper.$(config).ex : libs
	$(fcc) $(fcflag) -g -DCH_LANG_CC -o  ftestwrapper.$(config).ex testwrapper.F90 -L ./ -lBisicles$(config) -lChomboLibs$(config) $(LIBFLAGS) -lstdc++ $(lmpicxx)



ftestwrapper:  ftestwrapper.$(config).ex


.PHONY : ctestwrapper ftestwrapper libs
