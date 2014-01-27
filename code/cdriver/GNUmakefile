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
src_dirs = ../src   

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

bisiclesObjects:
	cd ../exec2D; $(MAKE) nolink


o/$(config)/cwrapper.o: cwrapper.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DCH_LANG_CC -c -o  o/$(config)/cwrapper.o cwrapper.cpp 

libBisicles.a: ../exec2D/o/$(config) o/$(config)/cwrapper.o libChomboLibs.a
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


libs: bisiclesObjects libBisicles.a libChomboLibs.a

ctestwrapper.$(config).ex : bisiclesObjects libBisicles.a libChomboLibs.a
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DCH_LANG_CC -o  ctestwrapper.$(config).ex testwrapper.cpp -L ./ -lBisicles -lChomboLibs $(LIBFLAGS)   

ftestwrapper.$(config).ex : bisiclesObjects libBisicles.a libChomboLibs.a
	$(FC) -g -DCH_LANG_CC -o  ftestwrapper.$(config).ex testwrapper.F90 -L ./ -lBisicles -lChomboLibs $(LIBFLAGS) -lstdc++   


ctestwrapper: ctestwrapper.$(config).ex

ftestwrapper: ftestwrapper.$(config).ex


.PHONY : ctestwrapper ftestwrapper libs
