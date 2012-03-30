# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase := bisicles
#ebase := driver 

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2

TRAPFPE=FALSE
ifeq ($(TRAPFPE),TRUE)	
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += -DTRAP_FPE
endif

# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += #-DPAPI 
PAPILIBS = #/usr/local/lib/libpapi.a
PAPILIBS = #/usr/common/usg/papi/2.1/lib/libpapi.a -L/usr/lpp/pmtoolkit/lib -lpmapi

ifeq ($(TARMAC),TRUE)
TARMACDIR=$(PWD)/libtarmac
XTRALIBFLAGS += -L$(TARMACDIR) -ltarmac
XTRACPPFLAGS += -DTARMAC
endif

#  have a default for machines not defined here.
# make this a relative path so it will find the Chombo 
# on same level as BISICLES
#PWD = ./
CHOMBO_HOME = $(PWD)/../../../Chombo/lib
GLIMMER_HOME = $(PWD)/../../../gc1/parallel



include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config

#cases := $(shell echo case*)

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic  AMRTimeDependent AMRTools BoxTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

DYCORE_DIR = $(PWD)/libdycore
TARMAC_DIR = $(PWD)/libtarmac

SRC_DIR = ../src
UTIL_DIR = ../util



src_dirs = $(SRC_DIR) $(UTIL_DIR)  
src_dirs += $(DYCORE_DIR)
src_dirs += $(DYCORE_DIR)/BISICLES
src_dirs += $(DYCORE_DIR)/Ymir



include $(CHOMBO_HOME)/mk/Make.example

ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif

.PHONY: libdycore glimmer libtarmac

bisicles: copyFiles glimmer libdycore libtarmac libChomboLibs libBisicles glimmer simple_bisicles
	@echo "BISICLES-Glimmer is up to date"

copyFiles: 
	cd $(DYCORE_DIR); \
	$(MAKE) copyFiles

# note that we force f77 to be whatever we're using for FC, since otherwise
# the SLAP stuff has a habit of using g77 
glimmer: 
	@echo "building glimmer"; \
	cd $(GLIMMER_HOME); \
	$(MAKE) F77=$(FC)

libdycore: 
	@echo "building dycore library"; \
	cd $(DYCORE_DIR); \
	$(MAKE) FC=$(FC) CXX=$(CXX) BISICLES

libtarmac:
	@echo "building tarmac library"; \
	cd $(TARMAC_DIR); \
	$(MAKE) FC=$(FC) libtarmac.a

libChomboLibs: lib
	@echo "building libChomboLibs.a"; \
	mkdir tmp_lib_objects; \
	cd tmp_lib_objects; \
		pwd; \
		ar -x $(CHOMBO_HOME)/libamrelliptic$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtimedependent$(config).a; \
		ar -x $(CHOMBO_HOME)/libamrtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libboxtools$(config).a; \
		ar -x $(CHOMBO_HOME)/libbasetools$(config).a; 

	if test -f libChomboLibs.a; then rm libChomboLibs.a; fi
	ar -cvr libChomboLibs.a  tmp_lib_objects/*.o
	rm tmp_lib_objects/*.o
	rmdir tmp_lib_objects

libBisicles: $(objects) o/$(config)/bike_driver.o
	@echo "building libBisicles.a"; \
	if test -f libBisicles.a; then rm libBisicles.a; fi
	if test -f o/$(config)/driver.o; then \
           mv o/$(config)/driver.o o/$(config)/driver.o.tmp; \
        fi \

	ar -cvr libBisicles.a o/$(config)/*.o 

	if test -f o/$(config)/driver.o.tmp; then \
	   mv o/$(config)/driver.o.tmp o/$(config)/driver.o; \
	fi



simple_bisicles:
	@echo "building simple_bisicles"; \
	cd $(GLIMMER_HOME)/src/fortran; \
	$(MAKE) -f Makefile.dycore FC=$(FC) simple_bisicles; \
	$(MAKE) -f Makefile.dycore FC=$(FC) install

squeakyclean: 
	$(MAKE) realclean NODEPENDS=TRUE ; \
	cd $(DYCORE_DIR) && $(MAKE) clean; \
	cd $(TARMAC_DIR) && $(MAKE) clean


doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""

