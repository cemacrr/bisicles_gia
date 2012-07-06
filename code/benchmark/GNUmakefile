# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase := solverBenchmark

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

TARMACDIR=../interface/libtarmac
XTRALIBFLAGS=-L$(TARMACDIR) -ltarmac

#  have a default for machines not defined here.
# make this a relative path so it will find the Chombo 
# on same level as BISICLES
CHOMBO_HOME = ../../../Chombo/lib

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

SRC_DIR = ../src
UTIL_DIR = ../util 
src_dirs = $(SRC_DIR) $(UTIL_DIR) 
#src_dirs += /home/ranken/util/pglimdev/parallel/src/external_dycore

RUNCOMMAND = 
ifeq ($(MPI),TRUE)
RUNCOMMAND += mpirun -np 2
endif

include $(CHOMBO_HOME)/mk/Make.example

all: all-example $(TARMACDIR)/libtarmac.a

$(TARMACDIR)/libtarmac.a: $(TARMACDIR)/tarmac.f90
	cd $(TARMACDIR); \
	$(MAKE) all


ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif

driver: 
	cd ../exec2D; \
	$(MAKE) driver 

regression: all driver
	$(RUNCOMMAND) ../exec2D/driver$(config).ex inputs.regression
	$(RUNCOMMAND) $(ebase)$(config).ex inputs.regression

doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""
