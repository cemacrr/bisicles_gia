# -*- Mode: Makefile -*- 

#set compiler and CHOMBO_HOME
include ../mk/Make.defs

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := solverBenchmark



MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2



##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic  AMRTimeDependent AMRTools BoxTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .
USE_FAS = FALSE
ifeq ($(USE_FAS), TRUE)
#XTRACXXFLAGS = -DFAS_USE_NEW_PROL -DCH_HAS_FAS 
XTRACXXFLAGS = -DCH_USE_FAS 
SRC_DIR = ../src  ../../../AMRFAS/src
else
SRC_DIR = ../src  
endif
UTIL_DIR = ../util 
src_dirs = $(SRC_DIR) $(UTIL_DIR) 
#src_dirs += /home/ranken/util/pglimdev/parallel/src/external_dycore

RUNCOMMAND = 
ifeq ($(MPI),TRUE)
RUNCOMMAND += mpirun -np 2
endif

#include $(CHOMBO_HOME)/mk/Make.example
BISICLES_MK = ../mk/
include $(BISICLES_MK)/Make.example

all: all-example 


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
