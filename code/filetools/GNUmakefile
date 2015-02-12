
# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable



DIM=2

include ../mk/Make.defs

ebase :=  nctoamr amrtotxt amrtoplot flatten extract merge addbox amrtocf stats glfaces faces rescale sum

ifneq ($(PYTHON_INC),)
ebase += pythonf
endif

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

src_dirs =  ../src ../util  src


XTRACPPFLAGS=-I$(NETCDF_HOME)/include
XTRALIBFLAGS=-L$(NETCDF_HOME)/lib $(NETCDF_LIBS)
include $(CHOMBO_HOME)/mk/Make.example
#CPPFLAGS+=-I$(NETCDF_HOME)/include 
#LIBFLAGS+=-L$(NETCDF_HOME)/lib -lnetcdf -lhdf5_hl -lhdf5

ifneq ($(PYTHON_INC),)
CPPFLAGS+=-I$(PYTHON_INC)	
LIBFLAGS+=$(PYTHON_LIBS)
endif
