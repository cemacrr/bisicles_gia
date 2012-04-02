# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase :=  nctoamr amrtotxt flatten extract pythonf merge
#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2

#  have a default for machines not defined here.
# make this a relative path so it will find the Chombo 
# on same level as BISICLES
CHOMBO_HOME = ../../../Chombo/lib

include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config



#CPPFLAGS+=-I($(NETCDF_HOME))/include
#LDFLAGS+=-L($(NETCDF_HOME))/lib -lnetcdf
#cases := $(shell echo case*)

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTools BoxTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

src_dirs = src 

include $(CHOMBO_HOME)/mk/Make.example

UNAMEN := $(shell uname -n)
WHO := $(shell whoami)
#netcdf location
ifeq ($(UNAMEN),martins)
NETCDF_HOME = /usr/local/netcdf
endif

ifeq ($(UNAMEN),loren)
NETCDF_HOME =/home/loren/users/dmartin/util/netcdf/netcdf-4.1.2
endif


ifeq ($(UNAMEN),porthos.ggy.bris.ac.uk)
ifeq ($(WHO),ggslc)
NETCDF_HOME = /local/ggslc/BISICLES/netcdf/serial
endif
endif


ifeq ($(UNAMEN),geog-bsq92.ggy.bris.ac.uk)
NETCDF_HOME = /local/GlimmerBisicles/netcdf
endif

ifeq ($(UNAMEN),stephmbp)
NETCDF_HOME = /home/steph/LBLvisit/GlimmerBISICLES/netcdf/serial/
endif

CPPFLAGS+=-I$(NETCDF_HOME)/include 
LIBFLAGS+=-L$(NETCDF_HOME)/lib -lnetcdf

PYTHON_VERSION=2.4

ifeq ($(UNAMEN),loren)
PYTHON_VERSION=2.6
endif

ifeq ($(UNAMEN),martins)
PYTHON_VERSION=2.6
CPPFLAGS+=-I/home/dan/install/visit/python/2.6.4/linux-x86_64_gcc-4.4/include/python2.6
LIBFLAGS+= -L/home/dan/install/visit/python/2.6.4/linux-x86_64_gcc-4.4/lib/
endif

ifeq ($(PYTHON_VERSION),2.6)
PYTHON_INC=/usr/include/python2.6
CPPFLAGS+=-I$(PYTHON_INC)
LIBFLAGS+= -L/usr/lib -lpython2.6 
else
PYTHON_INC=/usr/include/python2.4
CPPFLAGS+=-I$(PYTHON_INC)
LIBFLAGS+=-lpython2.4
endif
