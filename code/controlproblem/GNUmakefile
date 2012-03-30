# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase := testCGOptimize control
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
TARMACDIR=../interface/libtarmac
XTRALIBFLAGS += -L$(TARMACDIR) -ltarmac
XTRACPPFLAGS += -DTARMAC
endif

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

src_dirs = src $(SRC_DIR) $(UTIL_DIR)  
#src_dirs += /home/ranken/util/pglimdev/parallel/src/external_dycore

include $(CHOMBO_HOME)/mk/Make.example


ifeq ($(TARMAC),TRUE)
all:  $(TARMACDIR)/libtarmac.a all-example
$(TARMACDIR)/libtarmac.a: $(TARMACDIR)/tarmac.f90
	$(RM) $(patsubst %,%$(config).ex,$(ebase)); \
	cd $(TARMACDIR); \
	$(MAKE) all
endif

ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif

doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""
