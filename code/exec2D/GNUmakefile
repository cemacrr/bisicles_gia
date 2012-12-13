# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
include ../mk/Make.defs
ebase := driver stats glfaces faces

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2


ifeq ($(TARMAC),TRUE)
TARMACDIR=../interface/libtarmac
XTRALIBFLAGS += -L$(TARMACDIR) -ltarmac
XTRACPPFLAGS += -DTARMAC
endif
ifeq ($(TARMAC),TRUE)
all:  $(TARMACDIR)/libtarmac.a all-example
$(TARMACDIR)/libtarmac.a: $(TARMACDIR)/tarmac.f90
	$(RM) $(patsubst %,%$(config).ex,$(ebase)); \
	echo $(TARMACDIR)
	cd $(TARMACDIR); \
	$(MAKE) all
endif
##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

ifeq ($(FAS), TRUE)
XTRACXXFLAGS = -DFAS_USE_NEW_CF_INTERP -DCH_HAS_FAS 
#XTRACXXFLAGS = -DCH_HAS_FAS 
SRC_DIR = ../src  ../../../AMRFAS/src
else
SRC_DIR = ../src  
endif
UTIL_DIR = ../util 

src_dirs = $(SRC_DIR) $(UTIL_DIR)  
#src_dirs += /home/ranken/util/pglimdev/parallel/src/external_dycore

include $(CHOMBO_HOME)/mk/Make.example


ifneq ($(PYTHON_INC),)
CPPFLAGS+=-I$(PYTHON_INC) -DHAVE_PYTHON	
LIBFLAGS+=$(PYTHON_LIBS)
endif




ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif

nolink: $(shell $(CHOMBO_HOME)/mk/reverse $(LibNames)) $(objects)  
	@echo BISICLES Objects built

doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""
