# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

DIM=2
CHOMBO_HOME = ../../../../Chombo/lib
include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config

ebase := control
BUILD_DIR = ../../controlproblem

#set compiler and CHOMBO_HOME

all: 
	cd $(BUILD_DIR); $(MAKE) $(ebase)
	ln -sf $(BUILD_DIR)/$(ebase)$(config).ex ./$(ebase)$(config).ex

example-only: all

realclean:
	$(RM) *.ex *.hdf5 pout.* time.table
