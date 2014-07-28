# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
DIM=2

include ../mk/Make.defs

ebase := testlibamrfile
#set compiler and CHOMBO_HOME


##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames := AMRTools BoxTools

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

src_dirs = src 
#src_dirs += /home/ranken/util/pglimdev/parallel/src/external_dycore

include $(CHOMBO_HOME)/mk/Make.example


libamrfile.so: $(objects) testlibamrfile
	$(CXX) -o libamrfile.so -shared  $(objects)  $(CHOMBO_HOME)/libamrtools$(config).a $(CHOMBO_HOME)/libboxtools$(config).a $(CHOMBO_HOME)/libbasetools$(config).a $(HDFLIBFLAGS)

doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""


