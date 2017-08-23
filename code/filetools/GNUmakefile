
# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase :=  nctoamr amrtotxt amrtoplot flatten extract merge addbox amrtocf stats glfaces faces rescale sum pythonf
src_dirs = src
all:  clean libbike 

.NOTPARALLEL: all

include ../mk/Make.defs
 

libbike:
	cd ../lib; $(MAKE) lib

ifneq ($(NETCDF_INC), )
CPPFLAGS += -I/$(NETCDF_INC)
endif
LIBFLAGS := -L../lib/ -lBisicles$(config) $(NETCDF_LIBS)  $(LIBFLAGS)
