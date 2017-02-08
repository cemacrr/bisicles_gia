# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := solverBenchmark
all:  libbike 

include ../mk/Make.defs

libbike:
	cd ../lib; $(MAKE) lib

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)


