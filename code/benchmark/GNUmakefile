# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := solverBenchmark
include ../mk/Make.defs

lib:
	cd ../lib; $(MAKE) lib

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)
