# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := solverBenchmark
include ../mk/Make.defs

all: libbike

libbike:
	cd ../lib; $(MAKE) lib

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)
