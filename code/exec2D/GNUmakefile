# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase := driver
all:  libbike clean
driver:  libbike clean

.NOTPARALLEL: all

include ../mk/Make.defs

reallyclean: realclean libclean

libbike:
	cd ../lib; $(MAKE) lib

libclean:
	cd ../lib; $(MAKE) realclean

LIBFLAGS := -L../lib/ -lBisicles$(config)  $(LIBFLAGS)

.PHONY: libbike libclean
