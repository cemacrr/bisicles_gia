# -*- Mode: Makefile -*- 

#creates a library of column physics routine with a fortran77 interface 


FFLAGS=-O3
LIB=libtarmac.a
OBJS=tarmac.o


all: $(LIB) 

$(OBJS): %.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	ar crus $@ $< 

.PHONY: all clean

clean:
	\rm -f $(OBJS) $(LIB) *.mod
