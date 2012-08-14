FC = gfortran

include ../../code/mk/Make.defs

FFLAGS = -g -I. -I$(NETCDF_HOME)/include 
LIBFLAGS = -L$(NETCDF_HOME)/lib -lnetcdf 
LDFLAGS = -g

OBJS = mgrelax.o friction.o

all: $(OBJS)
	$(FC) $(FFLAGS) mgrelax.o friction.o $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FFLAGS) -c  $< -o $@

clean: 
	$(RM) $(OBJS) a.out

realclean:
	$(RM) $(OBJS) a.out *.mod
