FC = gfortran

UNAMEN := $(shell uname -n)
WHO := $(shell whoami)
#netcdf location
ifeq ($(UNAMEN),martins)
NETCDF_HOME = /usr/local/netcdf
endif

ifeq ($(UNAMEN),loren)
NETCDF_HOME =/home/loren/users/dmartin/util/netcdf/netcdf-4.1.2
endif


ifeq ($(WHO),ggslc)

ifeq ($(UNAMEN),porthos.ggy.bris.ac.uk)
NETCDF_HOME = /local/ggslc/BISICLES/netcdf/serial
endif

ifeq ($(UNAMEN),athos.ggy.bris.ac.uk)
NETCDF_HOME = /local/ggslc/BISICLES/netcdf/serial
endif

endif



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
