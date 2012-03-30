FC = gfortran

UNAMEN := $(shell uname -n)

#netcdf location
ifeq ($(UNAMEN),martins)
NETCDF_HOME = /usr/local/netcdf
endif

ifeq ($(UNAMEN),loren)
NETCDF_HOME =/home/loren/users/dmartin/util/netcdf/netcdf-4.1.2
endif

ifeq ($(UNAMEN),geog-bsq92.ggy.bris.ac.uk)
NETCDF_HOME = /local/GlimmerBisicles/netcdf
endif

ifeq ($(UNAMEN),stephmbp)
NETCDF_HOME = /home/steph/LBLvisit/GlimmerBISICLES/netcdf/serial/
endif

FFLAGS = -O2 -I. -I$(NETCDF_HOME)/include 
LIBFLAGS = -L$(NETCDF_HOME)/lib -lnetcdf 
LDFLAGS = -O2

OBJS = mgrelax.o pig.o

all: $(OBJS)
	$(FC) $(FFLAGS) mgrelax.o pig.o $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FFLAGS) -c  $< -o $@

clean: 
	$(RM) $(OBJS) a.out

realclean:
	$(RM) $(OBJS) a.out *.mod
