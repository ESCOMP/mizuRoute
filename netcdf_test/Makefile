#========================================================================
# Makefile to compile the netCDF test program
#========================================================================
#
# To troubleshoot your paths and setup, type 'make check'
# Change library path HDF, NCDF, LOCAL, SZIP, ZLIB

# Define the Fortran Compiler
FC=

FC_EXE=

MODE=debug

ifeq "$(FC)" "gnu"
  NCDF_PATH=
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
            -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif

ifeq "$(FC)" "intel"
  NCDF_PATH =
  LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif

ifeq "$(FC)" "nvidia"
 NCDF_PATH =
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
 INCNETCDF = -I$(NCDF_PATH)/include
endif

# define flags
ifeq "$(FC)" "gnu"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -fmax-errors=0 -ffree-line-length-none
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none
  endif
endif

ifeq "$(FC)" "intel"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -FR -auto -fltconsistency -fpe0
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -debug all -warn all -check all,nouninit -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency
  endif
endif

ifeq "$(FC)" "nvidia"
  ifeq "$(MODE)" "fast"
    FLAGS = -fast -O3 -Mdclchk $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Mbounds -Mlist -Minfo -Mdclchk
  endif
endif

# define program
SUBROUTINE = example_netcdf.f90
PROGRAM = test_netcdf.f90

# define executable
EX = test_netcdf.exe

#.SUFFIXES: .f .o .f90

# Compile
all: compile link clean

# compile
compile:
	$(FC_EXE) $(FLAGS) -c $(SUBROUTINE) $(PROGRAM) $(INCNETCDF)

# link routines
link:
	$(FC_EXE) *.o $(LIBNETCDF) -o $(EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
