#========================================================================
# Makefile to compile the routing model
#========================================================================
#
#
#
#========================================================================
# PART 0: User configure part 
#========================================================================
# Define fortran compiler - gfortran, ifort or pgf90
FC  =

# Define the compiler exe
FC_EXE =

# define open MP flags
isOpenMP =
FLAGS_OMP =
LIBOPENMP =
 
# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER =

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
  NCDF_PATH = /opt/local
  LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff
  INCNETCDF = -I$(NCDF_PATH)/include
  ifeq "$(isOpenMP)" "yes"
   FLAGS_OMP = -fopenmp
  endif
endif

ifeq "$(FC)" "ifort"
  NCDF_PATH = /opt/netcdf4-intel
  LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff
  INCNETCDF = -I$(NCDF_PATH)/include
  ifeq "$(isOpenMP)" "yes"
   FLAGS_OMP = -fopenmp
  endif
endif

ifeq "$(FC)" "pgf90"
 NCDF_PATH = /opt/netcdf4-pgf90
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
 INCNETCDF = -I$(NCDF_PATH)/include
endif

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS1 = -p -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=all $(FLAGS_OMP)
 FLAGS = -O3 -ffree-line-length-none -fmax-errors=0 $(FLAGS_OMP)
endif
ifeq "$(FC)" "ifort"
 FLAGS1 = -p -g -debug -warn all -check all -noerror_limit -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency $(FLAGS_OMP)
 FLAGS = -O3 -FR -auto -fltconsistency -fpe0 $(FLAGS_OMP)
endif
ifeq "$(FC)" "pgf90"
 FLAGS1 = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS = -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv
endif

#========================================================================
# PART 1: Define directory paths
#========================================================================
# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)build/src/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin

#========================================================================
# PART 2: Assemble all of the sub-routines
#========================================================================
# Define subroutines
#
# data types
DATATYPES = \
    nrtype.f90 \
    public_var.f90 \
    dataTypes.f90 \
    var_lookup.f90 \
    globalData.f90 \
    popMetadat.f90
#
# define utilities
UTILS = \
    nr_utility.f90 \
    ascii_util.f90 \
		time_utils.f90 \
    gamma_func.f90
#
# read/write files
IO = \
    remap.f90 \
    read_control.f90 \
    read_param.f90 \
    read_netcdf.f90 \
    write_netcdf.f90 \
    read_streamSeg.f90 \
    write_streamSeg.f90 \
    read_runoff.f90 \
    getAncillary.f90 \
    read_restart.f90 \
    write_restart.f90 \
    write_simoutput.f90

# perform model routing
QROUTE = \
    process_time.f90 \
    accum_runoff.f90 \
    basinUH.f90 \
    network_topo.f90 \
    pfafstetter.f90 \
    irf_route.f90 \
    process_ntopo.f90 \
    kwt_route.f90

# concatanate model subroutines
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(QROUTE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# PART 2: Define driver programs, and executables
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)route_runoff.f90

# Define the executable
EXE = route_runoff.exe

#========================================================================
# PART 4: Compile the puppy
#========================================================================
# tasks
all: compile install clean

# compile
compile:
	$(FC_EXE) $(FLAGS) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(EXE)

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod
						rm -f *__genmod.f90
#
# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(EXE) $(EXE_PATH)
