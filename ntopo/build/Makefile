#========================================================================
# Makefile to compile the network topology program
#========================================================================
#
#
#
#========================================================================
# PART 0: User configure part 
#========================================================================
# Define fortran compiler - gfortran, ifort or pgf90
FC  = COMPILERNAME

# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER = PARENTDIRECTORY 

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
  NCDF_PATH=/opt/netcdf4-gcc
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
  		-L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq "$(FC)" "ifort"
 NCDF_PATH = /opt/netcdf4-intel
 LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
		-L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
 INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq "$(FC)" "pgf90"
 NCDF_PATH = /opt/netcdf4-pgf90
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
 INCNETCDF = -I$(NCDF_PATH)/include
endif

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS = -Wall -g -ffree-line-length-none
endif
ifeq "$(FC)" "ifort"
 FLAGS1 = -warn all -noerror_limit -FR -O0 -auto -WB -traceback -g -fltconsistency
# FLAGS = -static -debug -warn all -check all,noarg_temp_created -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
 FLAGS = -debug -warn all -check all,noarg_temp_created -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
endif
ifeq "$(FC)" "pgf90"
 FLAGS1 = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv
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
# PART 1: Assemble all of the sub-routines
#========================================================================
# Define subroutines
# data types
DATATYPES = \
		nrtype.f90 \
		dataTypes.f90 \
		var_lookup.f90 \
		nhru2basin.f90 \
		reachparam.f90 \
		lake_param.f90 

# define utilities
UTILS = \
		nr.f90 \
		nrutil.f90 \
		ascii_util.f90

# read/write files
IO = \
		read_streamSeg.f90 \
		write_ntopo.f90

# perform model routing
QROUTE = \
		network_topo.f90 \

# concatanate model subroutines
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(QROUTE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# PART 2: Define driver programs, and executables
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)process_river_topology.f90

# Define the executable
EXE = process_river_topology.exe

#========================================================================
# PART 3: Compile the puppy
#========================================================================
#.SUFFIXES: .f .o .f90

# tasks
all: compile install clean

# compile
compile:
	$(FC) $(FLAGS) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(EXE)

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod
#
# Copy the executable to the bin directory
install:
		mv $(EXE) $(EXE_PATH)
