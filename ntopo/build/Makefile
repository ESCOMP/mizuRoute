#========================================================================
# Makefile to compile the network topology program
#========================================================================
#
#
#
#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
F_MASTER = /home/mizukami/hydro_nm/nHRU_routing/ntopo.v1/

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)build/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin
#
#========================================================================
# PART 1: Assemble all of the FUSE sub-routines
#========================================================================
# Define subroutines
#
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
#TEMP_MODSUB = $(DATATYPES) $(UTILS) $(NUMREC) $(IO) $(QROUTE)
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(QROUTE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# PART 2: Define the libraries, driver programs, and executables
#========================================================================

# Define the Fortran Compiler
FC  = ifort

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
  HDF_PATH=/opt/hdf5
  ZLIB=/opt/zlib
  NCDF_PATH=/opt/netcdf-4.3.0+gcc-4.7.2
  SZIP=/opt/szip
  LIBNETCDF = -Wl,-rpath,$(HDF_PATH)/lib:$(ZLIB)/lib:$(NCDF_PATH)/lib:$(SZIP)/lib \
  		-L$(NCDF_PATH)/lib -lnetcdff -lnetcdf \
      -L$(HDF_PATH)/lib -lhdf5_hl -lhdf5 \
      -L$(ZLIB) -L$(SZIP)-lsz -lz -lcurl
  INCNETCDF = -I$(NCDF_PATH)/include -I/$(HDF_PATH)/include -I$(SZIP)/include -I$(ZLIB)/include
endif
ifeq "$(FC)" "ifort"
 HDF_PATH  = /usr/local/hdf5
 NCDF_PATH = /usr/local/netcdf-4.3.0+ifort-12.1
 LOCL_PATH = /usr/local
 SZIP = /opt/szip
 ZLIB = /opt/zlib
 XTRA_LIBS = -L$(LOCL_PATH) -L$(SZIP)/lib -L$(ZLIB)/lib -lsz -lz -lcurl
 LIBNETCDF = -Wl,-rpath,$(HDF_PATH)/lib:$(ZLIB)/lib:$(NCDF_PATH)/lib:$(SZIP)/lib \
		-L$(NCDF_PATH)/lib -lnetcdff -lnetcdf -L$(HDF_PATH)/lib -lhdf5_hl -lhdf5 $(XTRA_LIBS)
 INCNETCDF = -I$(NCDF_PATH)/include -I/$(HDF_PATH)/include -I$(SZIP)/include -I$(ZLIB)/include
endif
ifeq "$(FC)" "pgf90"
 LIBNETCDF = -L/usr/local/netcdf4-pgi/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz
 INCNETCDF = -I/usr/local/netcdf4-pgi/include
endif

# Define the driver routine
DRIVER = $(F_KORE_DIR)process_river_topology.f90

# Define the executable
EXE = process_river_topology.exe

#========================================================================
# PART 3: Compile the puppy
#========================================================================

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

#FLAGS = -O3 -W -v

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
