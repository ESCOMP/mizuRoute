#========================================================================
# Makefile to compile the routing model
#========================================================================
#
# Need to edit FC, FC_EXE, EXE, isOpenMP (if OMP is used), F_MASTER, NCDF_PATH
# Feel free to modify compiler flag (FLAGS)
# Do not leave space after your edits
#
#========================================================================
# User configure part
#========================================================================
# Define fortran compiler - gnu, intel or pgi
FC  =

# Define the compiler exe, e.g., gnu=>gfortran, intel=>ifort, pgi=>pgf90
FC_EXE =

# Define the executable
EXE =

# Define optional setting
# fast:      Enables optimizations
# debug:     Minimum debug options, still
# profile:   Enables profiling
MODE = debug

# define open MP option (put yes to activate OMP)
isOpenMP =

# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER =

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gnu"
 NCDF_PATH =
else ifeq "$(FC)" "intel"
 NCDF_PATH =
else ifeq "$(FC)" "pgi"
 NCDF_PATH =
else
 $(error FC is not set correctly)
endif

LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
            -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

#========================================================================
# generate a Git version string
#========================================================================
VERSION = $(shell git tag | tail -n 1)
GITBRCH = $(shell git describe --long --all --always | sed -e's/heads\///')
GITHASH = $(shell git rev-parse HEAD)

#========================================================================
# Define flags
#========================================================================
FLAGS_OMP =
LIBOPENMP =
ifeq "$(isOpenMP)" "yes"
 ifeq "$(FC)" "pgi"
  FLAGS_OMP = -mp
 endif
 ifeq "$(FC)" "gnu"
  FLAGS_OMP = -fopenmp
 endif
 ifeq "$(FC)" "intel"
  FLAGS_OMP = -qopenmp
 endif
endif

ifeq "$(FC)" "gnu"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -fmax-errors=0 -ffree-line-length-none $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
endif

ifeq "$(FC)" "intel"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -FR -auto -fltconsistency -fpe0 $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -debug all -warn all -check all -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -debug -warn all -check all -noerror_limit -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency $(FLAGS_OMP) -cpp -DVERSION=\"$(VERSION)\" -DBRANCH=\"$(GITBRCH)\" -DHASH=\"$(GITHASH)\"
  endif
endif

ifeq "$(FC)" "pgi"
  ifeq "$(MODE)" "fast"
    FLAGS = -fast -O3 -Mdclchk $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Mbounds -Mlist -Minfo -Mdclchk $(FLAGS_OMP)
  endif
endif

#========================================================================
# Define directory paths
#========================================================================
# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)build/src/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin

#========================================================================
# Assemble all of the sub-routines
#========================================================================
# Define subroutines
#
# define utilities
UTILS = \
    nrtype.f90 \
    nr_utility.f90 \
    ascii_util.f90 \
    ncio_utils.f90 \
    gamma_func.f90
# data types
DATATYPES = \
    public_var.f90 \
    dataTypes.f90 \
    var_lookup.f90 \
    time_utils.f90 \
    datetime_data.f90 \
    csv_data.f90 \
    gageMeta_data.f90 \
    obs_data.f90 \
    globalData.f90 \
    popMetadat.f90 \
    allocation.f90
# initialization
INIT = \
    network_topo.f90 \
    process_param.f90 \
    process_ntopo.f90 \
    pfafstetter.f90 \
    domain_decomposition.f90
# read/write files
IO = \
    remap.f90 \
    read_control.f90 \
    read_param.f90 \
    read_streamSeg.f90 \
    write_streamSeg.f90 \
    read_runoff.f90 \
    get_basin_runoff.f90 \
    read_remap.f90 \
    read_restart.f90 \
    write_restart.f90 \
    write_simoutput.f90
# CORE
CORE = \
    model_finalize.f90 \
    data_assimilation.f90 \
    accum_runoff.f90 \
    basinUH.f90 \
    irf_route.f90 \
    kwt_route.f90 \
    dfw_route.f90 \
    kw_route.f90 \
    mc_route.f90 \
    main_route.f90 \
    model_setup.f90

# concatanate model subroutines
TEMP_MODSUB = $(UTILS) $(DATATYPES) $(INIT) $(IO) $(CORE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# Define driver programs
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)route_runoff.f90

#========================================================================
# Compile the puppy
#========================================================================
# tasks
all: compile install clean

# compile
compile:
	$(FC_EXE) $(FLAGS) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(EXE)
	@echo "Succesfully compiled"

# Remove object files
clean:
						rm -f *.o
						rm -f *.lst
						rm -f *.mod
						rm -f *__genmod.f90
						@echo "Successfully cleaned object files"
#
# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(EXE) $(EXE_PATH)
	@echo "Successfully installed"
