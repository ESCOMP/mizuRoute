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

# Define the compiler exe
# For MPI use, mpifort, mpif77, mpif90, mpif08
FC_EXE =

# Define the executable
EXE =

# Define optional setting
# fast:      Enables optimizations
# debug:     Minimum debug options, still
# profile:   Enables profiling
MODE = debug

# define open MP flags (yes if used)
isOpenMP = yes
FLAGS_OMP =
LIBOPENMP =

# define if should use external librarys (yes if used)
# external librarys are: mpi-serial, Parallel-I/O, and GPTL timing library
# define path to PNETCDF if you want to use it
isMPISERIAL =
isPIO = yes
isGPTL = yes

# If isPIO==yes then you can give hints to PIO about the local system
# gpfs is the specific filesystem used on cheyenne
# (can be gpfs or lustre) Only use it if you know for sure about your filesystem.
# If not make this setting blank
PIO_FILESYSTEM_HINTS = gpfs

# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER =

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gnu"
 NCDF_PATH =
 PNETCDF_PATH =
endif
ifeq "$(FC)" "intel"
 NCDF_PATH =
 PNETCDF_PATH =
endif
ifeq "$(FC)" "pgi"
 NCDF_PATH =
 PNETCDF_PATH =
endif

LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
            -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

ifeq "$(isOpenMP)" "yes"
 ifeq "$(FC)" "pgi"
  FLAGS_OMP = -mp
 else
  FLAGS_OMP = -fopenmp
 endif
endif

#========================================================================
# Define flags
#========================================================================
ifeq "$(FC)" "gnu"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -fmax-errors=0 -ffree-line-length-none $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -Wall -fmax-errors=0 -fbacktrace -fbounds-check -ffpe-trap=zero -ffree-line-length-none $(FLAGS_OMP)
  endif
endif

ifeq "$(FC)" "intel"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -FR -auto -fltconsistency -fpe0 $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -debug all -warn all -check all -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -debug -warn all -check all -noerror_limit -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency $(FLAGS_OMP)
  endif
endif

ifeq "$(FC)" "pgi"
  ifeq "$(MODE)" "fast"
    FLAGS = -fast -O3 -Mdclchk $(FLAGS_OMP)
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -traceback -Mbounds -Mlist -Minfo -Mdclchk $(FLAGS_OMP)
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

# External libraries (if used)
EXTLIBS =
EXTINCLUDES =

#========================================================================
# Assemble all of the sub-routines
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
    popMetadat.f90 \
    allocation.f90
# define utilities
UTILS = \
    nr_utility.f90 \
    ascii_util.f90 \
    time_utils.f90 \
    ncio_utils.f90 \
    mpi_utils.f90 \
    pio_utils.f90 \
    model_utils.f90 \
    gamma_func.f90
# read/write files
IO = \
    remap.f90 \
    read_control.f90 \
    read_param.f90 \
    read_streamSeg.f90 \
    write_streamSeg.f90 \
    standalone/read_runoff.f90 \
    standalone/get_basin_runoff.f90 \
    standalone/read_remap.f90 \
    read_restart.f90 \
    write_restart_pio.f90 \
    write_simoutput_pio.f90
# initialization
INIT = \
    process_time.f90 \
    network_topo.f90 \
    process_param.f90 \
    process_ntopo.f90 \
    pfafstetter.f90 \
    domain_decomposition.f90
# CORE
CORE = \
    accum_runoff.f90 \
    basinUH.f90 \
    irf_route.f90 \
    kwt_route.f90 \
    kwe_route.f90 \
    main_route.f90 \
    mpi_process.f90 \
    init_model_data.f90 \
    standalone/model_setup.f90

# concatanate model subroutines
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(INIT) $(CORE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# Define driver programs
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)standalone/route_runoff.f90

# tasks
all: compile install clean

#========================================================================
# External libaries that might need to be built
#========================================================================
LIBDIR       = $(F_MASTER)build/lib
PIOLIBDIR    = $(LIBDIR)/piolib/lib
PIOINCDIR    = $(LIBDIR)/piolib/include
PIOLIB       = $(PIOLIBDIR)/libpiof.a $(PIOLIBDIR)/libpioc.a
MPISERLIBDIR = $(LIBDIR)/mpi-seriallib
MPISERLIB    = $(MPISERLIBDIR)/libmpi-serial.a

ifeq "$(isMPISERIAL)" "yes"
  EXTINCLUDES += -I$(MPISERLIBDIR)
  EXTLIBS += $(MPISERLIB)
  LDFLAGS += -L$(MPISERLIBDIR) -lmpi-serial
  ifdef PNETCDF_PATH
    undefine PNETCDF_PATH
  endif
  MPISERIAL = $(MPISERLIB)
else
  MPISERIAL =
endif

ifeq "$(isPIO)" "yes"
  EXTINCLUDES += -I$(PIOINCDIR)
  EXTLIBS += $(PIOLIB)
  LDFLAGS += -L$(PIOLIBDIR) -lpiof -lpioc

  ifeq "$(isGPTL)" "yes"
     LDFLAGS += -lgptl
  endif
endif

ifdef PNETCDF_PATH
        LDFLAGS    += -L$(PNETCDF_PATH)/lib -lpnetcdf
endif

$(PIOLIB): $(MPISERIAL)
	cd $(LIBDIR); \
	$(MAKE) $(MFLAGS) F_MASTER=$(F_MASTER) FC=$(FC) FC_EXE=$(FC_EXE) FLAGS="$(FLAGS)" \
	PIO_FILESYSTEM_HINTS=$(PIO_FILESYSTEM_HINTS) PNETCDF_PATH=$(PNETCDF_PATH) isPIO=yes isGPTL=$(isGPTL) \
	MPIINC_PATHLIST=$(CPATH) NCDF_PATH=$(NCDF_PATH) $(PIOLIB)

$(MPISERLIB):
	cd $(LIBDIR); \
	$(MAKE) $(MFLAGS) F_MASTER=$(F_MASTER) FC=$(FC) FC_EXE=$(FC_EXE) FLAGS="$(FLAGS)" isMPISERIAL=yes $(MPISERLIB)

# Clean external libs
cleanlibs: cleanpiolib cleanmpiseriallib
	@echo "Clean the pio and mpi-serial external libraries"
.PHONY : cleanlibs

cleanpiolib:
	cd $(LIBDIR); \
	$(MAKE) $(MFLAGS) F_MASTER=$(F_MASTER) cleanpiolib
.PHONY : cleanpiolib

cleanmpiseriallib:
	cd $(LIBDIR); \
	$(MAKE) $(MFLAGS) F_MASTER=$(F_MASTER) cleanmpiseriallib
.PHONY : cleanmpiseriallib

#========================================================================
# Compile the puppy
#========================================================================
# compile
compile: $(EXTLIBS)
	$(FC_EXE) $(FLAGS) $(EXTINCLUDES) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) $(LDFLAGS) -o $(EXE)

# Remove object files
clean:
						rm -f *.o
						rm -f *.lst
						rm -f *.mod
						rm -f *__genmod.f90
#
# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(EXE) $(EXE_PATH)
