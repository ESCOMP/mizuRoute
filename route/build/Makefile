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

# define open MP flags (yes if used)
isOpenMP =
FLAGS_OMP =
LIBOPENMP =

# define if should use external librarys (yes if used)
# external librarys are: mpi-serial, Parallel-I/O, and GPTL timing library
# define path to PNETCDF if you want to use it
isMPISERIAL =
isPIO =
isGPTL = 
 
PNETCDF_PATH =

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
ifeq "$(FC)" "gfortran"
  NCDF_PATH = 
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
              -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
  ifeq "$(isOpenMP)" "yes"
   FLAGS_OMP = -fopenmp
  endif
endif

ifeq "$(FC)" "ifort"
  NCDF_PATH = 
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
		          -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
  ifeq "$(isOpenMP)" "yes"
   FLAGS_OMP = -fopenmp
  endif
endif

ifeq "$(FC)" "pgf90"
 NCDF_PATH = 
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
 INCNETCDF = -I$(NCDF_PATH)/include
endif

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS = -p -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=all $(FLAGS_OMP)
 FLAGS1 = -O3 -ffree-line-length-none -fmax-errors=0 $(FLAGS_OMP)
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

# External libraries (if used)
EXTLIBS =
EXTINCLUDES =

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
    popMetadat.f90 \
    allocation.f90 
#
# define utilities
UTILS = \
    nr_utility.f90 \
    ascii_util.f90 \
    time_utils.f90 \
    mpi_utils.f90 \
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
    get_basin_runoff.f90 \
    read_remap.f90 \
    read_restart.f90 \
    write_restart.f90 \
    write_simoutput.f90

INITIAL = \
    process_time.f90 \
    network_topo.f90 \
    process_param.f90 \
    process_ntopo.f90 \
    pfafstetter.f90 \
    domain_decomposition.f90

# perform model routing
QROUTE = \
    accum_runoff.f90 \
    basinUH.f90 \
    irf_route.f90 \
    kwt_route.f90 \
    main_route.f90 \
    mpi_process.f90 \
    model_setup.f90

# concatanate model subroutines
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(INITIAL) $(QROUTE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# PART 3: Define driver programs, and executables
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)route_runoff.f90

# Define the executable
EXE = 

# tasks
all: compile install clean

#========================================================================
# PART 3.5: External libaries that might need to be built
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
# PART 4: Compile the puppy
#========================================================================
# compile
compile: $(EXTLIBS)
	$(FC_EXE) $(FLAGS) $(EXTINCLUDES) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) $(LDFLAGS) -o $(EXE)

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
