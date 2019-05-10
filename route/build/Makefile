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
isMPISERIAL =
isPIO =
 
PNETCDF_PATH =

# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER =


# External libraries (if used)
EXTLIBS =
EXTINCLUDES =

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
PIO_FILESYSTEM_HINTS = gpfs
CPPDEFS      = -DPIO1 -D_NO_MPI_RSEND
CIMEROOT     = $(F_MASTER)../cime
MCTDIR       = $(CIMEROOT)/src/externals/mct
MPISERIALDIR = $(MCTDIR)/mpi-serial
PIO1DIR      = $(CIMEROOT)/src/externals/pio1
LIBDIR       = $(F_MASTER)build/lib
PIOLIBDIR    = $(F_MASTER)build/lib/piolib
PIOLIB       = $(PIOLIBDIR)/pio/libpio.a
PIOLIBMAKE   = $(PIOLIBDIR)/Makefile
MPISERLIBDIR = $(F_MASTER)build/lib/mpi-seriallib
MPISERLIB    = $(MPISERLIBDIR)/libmpi-serial.a
MPISERCONF   = $(MPISERLIBDIR)/Makefile.conf
INCLDIR      = -I.
CMAKE_OPTS   = 

ifeq "$(isMPISERIAL)" "yes"
  EXTINCLUDES += -I$(MPISERLIBDIR)
  EXTLIBS += $(MPISERLIB)
  MPISERIAL = $(MPISERLIB)
  LDFLAGS += -L$(MPISERLIBDIR) -lmpi-serial
  CPPDEFS += -DNO_MPI2
  INCLDIR += -I$(MPISERLIBDIR)
  ifdef PNETCDF_PATH
    undefine PNETCDF_PATH
  endif
  CMAKE_OPTS += -DMPI_C_INCLUDE_PATH=$(PIOLIBDIR) \
      -DMPI_Fortran_INCLUDE_PATH=$(PIOLIBDIR) \
      -DMPI_C_LIBRARIES=$(MPISERLIB) \
      -DMPI_Fortran_LIBRARIES=$(MPISERLIB)
else
  INCLDIR += -I$(MPI_ROOT)/include
  CPPDEFS += -DHAVE_MPI
  MPISERIAL = 
endif

ifeq "$(isPIO)" "yes"
  EXTLIBS += $(PIOLIB)
  LDFLAGS += -L$(PIOLIBDIR)/pio -lpio
endif

CMAKE_OPTS += -D CMAKE_Fortran_FLAGS:STRING="$(FLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_C_FLAGS:STRING="$(CFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
              -D PIO_ENABLE_TESTS:BOOL=OFF \
              -D USER_CMAKE_MODULE_PATH:LIST="$(CIMEROOT)/src/CMake;$(CIMEROOT)/src/externals/pio2/cmake"
CMAKE_OPTS += -D NetCDF_PATH:PATH=$(NCDF_PATH)




ifdef PNETCDF_PATH
        CMAKE_OPTS += -D PnetCDF_PATH:STRING="$(PNETCDF_PATH)"
else
        CMAKE_OPTS += -D WITH_PNETCDF:LOGICAL=FALSE -D PIO_USE_MPIIO:LOGICAL=FALSE
endif
ifdef PIO_FILESYSTEM_HINTS
        CMAKE_OPTS += -D PIO_FILESYSTEM_HINTS:STRING="$(PIO_FILESYSTEM_HINTS)"
endif
CMAKE_OPTS += -D GENF90_PATH=$(CIMEROOT)/src/externals/genf90

$(MPISERCONF):
	cp -p $(MPISERIALDIR)/*.h $(MPISERLIBDIR)/.
	cd $(MPISERLIBDIR) && $(MPISERIALDIR)/configure FC=$(FC) MPIFC=$(FC_EXE) FCFLAGS="$(FLAGS) \
                                              -I$(MPISERLIBDIR)" MCT_PATH=$(MPISERIALDIR) --srcdir $(MPISERIALDIR)
	ln -sf $(MPISERIALDIR)/Makefile $(MPISERLIBDIR)/Makefile

$(MPISERLIB): $(MPISERCONF)
	$(MAKE) $(MFLAGS) -C $(MPISERLIBDIR) SRCDIR=$(MCTDIR)

cleanmpiseriallib:
	$(MAKE) $(MFLAGS) -C $(MPISERLIBDIR) SRCDIR=$(MCTDIR) clean
	rm $(MPISERCONF)

cleanpiolib:
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) clean
	rm $(PIOLIBDIR)/CMakeCache.txt
	rm $(PIOLIBMAKE)

$(PIOLIBMAKE): $(MPISERIAL)
	cd $(PIOLIBDIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) $(PIO1DIR)

$(PIOLIB): $(PIOLIBMAKE)
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR)

extlibs:
	@echo "CIMEROOT:  "       $(CIMEROOT)
	@echo "libdir:    "       $(LIBDIR)
	@echo "CPP-DEFS:  "       $(CPPDEFS)
	@echo "mpi-serial:"       $(MPISERIALDIR)
	@echo "pio1dir:   "       $(PIO1DIR)
	@echo "mpi-seriallibdir:" $(MPISERLIBDIR)
	@echo "piolibdir: "       $(PIOLIBDIR)
	@echo "mpi-seriallib:   " $(MPISERLIB)
	@echo "mpi-serialconf:  " $(MPISERCONF)
	@echo "piolib:    "       $(PIOLIB)
	@echo "piolib:Make"       $(PIOLIBMAKE)
	@echo "AR:          "     $(AR)
	@echo "MAKE:        "     $(MAKE)
	@echo "MFLAGS:      "     $(MFLAGS)
	@echo "FLAGS:       "     $(FLAGS)
	@echo "LDFLAGS:     "     $(LDFLAGS)
	@echo "CMAKE_OPTS:  "     $(CMAKE_OPTS)
	@echo "CMAKE_ENV_VARS:"   $(CMAKE_ENV_VARS)



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
