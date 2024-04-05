#========================================================================
# Makefile for the external libaries that might need to be built
# These are all under the top level "cime" directory.
#========================================================================

#=============================================
# First check that important settings are set
#=============================================

#======================================
# Library definitions
#=-====================================
CPPDEFS      = -DUSE_MPI_ISEND_FOR_FC
PIO2DIR      = $(F_MASTER)../libraries/parallelio
LIBDIR       = $(F_MASTER)build/lib
PIOLIBDIR    = $(F_MASTER)build/lib/piolib
PIOLIB       = $(PIOLIBDIR)/lib/libpiof.a $(PIOLIBDIR)/lib/libpioc.a
PIOLIBMAKE   = $(PIOLIBDIR)/Makefile
MPIINC_PATHLIST  = $(CPATH)
INCLDIR      = -I.
CMAKE_OPTS   =

all: $(PIO2DIR)

$(PIO2DIR):
	@echo "PIO2DIR does not exist: $(PIO2DIR)"

ifeq "$(isPIO)" "yes"
  INCLDIR += -I$(PIOLIBDIR)/include
  CPPDEFS += -DPIO2
endif

# Change colon delimited list to space with a -I in front of each
ifdef MPIINC_PATHLIST
  MPIINC   = $(subst :, -I,$(MPIINC_PATHLIST))
  INCLDIR += -I$(MPIINC)
endif
CPPDEFS += -DHAVE_MPI
CMAKE_OPTS += -DCMAKE_Fortran_COMPILER=mpif90

ifeq "$(isGPTL)" "yes"
  PIO_ENABLE_TIMING = ON
  PIOLIB       += $(PIOLIBDIR)/lib/libgptl.a
else
  PIO_ENABLE_TIMING = OFF
endif

CMAKE_OPTS += -Wno-dev -D CMAKE_Fortran_FLAGS:STRING="$(FLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_C_FLAGS:STRING="$(CFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
              -D PIO_ENABLE_TESTS:BOOL=OFF \
              -D CMAKE_INSTALL_PREFIX:PATH=$(PIOLIBDIR) \
              -D PIO_ENABLE_TIMING:BOOL=$(PIO_ENABLE_TIMING) \
              -D USER_CMAKE_MODULE_PATH:LIST="$(PIO2DIR)/cmake"
CMAKE_OPTS += -D NetCDF_PATH:PATH=$(NCDF_PATH)

ifdef PNETCDF_PATH
        CMAKE_OPTS += -D PnetCDF_PATH:STRING="$(PNETCDF_PATH)"
        CPPDEFS    += -D_PNETCDF
        INCLDIR    += -I$(PNETCDF_PATH)/include
else
        CMAKE_OPTS += -D WITH_PNETCDF:LOGICAL=FALSE -D PIO_USE_MPIIO:LOGICAL=FALSE
endif
ifdef PIO_FILESYSTEM_HINTS
        CMAKE_OPTS += -D PIO_FILESYSTEM_HINTS:STRING="$(PIO_FILESYSTEM_HINTS)"
endif
CMAKE_OPTS += -D GENF90_PATH=$(PIO2DIR)/scripts

MODEARGS :=

ifeq "$(MODE)" "debug"
   CMAKE_OPTS += -DCMAKE_BUILD_TYPE=Debug
endif

cleanlibs: cleanpiolib
	@echo "Clean the libraries"
.PHONY : cleanlibs

cleanpiolib:
	rm -rf $(PIOLIBDIR)/CMake* $(PIOLIBDIR)/*.cmake $(PIOLIBDIR)/*.h $(PIOLIB) $(PIOLIBDIR)/include/* $(PIOLIBMAKE)
.PHONY : cleanpiolib

$(PIOLIBMAKE):
	cd $(PIOLIBDIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) $(MODEARGS) $(PIO2DIR)

$(PIOLIB): $(PIOLIBMAKE)
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) piof
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) pioc
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) install

db_extlibs:
	@echo "libdir:    "       $(LIBDIR)
	@echo "CPP-DEFS:  "       $(CPPDEFS)
	@echo "pio2dir:   "       $(PIO2DIR)
	@echo "piolibdir: "       $(PIOLIBDIR)
	@echo "piolib:    "       $(PIOLIB)
	@echo "piolib:Make"       $(PIOLIBMAKE)
	@echo "MPIINC_PATHLIST: " $(MPIINC_PATHLIST)
	@echo "INCLDIR      "     $(INCLDIR)
	@echo "AR:          "     $(AR)
	@echo "MAKE:        "     $(MAKE)
	@echo "MFLAGS:      "     $(MFLAGS)
	@echo "FLAGS:       "     $(FLAGS)
	@echo "CMAKE_OPTS:  "     $(CMAKE_OPTS)
	@echo "CMAKE_ENV_VARS:"   $(CMAKE_ENV_VARS)
.PHONY : db_extlibs

