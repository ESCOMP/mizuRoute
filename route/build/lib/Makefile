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
CIMEROOT     = $(F_MASTER)../cime
MCTDIR       = $(CIMEROOT)/src/externals/mct
MPISERIALDIR = $(MCTDIR)/mpi-serial
PIO2DIR      = $(CIMEROOT)/src/externals/pio2
LIBDIR       = $(F_MASTER)build/lib
PIOLIBDIR    = $(F_MASTER)build/lib/piolib
PIOLIB       = $(PIOLIBDIR)/lib/libpiof.a $(PIOLIBDIR)/lib/libpioc.a
PIOLIBMAKE   = $(PIOLIBDIR)/Makefile
MPISERLIBDIR = $(F_MASTER)build/lib/mpi-seriallib
MPISERLIB    = $(MPISERLIBDIR)/libmpi-serial.a
MPISERCONF   = $(MPISERLIBDIR)/Makefile.conf
MPIINC_PATHLIST  = $(CPATH)
INCLDIR      = -I.
CMAKE_OPTS   = -DCMAKE_Fortran_COMPILER=mpif90

all: $(CIMEROOT) $(MCTDIR) $(MPISERIALDIR) $(PIO2DIR)

$(CIMEROOT):
	@echo "Run manage_externals/checkout_externals in the root directory (cd ../../..) first in order to create the cime root directory"
	@echo "CIMEROOT being looked for: $(CIMEROOT)"

$(MCTDIR):
	@echo "MCTDIR does not exist: $(MCTDIR)"

$(MPISERIALDIR):
	@echo "MPISERIALDIR does not exist: $(MPISERIALDIR)"

$(PIO2DIR):
	@echo "PIO2DIR does not exist: $(PIO2DIR)"

ifeq "$(isPIO)" "yes"
  INCLDIR += -I$(PIOLIBDIR)/include
endif

ifeq "$(isMPISERIAL)" "yes"
  MPISERIAL = $(MPISERLIB)
  CPPDEFS += -DNO_MPI2
  INCLDIR += -I$(MPISERIALDIR) -I$(MPISERLIBDIR)
  ifdef PNETCDF_PATH
    undefine PNETCDF_PATH
  endif
  CMAKE_OPTS += -DMPI_C_INCLUDE_PATH=$(PIOLIBDIR)/include \
      -DMPI_Fortran_INCLUDE_PATH=$(PIOLIBDIR)/include \
      -DMPI_C_LIBRARIES=$(MPISERLIB) \
      -DMPI_Fortran_LIBRARIES=$(MPISERLIB)
else
  # Change colon delimited list to space with a -I in front of each
  MPIINC   = $(subst :, -I,$(MPIINC_PATHLIST))
  INCLDIR += -I$(MPIINC)
  CPPDEFS += -DHAVE_MPI
  MPISERIAL = 
endif

ifeq "$(isGPTL)" "yes"
  PIO_ENABLE_TIMING = ON
  PIOLIB       += $(PIOLIBDIR)/lib/libgptl.a
else
  PIO_ENABLE_TIMING = OFF
endif

CMAKE_OPTS += -D CMAKE_Fortran_FLAGS:STRING="$(FLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_C_FLAGS:STRING="$(CFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS) $(CPPDEFS) $(INCLDIR)" \
              -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
              -D PIO_ENABLE_TESTS:BOOL=OFF \
              -D CMAKE_INSTALL_PREFIX:PATH=$(PIOLIBDIR) \
              -D PIO_ENABLE_TIMING:BOOL=$(PIO_ENABLE_TIMING) \
              -D USER_CMAKE_MODULE_PATH:LIST="$(CIMEROOT)/src/CMake;$(CIMEROOT)/src/externals/pio2/cmake"
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
CMAKE_OPTS += -D GENF90_PATH=$(CIMEROOT)/src/externals/genf90

$(MPISERCONF):
	cp -p $(MPISERIALDIR)/*.h $(MPISERLIBDIR)/.
	cd $(MPISERLIBDIR) && $(MPISERIALDIR)/configure FC=$(FC) MPIFC=$(FC_EXE) FCFLAGS="$(FLAGS) \
                                              -I$(MPISERIALDIR)" MCT_PATH=$(MPISERIALDIR) --srcdir $(MPISERIALDIR)
	ln -sf $(MPISERIALDIR)/Makefile $(MPISERLIBDIR)/Makefile

$(MPISERLIB): $(MPISERCONF)
	$(MAKE) $(MFLAGS) -C $(MPISERLIBDIR) SRCDIR=$(MCTDIR)

cleanlibs: cleanmpiseriallib cleanpiolib
	@echo "Clean the libraries"
.PHONY : cleanlibs

cleanmpiseriallib:
	rm $(MPISERLIBDIR)/*.h $(MPISERLIB) $(MPISERCONF) $(MPISERLIBDIR)/Makefile
.PHONY : cleanmpiseriallib

cleanpiolib:
	rm -rf $(PIOLIBDIR)/CMake* $(PIOLIBDIR)/*.cmake $(PIOLIBDIR)/*.h $(PIOLIB) $(PIOLIBDIR)/include/* $(PIOLIBMAKE)
.PHONY : cleanpiolib

$(PIOLIBMAKE): $(MPISERIAL)
	cd $(PIOLIBDIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) $(PIO2DIR)

$(PIOLIB): $(PIOLIBMAKE)
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) piof
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) pioc
	$(MAKE) $(MFLAGS) -C $(PIOLIBDIR) install

db_extlibs:
	@echo "CIMEROOT:  "       $(CIMEROOT)
	@echo "libdir:    "       $(LIBDIR)
	@echo "CPP-DEFS:  "       $(CPPDEFS)
	@echo "mpi-serial:"       $(MPISERIALDIR)
	@echo "pio2dir:   "       $(PIO2DIR)
	@echo "mpi-seriallibdir:" $(MPISERLIBDIR)
	@echo "piolibdir: "       $(PIOLIBDIR)
	@echo "mpi-seriallib:   " $(MPISERLIB)
	@echo "mpi-serialconf:  " $(MPISERCONF)
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

