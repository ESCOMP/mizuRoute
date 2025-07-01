.. _Build_model:

Build model
============

This version is used as stand-alone routing model, and also coupled with NCAR Community Terrestrial System Model (CTSM). This page is for building for stand-alone use

- A copy of the mizuRoute source code from `this repo <https://github.com/ESCOMP/mizuRoute>`_. You have a few options:

  - If you just want to use the latest release of mizuRoute, then simply look for the `latest release <https://github.com/ESCOMP/mizuRoute/releases>`_;
  - If you want the latest and greatest (and potentially erroneous), download a copy of the `cesm-coupling branch <https://github.com/ESCOMP/mizuRoute/tree/cesm-coupling>`_ (or clone it);

Dependencies
------------------------------------------

To compile mizuRoute, you will need:

- **Fortran compiler**: We recommend using the intel Fortran compiler or the gcc compiler, the latter of which is freely available. Since mizuRoute does not use any compiler-specific extensions, you should be able to compile mizuRoute with other Fortran compilers.

..

- **MPI (message passing interface) library**: `OpenMPI <https://www.open-mpi.org/>`_ is freely available and has been tested in mizuRoute. 

..

- **NetCDF libraries**: `NetCDF <http://www.unidata.ucar.edu/software/netcdf/>`_ or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. 

  Most \*nix package managers include a NetCDF port. Note that you need to ensure that:

  - You have NetCDF version 4.x;
  - The NetCDF libraries are compiled with the same compiler as you plan to use for compiling mizuRoute; and
  - You will have to have both NetCDF Fortran library installed (``libnetcdff.*``) and the C-version (``libnetcdf.*``).

- **PnetCDF libraries (optional)**: `PnetCDF <https://parallel-netcdf.github.io/>`_ is a parallel I/O library for accessing Unidata's NetCDF library.

..

- **CMake**: 


For Mac or Linux user, you may consider using `Homebrew <https://brew.sh/>`_ free and open-source software package management to obration all the libraries. The series of commands can be as simple as:

::

    brew install gcc
    brew install open-mpi
    brew install netcdf
    brew install pnetcdf
    brew install cmake

For windows user, .....

For HPC or larger cluster, please consult with system administrator(s) since the libaries are typically managed by them. 


Compilation
------------------------------------------

First, understand how mizuRoute directory is structured. There are hidden directories and files, but users mostly do not have to bother them.

     mizuRoute directory structure::

         mizuRoute/
         ├── route/
         |    ├── build/
         |    ├── bin/
         |    ├── input/
         |    ├── settings/
         |    └── ancillary_data/
         ├── bin/
         ├── cime_config/
         ├── docs/
         ├── libraries/
         ├── netcdf_test/
         ├── LICENSE
         ├── readme.md
         └── EADME_EXTERNALS.rst


Once you have all the external libraries, you can try compiling mizuRoute using the following steps for using the ``Makefile`` under ``mizuRoute/route/build``.
Here assume that all the libaries are installed using homebrew 🍺, and gnu compiler is used.

1. Navigate to your local copy of the mizuRoute top directory. 

..

2. Obtain `ParallelIO <https://github.com/NCAR/ParallelIO>`_ through git-fleximod tool that is already installed under ``mizuRoute/bin``. 

     Execute::

         mk -p libraries/parallelio
         ./bin/git-fleximod -g .gitmodules update

     See mizuRoute/README_EXTERNALS.rst for more details. 

3. Go to the ``route/build`` subdirectory.

..

4. Specify a number of environment variables that are used by the build process. 
   If you are using the ``bash`` shell, then you would set these environment variables with ``export``, e.g., ``export FC=gnu``.
   You will need to set the following:

   - ``BLDDIR``: This is the parent directory of the ``build`` directory.

     Then you would set by executing::

         export BLDDIR=`pwd`/../

   - ``FC``: Define the compiler family. This is only used to determine the compiler flags.

     if you use gnu compiler::

         export FC=gnu

     if you use intel compiler::

         export FC=intel

   - ``FC_EXE``: This is the actual compiler executable that is invoked. For most of compiler family, it is likely ``mpif90`` regardless of which compiler families.

     Example::

         export FC_EXE=mpif90

   - ``NCDF_PATH``: This is the path to the top level of NetCDF library directory. The directory typically contains ``bin include lib`` subdirectories. 

     Example (if netCDF is installed with homebrew)::

         export NCDF_PATH=/opt/homebrew/

   - ``PNETCDF_PATH`` (optional): This is also the path to top level of the PnetCDF directory. 

     Example (if pnetcdf is intalled with homebrew)::

         export PNETCDF_PATH=/opt/homebrew/


5. Once you have set up the environmental variables above, use the following command.

     ::
     
         make FC=$FC FC_EXE=$FC_EXE F_MASTER=$BLDDIR NCDF_PATH=$NCDF_PATH PNETCDF_PATH=$PNETCDF_PATH EXE=route_runoff 

     If the code compiles successfully, then the last line of output from the make process will tell you where the mizuRoute executable is installed (it goes into ``mizuRoute/route/bin``). 


NOTE:

   - You may add the variables directly in the ``Makefile``, rather than setting them as environment variables. They are located under ``User configure part``. 
     if you do that, you will just execute ``make`` (make sure to define ``EXE=<mizuRoute executable name>``)

..

   - To find netCDF and pnetCDF pathes, the following command might help.

     ::

         find / -type f \( -name "libnetcdf*.a*" \) -print

   - Often, netCDF-fortran and netCDF (c-version) libraries are located in separate location. If so, set variables ``NCDF_FORTRAN_PATH`` and ``NCDF_C_PATH``

     ::

        export NCDF_FORTRAN_PATH=<path_to_netcdf-fortran>
        export NCDF_C_PATH=<path_to_netcdf>
        make FC=$FC FC_EXE=$FC_EXE F_MASTER=$BLDDIR NCDF_C_PATH=$NCDF_C_PATH NCDF_FORTRAN_PATH=$NCDF_FORTRAN_PATH PNETCDF_PATH=$PNETCDF_PATH EXE=route_runoff
