.. _Build_model:

Build model
============

Dependencies
------------------------------------------

To compile mizuRoute, you will need:

- a Fortran compiler. We recommend using the intel Fortran compiler or the GNU Fortran compiler, the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile mizuRoute with other Fortran compilers such as NAG as well.

  If you do not have a Fortran compiler, you can install ``gfortran`` for free. The easiest way is to use a package manager. Note that ``gfortran`` is installed as part of the ``gcc`` compiler suite.

- a mpi (message passing interface) library. OpenMPI is freely available and has been tested in mizuRoute. 

- the NetCDF libraries. `NetCDF <http://www.unidata.ucar.edu/software/netcdf/>`_ or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. They are widely used in the hydrometeorological community and eventually almost all SUMMA I/O will use NetCDF. Most \*nix package managers include a NetCDF port. Note that you need to ensure that:

  - You have NetCDF version 4.x;
  - The NetCDF libraries are compiled with the same compiler as you plan to use for compiling mizuRoute; and
  - You will have to have both NetCDF Fortran library installed (``libnetcdff.*``) and the C-version (``libnetcdf.*``).

- a copy of the mizuRoute source code from `this repo <https://github.com/ESCOMP/mizuRoute>`_. You have a number of options:

  - If you just want to use the latest stable release of SUMMA, then simply look for the `latest release <https://github.com/ESCOMP/mizuRoute/releases>`_;
  - If you want the latest and greatest (and potentially erroneous), download a copy of the `cesm-coupling branch <https://github.com/ESCOMP/mizuRoute/tree/cesm-coupling>`_ (or clone it);

Compilation
------------------------------------------

To compile mizuRoute, please follow the following steps. 



