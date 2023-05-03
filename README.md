[![Documentation Status](https://readthedocs.org/projects/mizuroute/badge/?version=cesm-coupling)](https://mizuroute.readthedocs.io/en/cesm-coupling/?badge=cesm-coupling)

# mizuRoute
Stand-alone mizuRoute is a post-processor of runoff output from a hydrologic model or Land surface model to produce streamflow estimates in the river network provided by a user. The tool was developed for the large scale, network based river routing (e.g., river network over contiguous United States), but works for gridded river network as well.

Technical documentation is now being built on [readthedocs](https://mizuroute.readthedocs.io/en/cesm-coupling/)

# To get started
1. Obtaining mizuRoute package. Just to use the tool, download the package by clicking “Download Zip” button on right column. 

2. Fortran compiler. We have successfully used the intel Fortran compiler (ifort), the GNU Fortran compiler (gfortran, version 4.8 or higher), and PGI fortran compiler (pgf90), the latter two of which are freely available. Since we do not use any compiler-specific extensions, mizuRoute should be complied with any Fortran compilers. If the user does not have a Fortran compiler, [gfortran](https://gcc.gnu.org/wiki/GFortran) can be installed for free. The easiest way is to use a package manager. Which package manager depends on your machine flavor. 

3. NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format, is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. All the mizuRoute I/Ol use NetCDF. The user needs to ensure that:
NetCDF version 4.x is installed in the linux-like machine.
NetCDF Fortran library is installed (libnetcdff.*) and not just the C-version.
The NetCDF libraries are compiled with the same compiler as you plan to compile mizuRoute
The user may use netcdf test code to check if NetCDF libraries are properly installed, 

4. Compiling the source code (rive network preprocessor and routing program). Once you have all the above, you can compile mizuRoute source codes using the following steps: Navigate to your local copy of the mizuRoute directory and go to the build subdirectory.
 
    1. Edit F_MASTER and FC (to your desired compiler). You may also need to set NCDF_PATH and you may need to add some extra entries if you are using a different Fortran compiler or your setup is different (if someone wants to contribute an actual configure script that would be great).

    2. Type make under directory where Makefile is located. If all goes well, this will create  the executable runoff_route.exe (or process_river_topology.exe) to the bin directory. You may get some warnings (depending on your compiler settings), but you should not get any errors.

    3. Pay attention to the make output. You may need to set some environment variables (LD_LIBRARY_PATH in particular) to support dynamic linking;

    4. Ready to run the executables.

If you get this far then mizuRoute is built correctly and functional. Now, the user will have to generate input data, runoff data, river network topology and runoff mapping data (depending on input runoff option). Please refer to [readthedocs](https://mizuroute.readthedocs.io/en/develop/) to learn more about mizuRoute input data. 

The user are encouraged to start with example data to get familiarize the process. testCase are being now developed and posted separately.

For real application, getting river network data (netCDF) might be time consuming because this most likely requires GIS process and convert shapefile to netCDF). 
