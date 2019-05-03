[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.56043.svg)](http://dx.doi.org/10.5281/zenodo.56043)
# mizuRoute
mizuRoute is a stand-alone, post-processor of runoff output from a hydrologic model to produce streamflow estimates in the river network. The tool was developed for the large scale river routing (e.g., river network over contiguous United States), but works for small-size headwater basin as well.

# Contents in repository
This repository contains the source codes (Fortran90) and pre-process scripts (python and bash scripts) for mizuRoute. The repository also includes the user-manual, example dataset for a user to test mizuRoute, and netcdf test code (to make sure thtat netcdf library is loaded in the machine correctly).

The Fortran 90 source code of mizuRoute consists of two parts - 1) river network preprocessor and 2) routing program. The river network preprocessor, stored in ntopo, is to augment the basic information on river segment connectivity derived from GIS to facilitate subsequent execution of river routing computation. The routing program, stored in route, executes hillslope and river routing computation.

# To get started
1. Obtaining mizuRoute package. Just to use the tool, download the package by clicking “Download Zip” button on right column. 

2. Fortran compiler. We have successfully used the intel Fortran compiler (ifort), the GNU Fortran compiler (gfortran, version 4.8 or higher), and PGI fortran compiler (pgf90), the latter two of which are freely available. Since we do not use any compiler-specific extensions, mizuRoute should be complied with any Fortran compilers. If the user does not have a Fortran compiler, [gfortran](https://gcc.gnu.org/wiki/GFortran) can be installed for free. The easiest way is to use a package manager. Which package manager depends on your machine flavor. 

3. NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format, is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. All the mizuRoute I/Ol use NetCDF. The user needs to ensure that:
NetCDF version 4.x is installed in the linux-like machine.
NetCDF Fortran library is installed (libnetcdff.*) and not just the C-version.
The NetCDF libraries are compiled with the same compiler as you plan to compile mizuRoute
The user may use netcdf test code to check if NetCDF libraries are properly installed, 

4. Compiling the source code (rive network preprocessor and routing program). Once you have all the above, you can compile mizuRoute source codes using the following steps: Navigate to your local copy of the mizuRoute directory and go to the build subdirectory. The user will have to compile river network preprocessor and routing program separately.
 
    1. Edit F_MASTER and FC (to your desired compiler). You may also need to set NCDF_PATH and you may need to add some extra entries if you are using a different Fortran compiler or your setup is different (if someone wants to contribute an actual configure script that would be great).

    2. Type make under directory where Makefile is located. If all goes well, this will create  the executable runoff_route.exe (or process_river_topology.exe) to the bin directory. You may get some warnings (depending on your compiler settings), but you should not get any errors.

    3. Pay attention to the make output. You may need to set some environment variables (LD_LIBRARY_PATH in particular) to support dynamic linking;

    4. Ready to run the executables.

If you get this far then mizuRoute is installed correctly and functional. Now, the user will have to process runoff data and network topology data. Please refer to [user manual](docs/GMD_routing_v1_user_manual_20150831.pdf) to learn more about how to create mizuRoute input data  for your application. 

The user are encouraged to start with example data to get familiarize the process.
For real application, getting river network data (netCDF) might be time consuming part because this most likely requires GIS process and convert shapefile to netCDF). Please refer to section 2.1 what information is required. 
