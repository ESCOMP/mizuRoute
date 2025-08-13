[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.595402.svg)](https://doi.org/10.5281/zenodo.595402)
[![Documentation Status](https://readthedocs.org/projects/mizuroute/badge/?version=main)](https://mizuroute.readthedocs.io/en/main/?badge=main)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# mizuRoute
Stand-alone mizuRoute is a post-processor of runoff output from a hydrologic model or Land surface model to produce streamflow estimates in the river network provided by a user. 
The tool was developed with Message Passing Interface (MPI) for the large domain network based river routing (e.g., river network over contiguous United States or larger), but works for gridded river network as well.

# To get started
1. Obtaining mizuRoute package. Just to use the tool, download the package by clicking “Download Zip” button on right column. 

2. To compile the codes, and prepare for the input data, please refer to [User's Guide](https://mizuroute.readthedocs.io/en/main/)

A user is encouraged to start with example data to get familiarize the process. Link to testCase data are given in [testCase data](https://mizuroute.readthedocs.io/en/main/users_guide/testCase.html) in User's Guide.

# Prerequisite 
 1. Linux commands.
 2. Geographic Information System to develop and visualize river network data
 3. python or similar other languages to analyze/visualize data, prepare for the input
 4. Fortran (if a user desires to change the codes)


LICENSE
-------

mizuRoute is [licensed](LICENSE.txt) under Apache 2.0
