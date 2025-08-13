[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.595402.svg)](https://doi.org/10.5281/zenodo.595402)
[![Documentation Status](https://readthedocs.org/projects/mizuroute/badge/?version=main)](https://mizuroute.readthedocs.io/en/main/?badge=main)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

# mizuRoute
MizuRoute is a tool that takes runoff data from a hydrologic or land surface model to compute streamflow in a user-defined catchment-based (vector) river network.
A user can add Lakes into the river network, allowing the model to simulate discharge and volume in rivers and lakes (natural lakes or reservoir).

MizuRoute uses Message Passing Interface (MPI) for parallel computing to work on multi-decdal to a century simulation over large river networks (e.g., the river network over the contiguous United States).

MizuRoute can be run in two ways: 

1. Stand-alone mode - users needs provides their own runoff time series as a model input. 
2. CESM coupled model - MizuRoute runs as a river model component in Community Earth System Model (CESM).
Currently, mizuRoute is coupled to Community Terrestrial Systems Model ([CTSM](https://github.com/ESCOMP/CTSM)), the land-model component of CESM. 


# To get started for a stand-alone mode
1. Obtaining mizuRoute package. Just to use the tool, download the package by clicking “Download Zip” button on right column. 

2. To compile the codes, and prepare for the input data, please refer to [User's Guide](https://mizuroute.readthedocs.io/en/main/)

3. A user is encouraged to start with example data to get familiarize the process. Link to testCase data are given in [testCase data](https://mizuroute.readthedocs.io/en/main/users_guide/testCase.html) in User's Guide.

# To get started for ctsm coupling mode

User interested in using mizuRoute with CTSM is referred to CESM's user guide. Here, quick guide is provided.  

1. Obtain CTSM code from [github](https://github.com/ESCOMP/CTSM/tree/master)

2. Create the case

   ```bash
   cd cime/scripts
   ```

   ```bash
   ./create_newcase --case <testcase> --mach derecho --res f09_f09_rHDMAlk_mg17 -compset I2000Clm60SpMizGs 
   ```
   (`./create_newcase -help --` to get help on the script)

   # Setup the case

   ```bash
   cd <testcase>
   ```

   ```bash
   ./xmlchange id1=val1,id2=val2  # to make changes to any settings in the env_*.xml files
   ./case.setup
   ```
   (./case.setup -help -- to get help on the script)

   # Add any namelist changes to the `user_nl_*` files

   ```bash
   $EDITOR user_nl_*
   ```

   # Compile the code

   ```bash
   ./case.build
   ```

   # Submit the run

   ```bash
   ./case.submit
   ```

# Prerequisite 
 1. Linux commands.
 2. Geographic Information System to develop and visualize river network data
 3. python or similar other languages to analyze/visualize data, prepare for the input
 4. Fortran (if a user desires to change the codes)


LICENSE
-------

mizuRoute is [licensed](LICENSE.txt) under Apache 2.0
