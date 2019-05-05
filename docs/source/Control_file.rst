Control file
============

Control file is a simple text file, mainly defining model control such as simulation time, file name and locations, routing option etc. 
Variables in control file are read in the beginning of the code (see ``./build/src/read_control.f90``) and 
saved in fortran variable specified by tag (inside <> in table) and as public variables (see ``./build/src/public_var.f90``) . 
Some of such public varialbes have some default values, but most of them are not defined.
Those undefined variables need to be defined in control file.   
Other variables in supplement table have their default values but can be also included in control file to overwrite the values. 
The order of variables in the control file does not matter. However, grouping variables into similar themes is recommended for readibility. 

Minimum required variables depends on runoff input options

Example of control file is given in ``./route/settings/SAMPLE.control``

Some of rules

* Exclamation mark is for comment and after exclamation make is ignored for reading.
* Format: <tag>    variable    ! comments
* tag is Fortran variable name and cannot be changed and have to be enclosed by <>
* need ! after variable, otherwise getting error.


The following variables (not pre-defined in the code) need to be defined in control file.

+--------+------------------------+-------------------------------------------------------------------------------------------+
| option | tag                    | Description                                                                               |
+========+========================+===========================================================================================+
| 1,2,3  | <ancil_dir>            | Directory that contains ancillary data (river netowrk, remapping, and parameter namelist) |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <input_dir>            | Directory that contains runoff data                                                       |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <output_dir>           | Directory that contains runoff data                                                       |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <param_nml>            | Spatially constant parameter namelist (should be stored in <ancil_dir>                    |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <sim_start>            | time of simulation start. format: yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                       |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <sim_end>              | time of simulation end. format:  yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                        |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_ntopOld>        | name of input netCDF for River Network                                                    |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_sseg>           | dimension name for reach                                                                  |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_nhru>           | dimension name for RN_HRU                                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_qsim>           | netCDF name for HM_HRU runoff                                                             |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <vname_qsim>           | variable name for HM_HRU runoff                                                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <vname_time>           | variable name for time                                                                    |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 2      | <vname_hruid>          | variable name for HM_HRU ID                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 3      | <dname_xlon>           | dimension name for x, lat, or j dimension                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 3      | <dname_ylat>           | dimension name for y, lon, or i dimension                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_time>           | dimension name for time                                                                   |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_hruid>          | dimension name for HM_HRU                                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <units_qsim>           | units of runoff. e.g., mm/s                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dt_qsim>              | time interval of runoff data in second. e.g., 86400 sec for daily step                    |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <is_remap>             | Logical to indicate runoff needs to be remapped to RN_HRU. T or F                         |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <fname_remap>          | netCDF name of runoff remapping                                                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <vname_hruid_in_remap> | variable name for RN_HRUs                                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <vname_weight>         | variable name for areal weights of overlapping HM_HRUs                                    |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2    | <vname_qhruid>         | variable name for HM_HRU ID                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|     3  | <vname_i_index>        | variable name of ylat index (runoff input option 3)                                       |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|     3  | <vname_j_index>        | variable name of xlon index (runoff input option 3)                                       |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <vname_num_qhru>       | variable name for a numbers of overlapping HM_HRUs with RN_HRUs                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <dname_hru_remap>      | dimension name for HM_HRU (option 2)                                                      |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <dname_data_remap>     | dimension name for data                                                                   |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <restart_opt>          | option to use saved restart file T->yes, F->No (start with empty channels)                |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <route_opt>            | option for routing schemes 0-> both, 1->IRF, 2->KWT, otherwise error                      |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_output>         | output netCDF name for model simulation results                                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_state_in>       | input restart netCDF name.                                                                | 
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_state_out>      | output restart netCDF name.                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+

Variables that have default values but can be overwritten 

+------------------------+------------------------+--------------------------------------------------------------------------+
| tag                    | Default values         | Description                                                              |
+========================+========================+==========================================================================+
| <ntopWriteOption>      | F                      | logical to indicate if augmented river network is written. T or F        |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <seg_outlet>           | -999                   | outlet reach ID for subsetted river basin. See note 1                    |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <fname_ntopNew>        | <fname_ntopOld>_aug.nc | output netCDF name for augmented river network. See note 1               |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <newFileFrequency>     | annual                 | frequency for new output files (day, month, or annual)                   |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <hydGeometryOption>    | 1                      | option for hydraulic geometry calculations (0=read from file, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <topoNetworkOption>    | 1                      | option for network topology calculations (0=read from file, 1=compute)   |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <computeReachList>     | 1                      | option to compute list of upstream reaches (0=do not compute, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+

1. if <ntopWriteOption> is T or <seg_outlet> is not -999, river network is output into <fname_ntopNew> and program stops 



Examples

Option 1 - HM_HRU runoff input::

  ! ****************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ******************************************************************************
  ! ****************************************************************************************************************************
  ! ****************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!" 
  !
  ! *************************************************************************************************************************
  ! PART 1: DEFINE DIRECTORIES 
  ! --------------------------
  <ancil_dir>    ./ancillary_data/                    ! directory containing ancillary data (river network, remapping netCDF)
  <input_dir>    ./input/                             ! directory containing input data (runoff netCDF)
  <output_dir>   ./output/                            ! directory containing output data
  ! *************************************************************************************************************************
  ! PART 2: DEFINE TIME PERIOD OF THE SIMULATION
  ! --------------------------------------------
  <sim_start>         1950-1-1                        ! time of simulation start (year-month-day)
  <sim_end>           1951-1-1                        ! time of simulation end (year-month-day)
  ! **************************************************************************************************************************
  ! PART 3: DEFINE FINE NAME AND DIMENSIONS
  ! ---------------------------------------
  <fname_ntopOld>    ntopo_entire.nc                  ! name of netCDF containing river segment data 
  <dname_sseg>       seg                              ! dimension name of the stream segments
  <dname_nhru>       hru                              ! dimension name of the HRUs
  ! **************************************************************************************************************************
  ! PART 4: DEFINE DESIRED VARIABLES FOR THE NETWORK TOPOLOGY
  ! ---------------------------------------------------------
  <seg_outlet>  -9999                                 ! reach ID of outlet streamflow segment. -9999 for all segments 
  ! **************************************************************************************************************************
  ! PART 5: DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>  runoff.HM_HRU.nc                      ! name of netCDF containing the HRU runoff
  <vname_qsim>  RUNOFF                                ! name of HRU runoff variable
  <vname_time>  time                                  ! name of time variable in the runoff file
  <vname_hruid> hru                                   ! name of runoff HRU id variable
  <dname_time>  time                                  ! name of time dimension 
  <dname_hruid> hru                                   ! name of the HRU dimension 
  <units_qsim>  mm/s                                  ! units of runoff
  <dt_qsim>     86400                                 ! time interval of the runoff
  ! **************************************************************************************************************************
  ! PART 6: DEFINE RUNOFF MAPPING FILE 
  ! ----------------------------------
  <is_remap>    F                                     ! logical to indicate runnoff needs to be mapped to river network HRU 
  ! **************************************************************************************************************************
  ! PART 7 DEFINE RUN CONTROL 
  ! ---------------------------
  <restart_opt> F                                     ! option to use restart file T->yes, F->No (start with empty channels) 
  <route_opt>   0                                     ! option for routing schemes 0-> both, 1->IRF, 2->KWT otherwise error 
  ! **************************************************************************************************************************
  ! PART 8: DEFINE OUTPUT FILE
  ! ---------------------------
  <fname_output>    flow_                             ! prefix of netCDF for the model output (netCDF name = flow_nomapping_yyyy.nc)
  <fname_state_in>  state.in.nc                       ! netCDF name for the model state input 
  <fname_state_out> state.out.nc                      ! netCDF name for the channel state output 
  ! **************************************************************************************************************************
  ! PART 10: Namelist file name 
  ! ---------------------------
  <param_nml>    param.nml.default                    ! spatially constant model parameters    
  ! **************************************************************************************************************************

