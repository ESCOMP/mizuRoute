Control file
============

Control file is simple text file, mainly defining model control such as simulation time, file name and locations, routing option etc. 
Variables in control file are read in the beginning of the code and saved in fortran variable named tag and as public variables. 
Some of such public varialbes have some default values, but most of them are not defined.
Those undefined variables need to be defined in control file.   
Other variables in supplement table do have their default values but can be also included in control file to overwrite the values. 
The order of variables in the control file is flexible. However, grouping variables into similar theme is recommended for readibility. 

Example of control file  is given in ./route/settings/SAMPLE.control

Some of rules

* Exclamation mark is for comment and after exclamation make is ignored for reading.
* Format: <tag>    variable    ! comments
* tag is Fortran variable name and cannot be changed and have to be enclosed by <>
* need ! after variable, otherwise getting error.


The following variables (not pre-defined in the code) need to be defined in control file.

+------------------------+-------------------------------------------------------------------------------------------+
| tag                    | Description                                                                               |
+========================+===========================================================================================+
| <ancil_dir>            | Directory that contains ancillary data (river netowrk, remapping, and parameter namelist) |
+------------------------+-------------------------------------------------------------------------------------------+
| <input_dir>            | Directory that contains runoff data                                                       |
+------------------------+-------------------------------------------------------------------------------------------+
| <output_dir>           | Directory that contains runoff data                                                       |
+------------------------+-------------------------------------------------------------------------------------------+
| <param_nml>            | spatially constant parameter namelist (should be stored in <ancil_dir>                    |
+------------------------+-------------------------------------------------------------------------------------------+
| <sim_start>            | time of simulation start. format: yyyy-mm-dd                                              |
+------------------------+-------------------------------------------------------------------------------------------+
| <sim_end>              | time of simulation end. format:  yyyy-mm-dd                                               |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_ntopOld>        | netCDF name for River Network                                                             |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_sseg>           | dimension name for reach                                                                  |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_nhru>           | dimension name for RN_HRU                                                                 |
+------------------------+-------------------------------------------------------------------------------------------+
| <ntopWriteOption>      | logical to indicate if augmented river network is written. T or F                         |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_ntopNew>        | output netCDF name for augmented river network. See note 1                                |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_qsim>           | netCDF name for HM_HRU runoff                                                             |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_qsim>           | variable name for HM_HRU runoff                                                           |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_time>           | variable name for time                                                                    |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_hruid>          | variable name for HM_HRU ID                                                               |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_xlon>           | dimension name for x, lat, or j dimension                                                 |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_ylat>           | dimension name for y, lon, or i dimension                                                 |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_time>           | dimension name for time                                                                   |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_hruid>          | dimension name for HM_HRU                                                                 |
+------------------------+-------------------------------------------------------------------------------------------+
| <units_qsim>           | units of runoff. e.g., mm/s                                                               |
+------------------------+-------------------------------------------------------------------------------------------+
| <dt_qsim>              | time interval of runoff data in second. e.g., 86400 sec for daily step                    |
+------------------------+-------------------------------------------------------------------------------------------+
| <is_remap>             | Logical to indicate runoff needs to be remapped to RN_HRU. T or F                         |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_remap>          | netCDF name of runoff remapping                                                           |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_hruid_in_remap> | variable name for RN_HRUs                                                                 |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_weight>         | variable name for areal weights of overlapping HM_HRUs                                    |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_qhruid>         | variable name for HM_HRU ID                                                               |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_num_qhru>       | variable name for a numbers of overlapping HM_HRUs with RN_HRUs                           |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_hru_remap>      | dimension name for HM_HRU (option 2)                                                      |
+------------------------+-------------------------------------------------------------------------------------------+
| <dname_data_remap>     | dimension name for data                                                                   |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_i_index>        | dimension name for ylat index (runoff input option 3)                                     |
+------------------------+-------------------------------------------------------------------------------------------+
| <vname_j_index>        | name of xlon index dimension (runoff input option 3)                                      |
+------------------------+-------------------------------------------------------------------------------------------+
| <restart_opt>          | option to use saved restart file T->yes, F->No (start with empty channels)                |
+------------------------+-------------------------------------------------------------------------------------------+
| <route_opt>            | option for routing schemes 0-> both, 1->IRF, 2->KWT, otherwise error                      |
+------------------------+-------------------------------------------------------------------------------------------+
| <seg_outlet>           | outlet reach ID for subsetted river basin. See note 1                                     |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_output>         | output netCDF name for model simulation results                                           |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_state_in>       | input restart netCDF name.                                                                | 
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_state_out>      | output restart netCDF name.                                                               |
+------------------------+-------------------------------------------------------------------------------------------+

Variables that have default values but can be overwritten 

+------------------------+-----------------+--------------------------------------------------------------------------+
| tag                    | Default values  | Description                                                              |
+========================+=================+==========================================================================+
| <ntopWriteOption>      | F               | logical to indicate if augmented river network is written. T or F        |
+------------------------+-----------------+--------------------------------------------------------------------------+
| <fname_ntopNew>        |                 | output netCDF name for augmented river network. See note 1               |
+------------------------+-----------------+--------------------------------------------------------------------------+
| <newFileFrequency>     | annual          | frequency for new output files (day, month, or annual)                   |
+------------------------+-----------------+--------------------------------------------------------------------------+
| <hydGeometryOption>    | 1               | option for hydraulic geometry calculations (0=read from file, 1=compute) |
+------------------------+-----------------+--------------------------------------------------------------------------+
| <topoNetworkOption>    | 1               | option for network topology calculations (0=read from file, 1=compute)   |
+------------------------+-----------------+--------------------------------------------------------------------------+
| <computeReachList>     | 1               | option to compute list of upstream reaches (0=do not compute, 1=compute) |
+------------------------+-----------------+--------------------------------------------------------------------------+

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

