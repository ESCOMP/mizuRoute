Control file
============

Control file is a simple text file, mainly defining model control such as simulation time, file name and locations, routing option etc. 
Variables in control file are read in the beginning of the code (see ``./build/src/read_control.f90``) and 
saved in fortran variable specified by tag (inside <> in table) and as public variables (see ``./build/src/public_var.f90``) . 
Some of such public varialbes have some default values, but most of them are not defined.
Those undefined variables need to be defined in control file.   
Other variables in supplement table have their default values but can be also included in control file to overwrite the values. 
The order of variables in the control file does not matter. However, grouping variables into similar themes is recommended for readibility. 

Minimum required variables depends on runoff input options.

Example of control file is given in ``./route/settings/SAMPLE.control`` or see Examples at bottom of this page.

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
| 1,2,3  | <dname_sseg>           | dimension name for reach in river network netCDF                                          |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_nhru>           | dimension name for RN_HRU in river network netCDF                                         |
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
| 1,2,3  | <units_qsim>           | units of input runoff. e.g., mm/s                                                         |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <dt_qsim>              | time interval of input runoff in second. e.g., 86400 sec for daily step                   |
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
|     3  | <vname_i_index>        | variable name of ylat index                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|     3  | <vname_j_index>        | variable name of xlon index                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <vname_num_qhru>       | variable name for a numbers of overlapping HM_HRUs with RN_HRUs                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <dname_hru_remap>      | dimension name for HM_HRU                                                                 |
+--------+------------------------+-------------------------------------------------------------------------------------------+
|   2,3  | <dname_data_remap>     | dimension name for data                                                                   |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_output>         | output netCDF name for model simulation results                                           |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_state_in>       | input restart netCDF name. No need to be specified unless <restart_opt> is T.             | 
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_state_out>      | output restart netCDF name.                                                               |
+--------+------------------------+-------------------------------------------------------------------------------------------+
| 1,2,3  | <route_opt>            | option for routing schemes 0-> both, 1->IRF, 2->KWT, otherwise error                      |
+--------+------------------------+-------------------------------------------------------------------------------------------+

Variables that have default values but can be overwritten 

+------------------------+------------------------+--------------------------------------------------------------------------+
| tag                    | Default values         | Description                                                              |
+========================+========================+==========================================================================+
| <ntopAugmentMode>      | F                      | logical to indicate river network augmention mode. See note 1.           |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <seg_outlet>           | -9999                  | outlet reach ID for subsetted river basin. See note 1                    |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <fname_ntopNew>        | <fname_ntopOld>_new.nc | output netCDF name for augmented river network. See note 1               |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <newFileFrequency>     | annual                 | frequency for new output files (day, month, or annual)                   |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <hydGeometryOption>    | 1                      | option for hydraulic geometry calculations (0=read from file, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <topoNetworkOption>    | 1                      | option for network topology calculations (0=read from file, 1=compute)   |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <computeReachList>     | 1                      | option to compute list of upstream reaches (0=do not compute, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <doesBasinRoute>       | 1                      | hillslope routing options. 0-> no (already routed), 1->IRF               |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <restart_opt>          | F                      | option to use restart file T->yes, F->No (start with empty channels)     |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <calendar>             | From runoff input      | calenar used in input runoff data. See note 2.                           |
+------------------------+------------------------+--------------------------------------------------------------------------+

#. River network subset mode. 

  * if <seg_outlet> is given, the river network topology and parameters read from <fname_ntopOld> are written in <fname_ntopNew> and the program stops. 
 
#. River network augmentation mode. 

  * All the computed river network topology and parameters are written in <fname_ntopNew> and the program stops. 

#. If runoff input data does not have calendar attribute, it can be specified. Make sure time variable in runoff data use either ``noleap``, ``standard``, ``gregorian``, or ``proleptic_gregorian``. case sensitive


Often case, river network data has different variable names than defaults. In this case, variable names can be speficied in control file as well.
See :doc:`River parameters <seg_hru_param>`.   


Control file examples
---------------------

These are examples for three cases of runoff input. These are just templates to start with. 
Users need to specify appropreate directories, netCDF variables/dimension names based on their data

Option 1 - runoff input is given at RN_HRU::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!" 
  !
  ! *************************************************************************************************************************
  ! DEFINE DIRECTORIES 
  ! --------------------------
  <ancil_dir>         ./ancillary_data/               ! directory containing ancillary data (river network, remapping netCDF)
  <input_dir>         ./input/                        ! directory containing input data (runoff netCDF)
  <output_dir>        ./output/                       ! directory containing output data
  ! *************************************************************************************************************************
  ! DEFINE TIME PERIOD OF THE SIMULATION
  ! --------------------------------------------
  <sim_start>         1950-1-1 00:00:00               ! time of simulation start (year-month-day hh:mm:ss)
  <sim_end>           1951-1-1 00:00:00               ! time of simulation end (year-month-day hh:mm:ss)
  ! **************************************************************************************************************************
  ! DEFINE FINE NAME AND DIMENSIONS
  ! ---------------------------------------
  <fname_ntopOld>     ntopo_entire.nc                 ! name of netCDF containing river segment data 
  <dname_sseg>        seg                             ! dimension name of the stream segments
  <dname_nhru>        hru                             ! dimension name of the RN_HRUs
  ! **************************************************************************************************************************
  ! DEFINE DESIRED VARIABLES FOR THE NETWORK TOPOLOGY
  ! ---------------------------------------------------------
  <seg_outlet>        -9999                           ! reach ID of outlet streamflow segment. -9999 for all segments 
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>        runoff.RN_HRU.nc                ! name of netCDF containing the runoff
  <vname_qsim>        RUNOFF                          ! variable name of HRU runoff
  <vname_time>        time                            ! variable name of time in the runoff file
  <vname_hruid>       hru                             ! variable name of runoff HRU ID
  <dname_time>        time                            ! dimension name of time
  <dname_hruid>       hru                             ! dimension name of HM_HRU
  <units_qsim>        mm/s                            ! units of runoff
  <dt_qsim>           86400                           ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE 
  ! ----------------------------------
  <is_remap>          F                               ! logical to indicate runnoff needs to be mapped to river network HRU 
  ! **************************************************************************************************************************
  ! DEFINE RUN CONTROL 
  ! ---------------------------
  <route_opt>         0                               ! option for routing schemes 0-> both, 1->IRF, 2->KWT otherwise error 
  <restart_opt>           F                                 ! option to use saved flow states T->yes, F->No 
  ! **************************************************************************************************************************
  ! DEFINE OUTPUT FILE
  ! ---------------------------
  <fname_output>      flow                            ! prefix of output netCDF name (netCDF name = flow_nomapping_yyyy.nc)
  <fname_state_in>    state.in.nc                     ! netCDF name for the model state input 
  <fname_state_out>   state.out.nc                    ! netCDF name for the channel state output 
  ! **************************************************************************************************************************
  ! Namelist file name 
  ! ---------------------------
  <param_nml>         param.nml.default               ! spatially constant model parameters    
  ! **************************************************************************************************************************

Option 2 - runoff input is given at HM_HRU::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!" 
  !
  ! *************************************************************************************************************************
  ! DEFINE DIRECTORIES 
  ! --------------------------
  <ancil_dir>             ./ancillary_data/                ! directory containing ancillary data (river network, remapping netCDF)
  <input_dir>             ./input/                         ! directory containing input data (runoff netCDF)
  <output_dir>            ./output/                        ! directory containing output data
  ! *************************************************************************************************************************
  ! DEFINE TIME PERIOD OF THE SIMULATION
  ! --------------------------------------------
  <sim_start>             1950-1-1 00:00:00                ! time of simulation start (year-month-day hh:mm:ss)
  <sim_end>               1951-1-1 00:00:00                ! time of simulation end (year-month-day hh:mm:ss)
  ! **************************************************************************************************************************
  ! DEFINE FINE NAME AND DIMENSIONS
  ! ---------------------------------------
  <fname_ntopOld>         ntopo_entire.nc                  ! name of netCDF containing river segment data 
  <dname_sseg>            seg                              ! dimension name of the stream segments
  <dname_nhru>            hru                              ! dimension name of the RN_HRUs
  ! **************************************************************************************************************************
  ! DEFINE DESIRED VARIABLES FOR THE NETWORK TOPOLOGY
  ! ---------------------------------------------------------
  <seg_outlet>            -9999                            ! reach ID of outlet streamflow segment. -9999 for all segments 
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>            runoff.HM_HRU.nc                 ! name of netCDF containing the HRU runoff
  <vname_qsim>            RUNOFF                           ! variable name of HRU runoff
  <vname_time>            time                             ! variable name of time in the runoff file
  <vname_hruid>           hru                              ! variable name of runoff HRU ID
  <dname_time>            time                             ! dimension name of time
  <dname_hruid>           hru                              ! dimension name of HM_HRU
  <units_qsim>            mm/s                             ! units of runoff
  <dt_qsim>               86400                            ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE 
  ! ----------------------------------
  <is_remap>              T                                 ! logical to indicate runnoff needs to be mapped to RN_HRU 
  <fname_remap>           spatialweights_HM_HRU_RN_HRU.nc   ! name of netCDF for HM_HRU-RN_HRU mapping data
  <vname_hruid_in_remap>  polyid                            ! variable name of RN_HRU in the mapping file
  <vname_weight>          weight                            ! variable name of areal weights of overlapping HM_HUs for each RN_HRU
  <vname_qhruid>          overlapPolyId                     ! variable name of HM_HRU ID
  <vname_num_qhru>        overlaps                          ! variable name of numbers of HM_HRUs for each RN_HRU
  <dname_hru_remap>       polyid                            ! dimension name of RN_HRU (in the mapping file)
  <dname_data_remap>      data                              ! dimension name of ragged HM_HRU
  ! **************************************************************************************************************************
  ! DEFINE RUN CONTROL 
  ! ---------------------------
  <route_opt>             0                                 ! option for routing schemes 0-> both, 1->IRF, 2->KWT otherwise error 
  <restart_opt>           F                                 ! option to use saved flow states T->yes, F->No 
  ! **************************************************************************************************************************
  ! DEFINE OUTPUT FILE
  ! ---------------------------
  <fname_output>          flow                              ! prefix of output netCDF name (netCDF name = flow_nomapping_yyyy.nc)
  <fname_state_in>        state.in.nc                       ! netCDF name for the model state input 
  <fname_state_out>       state.out.nc                      ! netCDF name for the channel state output 
  ! **************************************************************************************************************************
  ! Namelist file name 
  ! ---------------------------
  <param_nml>             param.nml.default                 ! spatially constant model parameters    
  ! **************************************************************************************************************************

Option 3 - runoff input is given at grid::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!" 
  !
  ! *************************************************************************************************************************
  ! DEFINE DIRECTORIES 
  ! --------------------------
  <ancil_dir>             ./ancillary_data/                ! directory containing ancillary data (river network, remapping netCDF)
  <input_dir>             ./input/                         ! directory containing input data (runoff netCDF)
  <output_dir>            ./output/                        ! directory containing output data
  ! *************************************************************************************************************************
  ! DEFINE TIME PERIOD OF THE SIMULATION
  ! --------------------------------------------
  <sim_start>             1950-1-1 00:00:00                ! time of simulation start (year-month-day hh:mm:ss)
  <sim_end>               1951-1-1 00:00:00                ! time of simulation end (year-month-day hh:mm:ss)
  ! **************************************************************************************************************************
  ! DEFINE FINE NAME AND DIMENSIONS
  ! ---------------------------------------
  <fname_ntopOld>         ntopo_entire.nc                  ! name of netCDF containing river segment data 
  <dname_sseg>            seg                              ! dimension name of the stream segments
  <dname_nhru>            hru                              ! dimension name of the RN_HRUs
  ! **************************************************************************************************************************
  ! DEFINE DESIRED VARIABLES FOR THE NETWORK TOPOLOGY
  ! ---------------------------------------------------------
  <seg_outlet>            -9999                            ! reach ID of outlet streamflow segment. -9999 for all segments 
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>            runoff.HM_HRU.nc                 ! name of netCDF containing the HRU runoff
  <vname_qsim>            RUNOFF                           ! variable name of HRU runoff
  <vname_time>            time                             ! variable name of time in the runoff file
  <dname_time>            time                             ! dimension name of time
  <dname_xlon>            lon                              ! dimension name of x(j)
  <dname_ylat>            lat                              ! dimension name of y(i)
  <units_qsim>            mm/s                             ! units of runoff
  <dt_qsim>               86400                            ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE 
  ! ----------------------------------
  <is_remap>              T                                 ! logical to indicate runnoff needs to be mapped to RN_HRU 
  <fname_remap>           spatialweights_HM_HRU_RN_HRU.nc   ! name of netCDF for HM_HRU-RN_HRU mapping data
  <vname_hruid_in_remap>  polyid                            ! variable name of RN_HRU in the mapping file
  <vname_weight>          weight                            ! variable name of areal weights of overlapping HM_HUs for each RN_HRU
  <vname_i_index>         i_index                           ! variable name of ylat index
  <vname_j_index>         j_index                           ! variable name of xlon index
  <vname_num_qhru>        overlaps                          ! variable name of numbers of HM_HRUs for each RN_HRU
  <dname_hru_remap>       polyid                            ! dimension name of RN_HRU (in the mapping file)
  <dname_data_remap>      data                              ! dimension name of ragged HM_HRU
  ! **************************************************************************************************************************
  ! DEFINE RUN CONTROL 
  ! ---------------------------
  <route_opt>             0                                 ! option for routing schemes 0-> both, 1->IRF, 2->KWT otherwise error 
  <restart_opt>           F                                 ! option to use saved flow states T->yes, F->No 
  ! **************************************************************************************************************************
  ! DEFINE OUTPUT FILE
  ! ---------------------------
  <fname_output>          flow                              ! prefix of output netCDF name (netCDF name = flow_nomapping_yyyy.nc)
  <fname_state_in>        state.in.nc                       ! netCDF name for the model state input 
  <fname_state_out>       state.out.nc                      ! netCDF name for the channel state output 
  ! **************************************************************************************************************************
  ! Namelist file name 
  ! ---------------------------
  <param_nml>             param.nml.default                 ! spatially constant model parameters    
  ! **************************************************************************************************************************
