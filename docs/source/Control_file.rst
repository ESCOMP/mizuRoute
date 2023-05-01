============
Control file
============

Control file is a text file, mainly defining model control such as simulation time, file name and locations, routing option etc.
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

+--------+------------------------+----------------------------------------------------------------------------------------------------+
| option | tag                    | Description                                                                                        |
+========+========================+====================================================================================================+
| 1,2,3  | <case_name>            | simulation case name. This used for output netCDF, and restart netCDF name                         |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <ancil_dir>            | Directory that contains ancillary data (river netowrk, remapping, and parameter namelist)          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <input_dir>            | Directory that contains runoff data                                                                |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <output_dir>           | Directory that contains runoff data                                                                |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <param_nml>            | Spatially constant parameter namelist (should be stored in <ancil_dir>                             |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <sim_start>            | time of simulation start. format: yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                                |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <sim_end>              | time of simulation end. format:  yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                                 |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_ntopOld>        | name of input netCDF for River Network                                                             |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_sseg>           | dimension name for reach in river network netCDF                                                   |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_nhru>           | dimension name for RN_HRU in river network netCDF                                                  |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_qsim>           | netCDF name for HM_HRU runoff                                                                      |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <vname_qsim>           | variable name for HM_HRU runoff                                                                    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <vname_time>           | variable name for time                                                                             |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2    | <vname_hruid>          | variable name for HM_HRU ID                                                                        |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2    | <dname_time>           | dimension name for time                                                                            |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|     3  | <dname_xlon>           | dimension name for x, lat, or j dimension                                                          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|     3  | <dname_ylat>           | dimension name for y, lon, or i dimension                                                          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <dname_hruid>          | dimension name for HM_HRU                                                                          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <units_qsim>           | units of input runoff. e.g., mm/s                                                                  |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <dt_qsim>              | time interval of simulation in second. e.g., 86400 sec for daily step                              |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <is_remap>             | Logical to indicate runoff needs to be remapped to RN_HRU. T or F                                  |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <fname_remap>          | netCDF name of runoff remapping                                                                    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <vname_hruid_in_remap> | variable name for RN_HRUs                                                                          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <vname_weight>         | variable name for areal weights of overlapping HM_HRUs                                             |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2    | <vname_qhruid>         | variable name for HM_HRU ID                                                                        |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|     3  | <vname_i_index>        | variable name of ylat index                                                                        |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|     3  | <vname_j_index>        | variable name of xlon index                                                                        |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <vname_num_qhru>       | variable name for a numbers of overlapping HM_HRUs with RN_HRUs                                    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <dname_hru_remap>      | dimension name for HM_HRU                                                                          |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <dname_data_remap>     | dimension name for data                                                                            |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <fname_output>         | output netCDF name for model simulation results                                                    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <outputFrequency>      | time frequency used for temporal aggregation of output variables                                   |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <newFileFrequency>     | frequency for new history files (daily, monthly, yearly, single)                                   |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
| 1,2,3  | <route_opt>            | routing schem options: 0->Sum, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW, otherwise error. see Note 1    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+

1. routing option: route_opt

  * it is possible to specify multiple options (e.g., 0125 -> run with SUM, IRF KWT and DW).


2 Restrictions related to history output: dt_qsim, outputFrequency and newFileFrequency

  * dt_qsim (simulation time step) must be less than 86400 sec (one day). Muskingum-Cunge method will run at much shorter time step. Other methods can run at this time step, but Diffusive wave routing produce the results with less errors at shorter time step.

  * dt_qsim can be different time step than input time step.

  * outputFrequency can be integer numeric (e.g, 1, 2 etc), which is interpreted as a number of simulation time steps for temporal aggregation of the history flux variables, or literal (daily, monthly yearly).
    The numeric outputFrequency can be used for sub-daily dt_qsim, and remainder of 86400 divided by numeric x dt_qsim must be zero. For example, if dt_qsim is 10800 sec (=3hr), accepted outputFrequency are
    1, 2, 4, 8 

  * newFileFrequency must be the same as or shorter than outputFrequency. For example, with monthly outputFrequency, newFileFrequency must be monthly, yearly or single
   
  * The abovementioned restrictions are check in the code, so any violations are notified as error and the program is terminated.

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
| <hydGeometryOption>    | 1                      | option for hydraulic geometry calculations (0=read from file, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <topoNetworkOption>    | 1                      | option for network topology calculations (0=read from file, 1=compute)   |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <computeReachList>     | 1                      | option to compute list of upstream reaches (0=do not compute, 1=compute) |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <doesBasinRoute>       | 1                      | hillslope routing options. 0-> no (already routed), 1->IRF               |
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

Restart options 
---------------------

mizuRoute does not write restart netCDF as default. The following control variables are used to control restart dropoff timing and use restart file for continuous run from the previous simulations.
The restart file is written at previous time step to the specified time. In other words, if ``Specified`` is used for <restart_write> and ``1981-01-01-00000`` is specified in <restart_date>, mizuRoute writes restart file
at ``1980-12-31 00:00:00`` for daily time step. The restart file name uses the time stamp at user specified timing. ``yearly``, ``monthly``, ``daily`` options also follow this convention.

The restart file name convension:  <case_name>.r.yyyy-mm-dd-sssss.nc


+---------------------+---------------------------------------------------------------------------------------------------------+
| tag                 | Description                                                                                             |
+=====================+=========================================================================================================+
| <restart_dir>       | directory for restart files. defualt is <output_dir>                                                    |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_write>     | restart ouput options. never (default), last, specified, yearly, monthly, daily.                        |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_date>      | restart time in yyyy-mm-dd (hh:mm:ss). required if <restart_write> = "Specified"                        |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_month>     | periodic restart month (default 1). Effective if <restart_write>="yearly"                               |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_day>       | periodic restart day (default 1). Effective if <restart_write>="yearly" or "monthly"                    |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_hour>      | periodic restart hour (default 0). Effective if <restart_write>="yearly", "monthly", or "daily"         |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <fname_state_in>    | input restart netCDF name. If not specified, simulation start with cold start                           |
+---------------------+---------------------------------------------------------------------------------------------------------+


History output variables
---------------------

The following variables, besides time, basinID (RN_hru ID) and reachID can be output in netCDF. Users can control which variables are output by setting <variable_name> to T or F in control file. All the variables are set to T by default.
The output file name includes a timie stamp at the first time step.

The output file name convension:  <case_name>.h.yyyy-mm-dd-sssss.nc


+------------------------+------------------------------------------------------------------------------------------------+
| output variables       | Descriptions                                                                                   |
+========================+================================================================================================+
| <basRunoff>            | runoff depth at RN_hru, remapped from HM_hru. See note 1 and 2.                                |
+------------------------+------------------------------------------------------------------------------------------------+
| <instRunoff>           | runoff volume [m3/s] at reach, converted by mulitplying basRunoff by RN_hru area . See note 2  |
+------------------------+------------------------------------------------------------------------------------------------+
| <dlayRunoff>           | runoff volume [m3/s] at reach, after hillslope routing instRunoff. see Note 2                  |
+------------------------+------------------------------------------------------------------------------------------------+
| <sumUpstreamRunoff>    | accumulated delayed runoff volume (dlyRunoff) over all upstream reaches.                       |
+------------------------+------------------------------------------------------------------------------------------------+
| <KWTroutedRunoff>      | runoff volume [m3/s] after Kinematic wave tracking (KWT) reach routing dlayRunoff. See note 3  |
+------------------------+------------------------------------------------------------------------------------------------+
| <IRFroutedRunoff>      | runoff volume [m3/s] after IRF reach routing dlayRunoff. See note 3                            |
+------------------------+------------------------------------------------------------------------------------------------+
| <KWroutedRunoff>       | runoff volume [m3/s] after KW (Kinematic Wave) reach routing dlayRunoff. See note 3            |
+------------------------+------------------------------------------------------------------------------------------------+
| <MCroutedRunoff>       | runoff volume [m3/s] after MC (Muskingum-Cunge) reach routing dlayRunoff. See note 3           |
+------------------------+------------------------------------------------------------------------------------------------+
| <DWroutedRunoff>       | runoff volume [m3/s] after DW (Diffusive wave) reach routing dlayRunoff. See note 3            |
+------------------------+------------------------------------------------------------------------------------------------+

1. The unit of runoff depth is the same as the unit used in runoff data.


2. If runoff depth from runoff data is already delayed by hill-slope routing outside mizuRoute, <doesBasinRoute> should be set to 0. In this case, runoff volume computed from basRunoff is populated in <dlayRunoff> and <instRunoff> is not output.


3. routed runoff corresponding to the scheme is not ouput if users deactivate a particular routing scheme with <route_opt> tag.


gauge data options
---------------------

mizuRoute can read gauge observed discharge data (in netCDF) along with gauge meta ascii data. To read gauge observation and gauge metadata, the following control variables need to be specified.
gauge meta ascii file is csv format, and  should include at least gauge id and corresponding reach id
gauge discharge data is used for data assimilation (current version does not include this at this moment)


+---------------------+---------------------------------------------------------------------------------------------------------+
| tag                 | Description                                                                                             |
+=====================+=========================================================================================================+
| <gageMetaFile>      | gauge meta file (two column csv format): gauge_id (non-numeric ID is accepted), seg_id                  |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <outputAtGage>      | logical value (T or F) to limit history variable output at gauge reaches.                               |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <fname_gageObs>     | gauge discharge data                                                                                    |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <vname_gageFlow>    | variable name for discharge [m3/s]                                                                      |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <vname_gageSite>    | variable name for gauge site name (character array)                                                     |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <vname_gageTime>    | variable name for time                                                                                  |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <dname_gageSite>    | dimension name for site                                                                                 |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <dname_gageTime>    | dimension name for time                                                                                 |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <strlen_gageSite>   | maximum gauge name string length                                                                        |
+---------------------+---------------------------------------------------------------------------------------------------------+



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
  ! DEFINE SIMULATION CONTROLS
  ! --------------------------------------------
  <case_name>             cameo_opt1                               ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00                      ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00                      ! time of simulation end.   year-month-day (hh:mm:ss)
  <fname_state_in>        cameo_opt1.mizuRoute.r.1950-1-1-00000.nc ! netCDF name for the model state input
  <restart_write>         specified                                ! restart write option. never, last, specified (need to specify date with <restart_date>
  <restart_date>          1950-08-31 00:00:00                      ! desired restart starting datetime
  <route_opt>             012345                                   ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <dt_qsim>               3600                                     ! 1 hour simulation
  <newFileFrequency>      daily                                    ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                                    ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
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
  ! DEFINE SIMULATION CONTROLS
  ! --------------------------------------------
  <case_name>             cameo_opt2                               ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00                      ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00                      ! time of simulation end.   year-month-day (hh:mm:ss)
  <fname_state_in>        cameo_opt2.mizuRoute.r.1950-1-1-00000.nc ! netCDF name for the model state input
  <restart_write>         specified                                ! restart write option. never, last, specified (need to specify date with <restart_date>
  <restart_date>          1950-08-31 00:00:00                      ! desired restart starting datetime
  <route_opt>             012345                                   ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <dt_qsim>               3600                                     ! 1 hour simulation
  <newFileFrequency>      daily                                    ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                                    ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
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
  ! DEFINE SIMULATION CONTROLS
  ! --------------------------------------------
  <case_name>             cameo_opt3                               ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00                      ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00                      ! time of simulation end.   year-month-day (hh:mm:ss)
  <fname_state_in>        cameo_opt3.mizuRoute.r.1950-1-1-00000.nc ! netCDF name for the model state input
  <restart_write>         specified                                ! restart write option. never, last, specified (need to specify date with <restart_date>
  <restart_date>          1950-08-31 00:00:00                      ! desired restart starting datetime
  <route_opt>             012345                                   ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <dt_qsim>               3600                                     ! 1 hour simulation
  <newFileFrequency>      daily                                    ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                                    ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
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
  ! Namelist file name 
  ! ---------------------------
  <param_nml>             param.nml.default                 ! spatially constant model parameters
  ! **************************************************************************************************************************
