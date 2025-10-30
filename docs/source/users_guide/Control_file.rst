.. _Control_file:

============
Control file
============

**Control file** is a plain-text configuration file that defins all model settings, including simulation periods, file paths, routing option and more.
Each setting is called **contriol variables** and is specified using a variable name enclosed in angle brackets (``<variable_name>``).
These control variables are read at the start of the simulation (see ``./route/build/src/read_control.f90``) and
stored in fortran's public variables defined in ``./route/build/src/public_var.f90``.

🔧Key Features:
* Variables can be defined in any order in the control file, but grouping variables thematically improve readibility.
* Control variables with default values can be also included to override those defaults.

📌Syntax Rules:

* Lines starting or after exclamation mark (``!``) are treated as comments and ignored.
* Each line must follow the Format: ``<control variable name>    value    ! comments``
* Control variable name must match the corresponding Fortran variable name **exactly**.
* A exclamation mark (``!``) must appear after the value even if you don't put any comment text; otherwise getting error.
* **DO NOT include empty lines**-this will cause runtime errors.

Example control file: See ``./route/settings/SAMPLE.control`` or scroll to the bottom of this page.

.. _Basic_routing_setup:

Basic routing setup
------------------------------------------

It is difficult to instruct exactly which control variables every user needs to include, becase configurations can vary widely depending on user's goal of simulation, e.g., turn on/off lakes, water management, floodplain etc.

For a simple river routing without lake or water management-with cold start and runoff input provided at the same catchment as the river network, a user can get started with following control variables.
Control variable with None in default value must be included and assigned the appropriate value.
Additional control variables, needed for more advanced simulations, are listed after this basic set of the control variables.

For lake model option, See :doc:`Lake model option <lake>`.

For water management option, See :doc:`Water management option <Input_files>`.

+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| control variable       | Default value   | Description                                                                                             |
+========================+=================+=========================================================================================================+
| <case_name>            | None            | simulation case name. This used for output netCDF, and restart netCDF name                              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <ancil_dir>            | None            | Directory that contains ancillary data (river netowrk, remapping, and parameter namelist)               |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <input_dir>            | None            | Directory that contains runoff data                                                                     |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <output_dir>           | None            | Directory where history netCDFs are output                                                              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <restart_dir>          | <output_dir>    | Directory where restart netCDFs are output                                                              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <param_nml>            | None            | Spatially constant parameter namelist. should be stored in <ancil_dir>                                  |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <sim_start>            | None            | time of simulation start. format: yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                                     |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <sim_end>              | None            | time of simulation end. format:  yyyy-mm-dd or yyyy-mm-dd hh:mm:ss                                      |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <fname_ntopOld>        | None            | name of input netCDF for River Network. see Note 1                                                      |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_sseg>           | None            | dimension name for reach in river network netCDF                                                        |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_nhru>           | None            | dimension name for RN_HRU in river network netCDF                                                       |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <fname_qsim>           | None            | either a single runoff netCDF name, multiple netCDFs with wild card, or a text file listing netCDFs     |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <vname_qsim>           | None            | variable name for HM_HRU runoff in runoff netCDF                                                        |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <vname_hruid>          | None            | variable name for HM_HRU ID in runoff netCDF(s)                                                         |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <vname_time>           | None            | variable name for time in runoff netCDF(s)                                                              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_hruid>          | None            | dimension name for HM_HRU in runoff netCDF(s)                                                           |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_time>           | None            | dimension name for time in runoff netCDF(s)                                                             |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_xlon>           | None            | dimension name for x, lat, or j dimension in runoff netCDF(s). Only runoff is given in grid.            |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dname_ylat>           | None            | dimension name for y, lon, or i dimension in runoff netCDF(s). Only runoff is given in grid.            |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <units_qsim>           | None            | units of input runoff. depth/time where depth: m or mm, time: second, sec or s ,hour, hr or h ,day or d |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dt_ro>                | None            | time interval of runoff input in second. e.g., 86400 sec for daily step                                 |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <ro_calendar>          | None            | calenar used in runoff input netCDF. See Note 2                                                         |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <ro_time_units>        | None            | time units used in runoff input netCDF. See Note 3                                                      |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <ro_time_stamp>        | 'start'         | time stamp used in runoff input - start (default), middle, or end, otherwise error                      |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <fname_output>         | None            | output netCDF name for model simulation results. See Note 4                                             |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <outputFrequency>      | None            | time frequency used for temporal aggregation of output variables                                        |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <newFileFrequency>     | None            | frequency for new history files (daily, monthly, yearly, single)                                        |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <histTimeStamp_offset> | 0               | time stamp offset [second] from a start of time step                                                    |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <route_opt>            | 0               | routing schem options: 0->Sum, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW, otherwise error. See Note 5         |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <floodplain>           | F               | logical to add floodplain to main channel. flood water is computed if floodplain is added.              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <dt_qsim>              | None            | time interval of simulation in second. e.g., 86400 sec for daily step                                   |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <hw_drain_point>       | 2               | how to add catchment runoff to reach for headwater HRUs. 1->top of reach, 2->bottom of reach (default)  |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <min_length_route>     | 0               | min. reach length for routing to be performed. pass-through if reach length less than this threshold    |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <doesBasinRoute>       | 1               | hillslope routing options. 0-> no (already routed), 1->IRF                                              |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <hydGeometryOption>    | 1               | option for hydraulic geometry calculations (0=read from file, 1=compute)                                |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <topoNetworkOption>    | 1               | option for network topology calculations (0=read from file, 1=compute)                                  |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+
| <computeReachList>     | 1               | option to compute list of upstream reaches (0=do not compute, 1=compute)                                |
+------------------------+-----------------+---------------------------------------------------------------------------------------------------------+

Terminologies: RN_HRU=River Network HRU (Hydrologic Response Unit or simply catchment), HM_HRU=hydrologic model (or forcing) HRU. HRU can be grid box. "Forcing" for river model means runoff, evaporation and precipitation for lake routing, solutes for solute transport

1. Often river network data has different variable names than defaults. In this case, variable names can be speficied in control file as well. See :doc:`River parameters <Riv>`.

2. Calendar in runoff input time should be read from netCDF, but If runoff input netCDF does not have calendar attribute, it can be specified. Make sure time variable in runoff data use either ``noleap``, ``standard``, ``gregorian``, or ``proleptic_gregorian``. case insensitive

3. Like Calendar, If runoff input netCDF does not have time unit attribute, it can be specified. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted. <unit> can be days, hours, minutes, seconds.

4. routing physics option: route_opt

  * it is possible to specify multiple options (e.g., 0125 -> run with SUM, IRF KWT and DW).

5 Restrictions related to history output: dt_qsim, outputFrequency and newFileFrequency

  * dt_qsim (simulation time step) must be less than 86400 sec (one day). Muskingum-Cunge method will run at much shorter time step. Other methods can run at this time step, but Diffusive wave routing produce the results with less errors at shorter time step.

  * dt_qsim can be different time step than input time step.

  * outputFrequency can be integer numeric (e.g, 1, 2 etc), which is interpreted as a number of simulation time steps for temporal aggregation of the history flux variables, or literal (daily, monthly yearly).
    The numeric outputFrequency can be used for sub-daily dt_qsim, and remainder of 86400 divided by numeric x dt_qsim must be zero. For example, if dt_qsim is 10800 sec (=3hr), accepted outputFrequency are
    1, 2, 4, 8

  * newFileFrequency must be the same as or shorter than outputFrequency. For example, with monthly outputFrequency, newFileFrequency must be monthly, yearly or single

  * The abovementioned restrictions are check in the code, so any violations are notified as error and the program is terminated.


.. _Control_file_basic_examples:

Control file basic examples
----------------------------

These are examples for three cases of runoff input. These are just templates to start with.
Users need to specify appropreate directories, netCDF variables/dimension names based on their data

Option 1 - runoff input is given at RN_HRU

::

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
  <dt_rof>            86400                           ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE
  ! ----------------------------------
  <is_remap>          F                               ! logical to indicate runnoff needs to be mapped to river network HRU
  ! **************************************************************************************************************************
  ! Namelist file name
  ! ---------------------------
  <param_nml>         param.nml.default               ! spatially constant model parameters
  ! **************************************************************************************************************************

Option 2 - runoff input is given at HM_HRU

::

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
  <dt_rof>                86400                            ! time interval of the runoff
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

Option 3 - runoff input is given at grid

::

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
  <dt_rof>                86400                            ! time interval of the runoff
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
