.. _Control_file:

============
Control file
============

**Control file** is a plain-text file that defins all model settings, including simulation periods, file paths, routing option and more.
Each setting is called **contriol keys** or **control variables** and is specified using a variable name enclosed in angle brackets (``<control_key_name>``).
These control keys are read at the start of the simulation (see ``./route/build/src/read_control.f90``) and
stored in fortran's public variables defined in ``./route/build/src/public_var.f90``.

🔧Features:

* Control keys can be defined in any order in the control file, but grouping variables thematically improve readibility.
* Control keys with default values can be also included to override those defaults.
* All the control key names match the corresponding Fortran variable names in ``./route/build/src/public_var.f90``.

📌Syntax Rules:

* Lines starting or after exclamation mark (``!``) are treated as comments and ignored.
* Each line must follow the Format: ``<control_key_name>    value    ! comments``. A exclamation mark (``!``) must appear after the value even if you don't put any comment text; otherwise getting error.
* **DO NOT include empty lines**-this will cause runtime errors.


.. _Basic_routing_setup:

Basic routing setup
------------------------------------------

It is difficult to tell exactly which control keys every user needs to include, becase configurations can vary widely depending on user's goal of simulation, e.g., turn on/off lakes, water management, floodplain etc.

For a river routing without lake or water management with cold start and runoff input provided at the same catchment as the river network, a user can get started with following control keys.
See :ref:`Control file basic examples <Control_file_basic_examples>`. Also, an example control file including additional control keys (for advanced features) are provided under ``./route/settings/SAMPLE.control``.
In the table below, control keys with None in default value must be defined with the appropriate value.
The netCDF variables (e.g.,  runoff and river data netCDFs) have their default names; if their netCDF variables do not match the default names, a user needs to define in the control file.
Please see :doc:`Input files <Input_files>`.

Additional control variables, needed for more advanced simulations, are listed after this basic set of the control keys.

For lake model option, See :doc:`Lake model option <lake>`.

For water management option, See :doc:`Water management option <Input_files>`.


.. list-table:: Control keys to setup basic routing
   :header-rows: 1
   :widths: 20 15 15 50
   :name: basic routing configuration

   * - Control key
     - Type
     - Default
     - Description
   * - ``<case_name>``
     - char
     - None
     - simulation case name. This used for output netCDF, and restart netCDF name
   * - ``<ancil_dir>``
     - char
     - None
     - Directory that contains ancillary data (river netowrk, remapping, and parameter namelist)
   * - ``<input_dir>``
     - char
     - None
     - Directory that contains runoff data
   * - ``<output_dir>``
     - char
     - None
     - Directory where history netCDFs are output
   * - ``<sim_start>``
     - char
     - None
     - time of simulation start. format: yyyy-mm-dd or yyyy-mm-dd hh:mm:ss
   * - ``<sim_end>``
     - char
     - None
     - time of simulation end. format:  yyyy-mm-dd or yyyy-mm-dd hh:mm:ss
   * - ``<fname_ntopOld>``
     - char
     - None
     - name of input netCDF for River Network. see Note 1
   * - ``<fname_qsim>``
     - char
     - None
     - either a single runoff netCDF name, multiple netCDFs with wild card, or a text file listing netCDFs. See Note 2
   * - ``<units_qsim>``
     - char
     - None
     - units of input runoff. depth/time where depth: m or mm, time: second, sec or s ,hour, hr or h ,day or d
   * - ``<dt_ro>``
     - real
     - None
     - time interval of runoff input in second. e.g., 86400 sec for daily step
   * - ``<ro_calendar>``
     - char
     - None
     - calenar used in runoff input netCDF. See Note 3
   * - ``<ro_time_units>``
     - char
     - None
     - time units used in runoff input netCDF. See Note 4
   * - ``<ro_time_stamp>``
     - char
     - start
     - time stamp used in runoff input - start (default), middle, or end, otherwise error
   * - ``<outputFrequency>``
     - char
     - None
     - time frequency used for temporal aggregation of output variables. See Note 5
   * - ``<newFileFrequency>``
     - char
     - None
     - frequency for new history files (daily, monthly, yearly, single). See Note 5
   * - ``<histTimeStamp_offset>``
     - real
     - 0
     - time stamp offset [second] from a start of time step
   * - ``<route_opt>``
     - int
     - 1
     - routing schem options: 0->Sum, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW, otherwise error. See :doc:`River routing configuration <Riv>`
   * - ``<dt_qsim>``
     - int
     - None
     - time interval of simulation in second. e.g., 86400 sec for daily step
   * - ``<param_nml>``
     - char
     - None
     - Spatially constant parameter namelist. should be stored in <ancil_dir>. See Note 6


1. NetCDF variable and dimension names in river data has their default names, but a user can use different names than defaults. In this case, variable or/and dimension names can be speficied in control file as well. See :ref:`River parameters <River_network_data>`.

2. The same as Note 1 except for runoff input netCDF. See :ref:`Runoff input files <Runoff_data>`.

3. Calendar in runoff input time should be read from netCDF, but If runoff input netCDF does not have calendar attribute for some reason, it can be specified. Make sure time variable in runoff data use either ``noleap``, ``standard``, ``gregorian``, or ``proleptic_gregorian``. case insensitive

4. Like Calendar, If runoff input netCDF does not have time unit attribute, it can be specified. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted. <unit> can be days, hours, minutes, seconds.

5. Restrictions related to history output: dt_qsim, outputFrequency and newFileFrequency

  * ``dt_qsim`` (simulation time step) must be less than 86400 sec (one day). Muskingum-Cunge method will run at much shorter time step. Other methods can run at this time step, but Diffusive wave routing produce the results with less errors at shorter time step.

  * ``dt_qsim`` can be different time step than input time step.

  * ``outputFrequency`` can be integer numeric (e.g, 1, 2 etc), which is interpreted as a number of simulation time steps for temporal aggregation of the history flux variables, or literal (daily, monthly yearly).
    The numeric outputFrequency can be used for sub-daily dt_qsim, and remainder of 86400 divided by numeric x ``dt_qsim`` must be zero. For example, if ``dt_qsim`` is 10800 sec (=3hr), accepted ``outputFrequency`` are
    1, 2, 4, 8

  * ``newFileFrequency`` must be the same (output netcdf contains only one time step) as or longer than ``outputFrequency``. For example, with monthly ``outputFrequency``, ``newFileFrequency`` must be monthly, yearly or single

  * The abovementioned restrictions are check in the code, so any violations are notified as error and the program is terminated.

6. Spatially constant parameters are provided in a namelist. See :ref:`Spatially-constant parameter namelist <namelist_file>`. Use the namelist provided in github as a template.


.. _Control_file_basic_examples:

Control file basic examples
----------------------------

These are examples for three cases of runoff input (see :ref:`Runoff file(s) <Runoff_data>` for runoff input cases). These examples are provided in testCase (:ref:`testCase <testCase_data>`).
For a example of more comple set of the control keys, see ``./route/settings/SAMPLE.control``


Option 1 - runoff input is given at river network HRUs

::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!"
  !    Do not leave lines without "!" at the beginning of the lines
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
  <case_name>             cameo_case1                 ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00         ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00         ! time of simulation end.   year-month-day (hh:mm:ss)
  <dt_qsim>               3600                        ! simulation time interval [sec]
  <newFileFrequency>      monthly                     ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                       ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
  ! ****************************************************************************************************************************
  ! ROUTE options
  ! --------------------------
  <route_opt>             012345                      ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <doesBasinRoute>        1                           ! basin routing?  0-> no, 1-> basin UH
  ! **************************************************************************************************************************
  ! DEFINE RIVER DATA
  ! ---------------------------------------
  <fname_ntopOld>     ntopo_nhdplus_cameo_pfaf.nc     ! name of netCDF containing river segment data
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>        RUNOFF_case1.nc                 ! name of netCDF containing the runoff
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


Option 2 - runoff input is given at hydrologic model HRUs (not the same as river network HRUs).
Runoff and mapping netCDF variables are shown just to show how netCDFs are structured, and commented out in this example, since they are the same as the default names.


::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!"
  !    Do not leave lines without "!" at the beginning of the lines
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
  <case_name>             cameo_case2                      ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00              ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00              ! time of simulation end.   year-month-day (hh:mm:ss)
  <dt_qsim>               3600                             ! simulation time interval [sec]
  <newFileFrequency>      monthly                          ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                            ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
  ! ****************************************************************************************************************************
  ! ROUTE options
  ! --------------------------
  <route_opt>             012345                           ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <doesBasinRoute>        1                                ! basin routing?  0-> no, 1-> basin UH
  ! **************************************************************************************************************************
  ! DEFINE RIVER DATA
  ! ---------------------------------------
  <fname_ntopOld>         ntopo_nhdplus_cameo_pfaf.nc      ! name of netCDF containing river segment data
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>            RUNOFF_case2.nc                  ! name of netCDF containing the HRU runoff
  !<vname_qsim>            runoff                           ! variable name of HRU runoff
  !<vname_time>            time                             ! variable name of time in the runoff file
  !<vname_hruid>           hru                              ! variable name of runoff HRU ID
  !<dname_time>            time                             ! dimension name of time
  !<dname_hruid>           hru                              ! dimension name of HM_HRU
  <units_qsim>            mm/s                             ! units of runoff
  <dt_rof>                86400                            ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE
  ! ----------------------------------
  <is_remap>              T                                          ! logical to indicate runnoff needs to be mapped to RN_HRU
  <fname_remap>           spatialweights_grid12km_nhdplus-cameo.nc   ! name of netCDF for HM_HRU-RN_HRU mapping data
  !<vname_hruid_in_remap>  polyid                                     ! variable name of RN_HRU in the mapping file
  !<vname_weight>          weight                                     ! variable name of areal weights of overlapping HM_HUs for each RN_HRU
  !<vname_qhruid>          overlapPolyId                              ! variable name of HM_HRU ID
  !<vname_num_qhru>        overlaps                                   ! variable name of numbers of HM_HRUs for each RN_HRU
  !<dname_hru_remap>       polyid                                     ! dimension name of RN_HRU (in the mapping file)
  !<dname_data_remap>      data                                       ! dimension name of ragged HM_HR
  ! **************************************************************************************************************************
  ! Namelist file name
  ! ---------------------------
  <param_nml>             param.nml.default                ! spatially constant model parameters
  ! **************************************************************************************************************************

Option 3 - runoff input is given at grid
Runoff and mapping netCDF variables are shown just to show how netCDFs are structured, and commented out in this example, since they are the same as the default names.

::

  ! *************************************************************************************************************************
  ! ***** DEFINITION OF MODEL CONTROL INFORMATION ***************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! Note: lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.
  !    lines starting with <xxx> are read till "!"
  !    Do not leave lines without "!" at the beginning of the lines
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
  <case_name>             cameo_case3                      ! simulation name - used for output netcdf name
  <sim_start>             1950-01-01 00:00:00              ! time of simulation start. year-month-day (hh:mm:ss)
  <sim_end>               1950-12-31 00:00:00              ! time of simulation end.   year-month-day (hh:mm:ss)
  <dt_qsim>               3600                             ! simulation time interval [sec]
  <newFileFrequency>      monthly                          ! history file frequency - daily, monthly, yearly or single
  <outputFrequency>       daily                            ! time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly
  ! *************************************************************************************************************************
  ! ROUTE options
  ! --------------------------
  <route_opt>             012345                           ! option for routing schemes 0-> SUM, 1->IRF, 2->KWT, 3->KW, 4->MC, 5->DW,  otherwise error
  <doesBasinRoute>        1                                ! basin routing?  0-> no, 1-> basin UH
  ! **************************************************************************************************************************
  ! DEFINE RIVER DATA
  ! ---------------------------------------
  <fname_ntopOld>         ntopo_nhdplus_cameo_pfaf.nc      ! name of netCDF containing river segment data
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF FILE
  ! ----------------------------------
  <fname_qsim>            RUNOFF_case3.nc                  ! name of netCDF containing the HRU runoff
  !<vname_qsim>            runoff                           ! variable name of HRU runoff
  !<vname_time>            time                             ! variable name of time in the runoff file
  !<dname_time>            time                             ! dimension name of time
  !<dname_xlon>            lon                              ! dimension name of x(j)
  !<dname_ylat>            lat                              ! dimension name of y(i)
  <units_qsim>            mm/s                             ! units of runoff
  <dt_rof>                86400                            ! time interval of the runoff
  ! **************************************************************************************************************************
  ! DEFINE RUNOFF MAPPING FILE
  ! ----------------------------------
  <is_remap>              T                                          ! logical to indicate runnoff needs to be mapped to RN_HRU
  <fname_remap>           spatialweights_grid12km_nhdplus-cameo.nc   ! name of netCDF for HM_HRU-RN_HRU mapping data
  !<vname_hruid_in_remap>  polyid                                     ! variable name of RN_HRU in the mapping file
  !<vname_weight>          weight                                     ! variable name of areal weights of overlapping HM_HUs for each RN_HRU
  !<vname_i_index>         i_index                                    ! variable name of ylat index
  !<vname_j_index>         j_index                                    ! variable name of xlon index
  !<vname_num_qhru>        overlaps                                   ! variable name of numbers of HM_HRUs for each RN_HRU
  !<dname_hru_remap>       polyid                                     ! dimension name of RN_HRU (in the mapping file)
  !<dname_data_remap>      data                                       ! dimension name of ragged HM_HRU
  ! **************************************************************************************************************************
  ! Namelist file name
  ! ---------------------------
  <param_nml>             param.nml.default                ! spatially constant model parameters
  ! **************************************************************************************************************************
