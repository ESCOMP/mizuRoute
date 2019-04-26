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


The following variables need to be defined in control file.

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
| <seg_outlet>           | outlet reach ID for subsetted river basin. See note 2                                     |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_output>         | output netCDF name for model simulation results                                           |
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_state_in>       | input restart netCDF name.                                                                | 
+------------------------+-------------------------------------------------------------------------------------------+
| <fname_state_out>      | output restart netCDF name.                                                               |
+------------------------+-------------------------------------------------------------------------------------------+

1. if <ntopWriteOption> is T or <seg_outlet> is not -999) 

2. 
