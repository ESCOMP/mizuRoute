.. _Output_files:

=================
Output Files
=================

.. River_Discharge

River Discharge and Lake and Reservoir Volume
---------------------------------------------

The following variables, besides time, basinID (RN_HRU ID) and reachID can be output in netCDF. Users can control which variables are output by setting <variable_name> to T or F in control file. All the variables are set to T by default.
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
| <KWTroutedRunoff>      | outflow [m3/s] from reach based on Kinematic wave tracking (KWT) reach routing. See note 3     |
+------------------------+------------------------------------------------------------------------------------------------+
| <IRFroutedRunoff>      | outflow [m3/s] from reach based on IRF reach routing. See note 3                               |
+------------------------+------------------------------------------------------------------------------------------------+
| <KWroutedRunoff>       | outflow [m3/s] from reach based on KW (Kinematic Wave) reach routing. See note 3               |
+------------------------+------------------------------------------------------------------------------------------------+
| <MCroutedRunoff>       | outflow [m3/s] from reach based on MC (Muskingum-Cunge) reach routing. See note 3              |
+------------------------+------------------------------------------------------------------------------------------------+
| <DWroutedRunoff>       | outflow [m3/s] from reach based on DW (Diffusive wave) reach routing. See note 3               |
+------------------------+------------------------------------------------------------------------------------------------+
| <KWTvolume>            | volume [m3] in reach based on Kinematic wave tracking (KWT) reach routing. See note 3          |
+------------------------+------------------------------------------------------------------------------------------------+
| <IRFvolume>            | volume [m3] in reach based on IRF reach routing. See note 3                                    |
+------------------------+------------------------------------------------------------------------------------------------+
| <KWvolume>             | volume [m3] in reach based on KW (Kinematic Wave) reach routing. See note 3                    |
+------------------------+------------------------------------------------------------------------------------------------+
| <MCvolume>             | volume [m3] in reach based on MC (Muskingum-Cunge) reach routing. See note 3                   |
+------------------------+------------------------------------------------------------------------------------------------+
| <DWvolume>             | volume [m3] in reach based on DW (Diffusive wave) reach routing. See note 3                    |
+------------------------+------------------------------------------------------------------------------------------------+
| <outputInflow>         | T -> output inflow [m3/s] to a reach for all the active routing methods                        |
+------------------------+------------------------------------------------------------------------------------------------+

1. The unit of runoff depth is the same as the unit used in runoff data.

2. If runoff depth from runoff data is already delayed by hill-slope routing outside mizuRoute, <doesBasinRoute> should be set to 0. In this case, runoff volume computed from basRunoff is populated in <dlayRunoff> and <instRunoff> is not output.

3. routed runoff corresponding to the scheme is not ouput if users deactivate a particular routing scheme with <route_opt>.

.. Augmented_River_Network_Topology

Augmented River Network Topology
--------------------------------

+------------------------+------------------------+--------------------------------------------------------------------------+
| control variable       | Default values         | Description                                                              |
+========================+========================+==========================================================================+
| <seg_outlet>           | -9999                  | outlet reach ID for subsetted river basin. See note 1                    |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <ntopAugmentMode>      | F                      | logical to indicate river network augmention mode. See note 2.           |
+------------------------+------------------------+--------------------------------------------------------------------------+
| <fname_ntopNew>        | <fname_ntopOld>_new.nc | output netCDF name for subsetted or augmented river network.             |
+------------------------+------------------------+--------------------------------------------------------------------------+

#. **River network subset mode:** If <seg_outlet> is given, network topology and parameters read from <fname_ntopOld> for just upstream of <seg_outlet> are written in <fname_ntopNew> and the program stops.

#. **River network augmentation mode:** All the computed river network topology and parameters are written in <fname_ntopNew> and the program stops.

.. _Restart_file_output:

Restart File
------------

mizuRoute does not write restart netCDF as default. The following control variables are used to control restart dropoff timing and use restart file for continuous run from the previous simulations.
The restart file is written at previous time step to the specified time. In other words, if ``Specified`` is used for <restart_write> and ``1981-01-01-00000`` is specified in <restart_date>, mizuRoute writes restart file
at ``1980-12-31 00:00:00`` for daily time step. The restart file name uses the time stamp at user specified timing. ``yearly``, ``monthly``, ``daily`` options also follow this convention.

The restart file name convension:  <case_name>.r.yyyy-mm-dd-sssss.nc

+---------------------+---------------------------------------------------------------------------------------------------------+
| control variable    | Description                                                                                             |
+=====================+=========================================================================================================+
| <restart_write>     | restart ouput options. never (default), last, specified, yearly, monthly, daily.                        |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <restart_dir>       | directory for restart files. defualt is <output_dir>                                                    |
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

