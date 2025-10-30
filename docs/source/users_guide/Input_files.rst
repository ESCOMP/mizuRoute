.. _Input_files:

=================
Input files
=================

mizuRoute expects 2 or 3 input data depending on how runoff data is provided for river routing (more needed if activating lake, water management, solute transport, gauge data).
This section describe input data required for river routing without those advanced features.

If runoff data is provided at each river network HRU (RN_HRU), river network data and runoff data are expected.
Otherwise, mizuRoute needs to remap runoff at hydrologic model HRU (HM_HRU) to river network HRU with areal weight averaging.
In this case, one additional data, remapping data, is required. All the data need to be stored in netCDF.

Basic netCDF requirement (variable, dimension, etc) are discussed below.
Dimension and variable names list mizuRoute default name but can be whatever.
If they are not default name, the variable names need to be speficied in :doc:`control file <Control_file>`.
Some of variables and dimensions (even if they are the same as default name) have to be specified in :doc:`control file <Control_file>`

.. _River_network_data:

River network file (required)
-----------------------------

River network netCDF holds river reach-reach topology, reach-hru topology, and river and hru physical parameters. The tables below list minimum requirement.
Full list of reach/hru physical parameters possibly included are :doc:`full list of river and hru physical parameters <Riv>`.

It is recommended that river network topology is built within mizuRoute instead of computing outside, while physically parameters are ideally provided per reach and hru.

Dimensions required

+------------+-----------------------------------------------------------+
| Dimension  | Description                                               |
+============+===========================================================+
| seg        | river reach                                               |
+------------+-----------------------------------------------------------+
| hru        | river network catchment or hru (hydrologic response unit) |
+------------+-----------------------------------------------------------+

Minimum variables required

+------------+------------+-----------+-------+---------------------------------------------+
| Variable   | Dimension  | Unit      | Type  | Description                                 |
+============+============+===========+=======+=============================================+
| segId      | seg        | ``-``     | int   | unique id of each stream segment            |
+------------+------------+-----------+-------+---------------------------------------------+
| HRUid      | hru        | ``-``     | int   | unique hru ID                               |
+------------+------------+-----------+-------+---------------------------------------------+
| downSegId  | seg        | ``-``     | int   | id of the downstream segment                |
+------------+------------+-----------+-------+---------------------------------------------+
| hruSegId   | hru        | ``-``     | int   | id of the stream segment the HRU flows into |
+------------+------------+-----------+-------+---------------------------------------------+
| area       | hru        | m2        | real  | hru area                                    |
+------------+------------+-----------+-------+---------------------------------------------+
| slope      | seg        | ``-``     | real  | slope of segment (elevation drop/length)    |
+------------+------------+-----------+-------+---------------------------------------------+
| length     | seg        | m         | real  | length of segment                           |
+------------+------------+-----------+-------+---------------------------------------------+

Negative or zero (<=0) value for downSegId and hruSegId is reserved for no downstream reach, meaning that such reach or hru does not flow into any reach.
(i.e., basin outlet). For this reason, segID is required to use positive integer value (>0).

.. _Runoff_data:

Runoff file(s) (required)
-------------------------

Runoff (total runoff) data can be provided as 1) 2D [time, RN_hru], 2) 2D [time, HM_hru] or 3) 3D [time, i, j].

* Option 1. runoff is given at each river network HRU
* Option 2. runoff is given at each hydrologic model HRU (non-grid)
* Option 3. runoff is given at grid

Dimensions

+--------+-----------+---------------------------------------------+
| Option | Dimension | Description                                 |
+========+===========+=============================================+
| 1,2,3  | time      | time dimension                              |
+--------+-----------+---------------------------------------------+
| 1      | RN_HRU    | river network catchment or HRU dimension    |
+--------+-----------+---------------------------------------------+
|   2    | HM_HRU    | hydrologic model catchment or HRU dimension |
+--------+-----------+---------------------------------------------+
|     3  | i         | x direction dimension                       |
+        +-----------+---------------------------------------------+
|        | j         | y direction dimension                       |
+--------+-----------+---------------------------------------------+

Variables

+--------+-----------+--------------+--------------------------------------+-------+-------------------------+
| Option | Variable  | Dimension    | Unit                                 | Type  | Description             |
+========+===========+==============+======================================+=======+=========================+
| 1,2,3  | time      | time         | [time-unit] since yyy-mm-dd 00:00:00 | real  | time                    |
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+
| 1      | RN_hruID  | RN_hru       | ``-``                                | int   | river network HRU ID    |
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+
|   2    | HM_hruID  | HM_hru       | ``-``                                | int   | hydrologic model HRU ID |
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+
| 1      | runoff    | time, RN_hru | [length-unit]/[time-unit]            | real  | total runoff            |
+--------+           +--------------+                                      +       +                         +
|   2    |           | time, HM_hru |                                      |       |                         |
+--------+           +--------------+                                      +       +                         +
|     3  |           | time, i, j   |                                      |       |                         |
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+

Attributes: Time variable need at least 2 attributes- *units* and *calendar*. Four types of calendar can be handled. These are noleap, standard, gregorian, and proleptic_gregorian.
Time unit format is shown in the table.

.. _Runoff_mapping_data:

Runoff mapping file (required depending on case)
------------------------------------------------

mizuRoute has a capability to remap forcing at different catchments or grid to catchment or grid defined in river network used for routing using weighted average. A user needs to provide a mapping file in netCDF.
See :ref:`Runoff mapping data <Runoff_mapping_data>` for mapping file structure.
Breifly, mapping can be either catchment (i.e., unstructure grid) to river network catchment (option 2) or grid to river network catchment (option 3). option 1 is forcing provided at the same catchment as the one in river network, in which case no mapping is required.

For runoff input options 2 and 3, runoff mapping data in netCDF must be provided so that weighted average runoff value is computed for each river network HRU.

+--------+-----------+---------------------------------------------+
| Option | Dimension | Description                                 |
+========+===========+=============================================+
| 2,3    | hru       | River network HRU                           |
+--------+-----------+---------------------------------------------+
| 2,3    | data      | Vectorized overlapping HRU (or grid boxes)  |
+--------+-----------+---------------------------------------------+

Control file keys for remapping.

+--------+------------------------+----------------------------------------------------------------------------------------------------+
| option | control variable       | Descriptions                                                                                       |
+========+========================+====================================================================================================+
|        | <is_remap>             | Logical to indicate runoff needs to be remapped to RN_HRU. set T to activate remapping option      |
+--------+------------------------+----------------------------------------------------------------------------------------------------+
|   2,3  | <fname_remap>          | netCDF name of runoff remapping                                                                    |
+--------+------------------------+----------------------------------------------------------------------------------------------------+

Required runoff mapping netCDF variables

+--------+------------------------+----------------------------------------------------------------------------------------------------+
| option | control variable       | Descriptions                                                                                       |
+========+========================+====================================================================================================+
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


+--------+------------+-----------+-------+-------+-----------------------------------------------+
| Option | Variable   | Dimension | Unit  | type  | Descriptions                                  |
+========+============+===========+=======+=======+===============================================+
| 2,3    | RN_hruId   | hru       | ``-`` | int   | River network HRU ID                          |
+--------+------------+-----------+-------+-------+-----------------------------------------------+
| 2,3    | nOverlaps  | hru       | ``-`` | int   | number of overlapping HM_HRUs for each RN_HRU |
+--------+------------+-----------+-------+-------+-----------------------------------------------+
| 2,3    | weihgt     | data      | ``-`` | real  | areal weight of overlapping HM_HRUs           |
+--------+------------+-----------+-------+-------+-----------------------------------------------+
| 2      | HM_hruId   | data      | ``-`` | int   | ID of overlapping HM_HRUs                     |
+--------+------------+-----------+-------+-------+-----------------------------------------------+
|   3    | i_index    | data      | ``-`` | int   | i direction index overlapping grid boxes      |
+        +------------+-----------+-------+-------+-----------------------------------------------+
|        | j_index    | data      | ``-`` | int   | j direction index overlapping grid boxes      |
+--------+------------+-----------+-------+-------+-----------------------------------------------+

Creating a mapping is basically GIS intersection of two geometries. The figure below visualizes intersection between runoff grid (option 3) and river network catchment (HRU) polygons.
This example (right bottom) shows river network HRU, c\ :sub:`k`\, has 11 overlapping grid boxes (nOverlaps=11 in a mapping netCDF. see table above) and corresponding weights (i.e., fractions of each overlapped grid boxes to total area of c\ :sub:`k`\) as well as i_index and j_index.
In a mapping netCDF, all 1D arrays of weights (and i_index and j_index) from each HRU are combined for a large single 1D array. The order of the arrays from each HRU must match the order of RN_hruId

.. image:: images/mapping_schematic.png
  :width: 600

There are a few tools available to create the netCDF with required data:

#. mizuRoute_remapping (https://github.com/ShervanGharari/mizuRoute_remapping)


.. _UnifiedASCII_file:

Unified ASCII parameter file (optional)
---------------------------------------


.. _Restart_file_input:

Restart file (optional)
-----------------------


.. _WaterManagement_file:

Water management file (optional)
--------------------------------



.. _GaugeData_file:

Gauge data file (optional)
--------------------------

mizuRoute can read gauge observed discharge data (in netCDF) along with gauge meta ascii data. To read gauge observation and gauge metadata, the following control variables need to be specified.
gauge meta ascii file is csv format, and  should include at least gauge id and corresponding reach id
gauge discharge data is used for data assimilation (current version does not include this at this moment)
Using gauge data, a user can output the simulation at gauge only output in addition to at the entire river network and/or direct insertion to modify discharge whenever observed discharge is available.

+---------------------+---------------------------------------------------------------------------------------------------------+
| control variable    | Description                                                                                             |
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

Direct insertion, the simplest data assimilation, can be  performed at a list of reaches in the metadata. Two parameters, <QerrTrend> and <ntsQmodStop>, are needed.
<QerrTrend> tells how bias computed at observation time at each reach evolves in the subsequent future <ntsQmodStop> time steps.
To activate direct insertion of observed discharge into simulated discharge, the following control variables need to be specified.

+---------------------+---------------------------------------------------------------------------------------------------------+
| control variable    | Description                                                                                             |
+=====================+=========================================================================================================+
| <qmodOption>        | activation of direct insertion. 0 -> do nothing, 1=> discharge direct insertion                         |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <ntsQmodStop>       | the number of time steps when flow correction stops                                                     |
+---------------------+---------------------------------------------------------------------------------------------------------+
| <QerrTrend>         | temporal discharge error trend. 1->constant, 2->linear, 3->logistic, 4->exponential                     |
+---------------------+---------------------------------------------------------------------------------------------------------+
