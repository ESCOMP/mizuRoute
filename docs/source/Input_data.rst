=================
Input data
=================

mizuRoute expects 2 or 3 input data depending on how runoff data is provided. 
If runoff data is provided at each river network HRU (RN_HRU), river network data and runoff data are expected.
Otherwise, mizuRoute needs to remap runoff at hydrologic model HRU (HM_HRU) to river network HRU with areal weight averaging. 
In this case, one additional data, remapping data, is required. All the data need to be stored in netCDF.

Basic netCDF requirement (variable, dimension, etc) are discussed below.
Dimension and variable names list mizuRoute default name but can be whatever. 
If they are not default name, the variable names need to be speficied in :doc:`control file <Control_file>`.
Some of variables and dimensions (even if they are the same as default name) have to be specified in :doc:`control file <Control_file>`

River network data
------------------

River network netCDF holds river reach-reach topology, reach-hru topology, and river and hru physical parameters. The tables below list minimum requirement.
Full list of reach/hru physical parameters possibly included are :doc:`full list of river and hru physical parameters <seg_hru_param>`. 

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

+------------+------------+-----------+-------+-----------------------------------------+
| Variable   | Dimension  | Unit      | Type  | Description                             |
+============+============+===========+=======+=========================================+
| segId      | seg        | ``-``     | int   | unique id of each stream segment        |
+------------+------------+-----------+-------+-----------------------------------------+
| HRUid      | hru        | ``-``     | int   | unique hru ID                           |
+------------+------------+-----------+-------+-----------------------------------------+
| downSegId  | seg        | ``-``     | int   | id of the downstream segment            |
+------------+------------+-----------+-------+-----------------------------------------+
| hruSegId   | hru        | ``-``     | int   | id of the stream segment below each HRU |
+------------+------------+-----------+-------+-----------------------------------------+
| area       | hru        | m2        | real  | hru area                                |
+------------+------------+-----------+-------+-----------------------------------------+
| slope      | seg        | ``-``     | real  | slope of segment                        |
+------------+------------+-----------+-------+-----------------------------------------+
| length     | seg        | m         | real  | length of segment                       |
+------------+------------+-----------+-------+-----------------------------------------+

Negative or zero (<=0) value for downSegId is reserved for no downstream reach, meaning that such reach or hru does not flow into any reach.
(i.e., basin outlet). For this reason, segID is required to use positive integer value (>0).

Runoff data
-----------

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
| 2      | HM_HRU    | hydrologic model catchment or HRU dimension | 
+--------+-----------+---------------------------------------------+
| 3      | i         | x direction dimension                       | 
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
| 2      | HM_hruID  | HM_hru       | ``-``                                | int   | hydrologic model HRU ID | 
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+
| 1      | runoff    | time, RN_hru | [length-unit]/[time-unit]            | real  | total runoff            |
+--------+           +--------------+                                      +       +                         +
| 2      |           | time, HM_hru |                                      |       |                         |
+--------+           +--------------+                                      +       +                         +
| 3      |           | time, i, j   |                                      |       |                         |
+--------+-----------+--------------+--------------------------------------+-------+-------------------------+

Attributes: Time variable need at least 2 attributes- *units* and *calendar*. Four types of calendar can be handled. These are noleap, standard, gregorian, and proleptic_gregorian.
Time unit format is shown in the table.

Runoff mapping data
-------------------

For runoff input options 2 and 3, runoff mapping data, also in netCDF format, is necessary to compute runoff value for each river network HRU

+--------+-----------+---------------------------------------------+
| Option | Dimension | Description                                 |
+========+===========+=============================================+
| 2,3    | hru       | River network HRU                           | 
+--------+-----------+---------------------------------------------+
| 2,3    | data      | Vectorized overlapping HRU (or grid boxes)  | 
+--------+-----------+---------------------------------------------+

Required runoff mapping netCDF variables 

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
| 3      | i_index    | data      | ``-`` | int   | i(x) direction index overlapping grid boxes   |
+        +------------+-----------+-------+-------+-----------------------------------------------+
|        | j_index    | data      | ``-`` | int   | j(y) direction index overlapping grid boxes   |
+--------+------------+-----------+-------+-------+-----------------------------------------------+

Hope this picture helps you understand mapping netCDF variables.

.. image:: images/mapping1.png
  :width: 400
