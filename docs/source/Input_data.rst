=================
Input data
=================
River network data
------------------

River network data should be stored in netCDF format.
dimension and varialbe names below use mizuRoute default name but can be whatever. 

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
| hruSegId   | hru        | ``-``     | int   | id of the stream segment below each HRU |
+------------+------------+-----------+-------+-----------------------------------------+
| area       | hru        | m2        | real  | hru area                                |
+------------+------------+-----------+-------+-----------------------------------------+
| slope      | seg        | ``-``     | real  | slope of segment                        |
+------------+------------+-----------+-------+-----------------------------------------+
| length     | seg        | m         | real  | length of segment                       |
+------------+------------+-----------+-------------------------------------------------+

Runoff data
-----------

Runoff (total runoff) data can be provided as 1) [time, RN_hru], 2) [time, HM_hru] or 3) [time, i, j].

* Option 1. runoff is given at each river network HRU 
* Option 2. runoff is given at each hydrologic HRU (non-grid) 
* Option 3. runoff is given by grid 

Dimension

+--------+-----------+---------------------------------------------+
| Option | Dimension | Description                                 |
+========+===========+=============================================+
| 1,2,3  | time      | time dimension                              | 
+--------+-----------+---------------------------------------------+
| 1      | RN_hru    | river network catchment or HRU dimension    | 
+--------+-----------+---------------------------------------------+
| 2      | HM_hru    | hydrologic model catchment or HRU dimension | 
+--------+-----------+---------------------------------------------+
| 3      | i         | x direction dimension                       | 
+        +-----------+---------------------------------------------+
|        | j         | y direction dimension                       | 
+--------+-----------+---------------------------------------------+

Variable

+-----------+----------------+-------+-------+---------------------+
| Variable  | Dimension      | Unit  | Type  | Description         |
+===========+================+=======+=======+=====================+
| runoff    | [time, RN_hru] | L/T   | real  | total runoff        |
|           | [time, HM_hru] |       |       |                     |
|           | [time, i, j]   |       |       |                     |
+-----------+----------------+-------+-------+---------------------+


Runoff mapping data
-------------------

For runoff input options 1 and 2, runoff mapping data, also in netCDF format is necessary to get runoff value at each river network HRU

+--------+-----------+---------------------------------------------+
| Option | Dimension | Description                                 |
+========+===========+=============================================+
| 2,3    | hru       | river network HRU                           | 
+--------+-----------+---------------------------------------------+
| 2,3    | data      | vectorized overlapping HRU (or grid boxes)  | 
+--------+-----------+---------------------------------------------+

minimum runoff mapping netCDF variables 

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
| 3      | i_index    | data      | ``-`` | int   | i direction index overlapping grid boxes      |
+        +------------+-----------+-------+-------+-----------------------------------------------+
|        | j_index    | data      | ``-`` | int   | j direction index overlapping grid boxes      |
+--------+------------+-----------+-------+-------+-----------------------------------------------+

