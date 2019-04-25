=================
Input data
=================
River network data
------------------

River network data should be stored in netCDF format.
dimension and varialbe names below use mizuRoute default name but can be whatever. 

Dimensions required

+------------+---------------------------------------------+
| dimension  | description                                 |
+============+=============================================+
| seg        | river reach                                 | 
+------------+---------------------------------------------+
| hru        | catchment or hru (hydrologic response unit) | 
+------------+---------------------------------------------+

Minimum variables required

+------------+------------+-----------+-------+-----------------------------------------+
| variables  | dimension  | units     | type  | descriptions                            |
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

runoff data
-----------


runoff mapping data
-------------------

Runoff mapping data, also in netCDF format is necessary if a hydrologic model use different hru than the one in river network data.
For example, typically hydrologic model use grid while river network use catchment. In this case, gridded runoff data need to be mapped to catchment runoff.

minimum runoff mapping netCDF variables 

+------------+------------+-----------+-------+-----------------------------------------+
| variables  | dimension  | units     | type  | descriptions                            |
+============+============+===========+=======+=========================================+
| segId      | polyid     | ``-``     | int   |                                         |
+------------+------------+-----------+-------+-----------------------------------------+
| HRUid      | polyid     | ``-``     | int   |                                         |
+------------+------------+-----------+-------+-----------------------------------------+
| hruSegId   | data       | ``-``     | int   |                                         |
+------------+------------+-----------+-------+-----------------------------------------+
| area       | data       | ``-``     | real  |                                         |
+------------+------------+-----------+-------+-----------------------------------------+
| slope      | seg        | ``-``     | real  |                                         |
+------------+------------+-----------+-------+-----------------------------------------+
| length     | seg        | ``-``     | real  |                                         |
+------------+------------+-----------+-------------------------------------------------+
