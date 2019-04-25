=================
input data
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

+------------+------------+-----------+-----------------------------------------+
| variables  | dimension  | units     | descriptions                            |
+============+============+===========+=========================================+
| segId      | seg        | ``-``     | unique id of each stream segment        |
+------------+------------+-----------+-----------------------------------------+
| HRUid      | hru        | ``-``     | unique hru ID                           |
+------------+------------+-----------+-----------------------------------------+
| hruSegId   | hru        | ``-``     | id of the stream segment below each HRU |
+------------+------------+-----------+-----------------------------------------+
| area       | hru        | m2        | hru area                                |
+------------+------------+-----------+-----------------------------------------+
| slope      | seg        | ``-``     | slope of segment                        |
+------------+------------+-----------+-----------------------------------------+
| length     | seg        | m         | length of segment                       |
+------------+------------+-----------+-----------------------------------------+
