River parameters
================

Full list of river parameters, both physical and topological ones, can be output in netCDF as river network augmentation mode. 
Those augmented variables can be read in from augmented network netCDF and variable names need to be specified in :doc:`control file <Control_file>` 

To read additional augmented network parameters, <hydGeometryOption> and <topoNetworkOption> needs to be turned on (specified as 0) in :doc:`control file <Control_file>`  

Names of the river network variables (both network topology and physical parameters) can be also speficied in :doc:`control file <Control_file>`,
if they are different than their default names. The format is 

<varname_PARAMETER_DEFAULT_NAME>   NEW_NAME    ! 


Dimensions

+------------+-----------------------------------------------------------+
| Dimension  | Description                                               |
+============+===========================================================+
| seg        | river reach                                               | 
+------------+-----------------------------------------------------------+
| hru        | river network catchment or hru (hydrologic response unit) | 
+------------+-----------------------------------------------------------+
| upSeg      | immediate upstream reaches                                | 
+------------+-----------------------------------------------------------+
| upHRU      | HRUs contributing to a reach                              | 
+------------+-----------------------------------------------------------+
| upAll      | all the upstream reaches                                  | 
+------------+-----------------------------------------------------------+

physical parameters
*******************

+---------------+------------+-----------+-------+-------------------------------------------------------+
| Variable      | Dimension  | Unit      | Type  | Description                                           |
+===============+============+===========+=======+=======================================================+
| width         | seg        | ``-``     | real  | channel width                                         |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| man_n         | seg        | ``-``     | real  | mannings n                                            |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| hruArea       | upHRU      | m2        | real  | area of each contributing HRU                         |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| weight        | upHRU      | ``-``     | real  | weight assigned to each HRU                           |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| basArea       | seg        | m2        | real  | total area of contributing HRUs                       |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| upsArea       | seg        | m2        | real  | area above the top of the reach. 0 if headwater       |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| totalArea     | seg        | m2        | real  | area above the bottom of the reach (bas + ups)        |
+---------------+------------+-----------+-------+-------------------------------------------------------+
| timeDelayHist | uh         | sec       | real  | time delay histogram for each reach (only UH routing) |
+---------------+------------+-----------+-------+-------------------------------------------------------+

Topology parameters
*******************

Extra or augmented river reach and hru topology are typically computed internally. It is recommended to compute instead of generating outside mizuRoute

Variables

+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| Variable        | Dimension  | Unit      | Type  | Description                                                    |
+=================+============+===========+=======+================================================================+
| segIndex        | seg        | ``-``     | int   | reach Index                                                    |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| downSegId       | seg        | ``-``     | int   | downstream reach ID                                            |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| downSegIndex    | seg        | ``-``     | int   | downstream reach index                                         |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| upSegIds        | upSeg      | ``-``     | int   | Immediate upstream reach IDs for each reach                    |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| upSegIndices    | upSeg      | ``-``     | int   | immediate upstream reach indices for each reach                |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| allUpSegIndices | upAll      | ``-``     | int   | all the upstream reach indices for each reach                  |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| rchOrder        | seg        | ``-``     | int   | routing processing order                                       |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| goodBasin       | upSeg      | ``-``     | int   | flag to indicate immediate upstream HRUs are good HRU (area>0) |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| HRUindex        | hur        | ``-``     | int   | RN_HRU index                                                   |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| hruSegIndex     | hur        | ``-``     | int   | index of the reach below each HRU                              |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| hruContribIx    | upHRU      | ``-``     | int   | indices of HRUs contributing flow to each reach                |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+
| hruContribId    | upHRU      | ``-``     | int   | IDs of HRUs contributing flow to each reach                    |
+-----------------+------------+-----------+-------+----------------------------------------------------------------+


