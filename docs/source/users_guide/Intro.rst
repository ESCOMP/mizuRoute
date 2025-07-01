.. _Introduction:

Introduction
============

Desired technical knowledge:

**Geographic Information System (GIS)**: River and catchment are represented by line and polygon respectively in mizuRoute. These "vector" spatial data (ESRI shapefile, Geopackage) are handled with GIS software. 
`QGIS <https://qgis.org/>`_ is a free GIS software that works on Windows, Mac and Linux. QGIS is intuitive, and convenient when interactively viewing spatial data (zoom, pan) and creating a map, and query of spatial elements.
GDAL library, and/or python packages e.g., geopanda, are also usefull and best used programically to manupulate the geometries and attributes.

**netCDF**: Though river spatial data usually originate from line and polygon vector spatial data, these vector data are eventually converted into netCDF for mizuRoute input. 
Output from mizuRoute is also netCDF file.
Basic technical knowledge on netCDF operations using netCDF commands (ncdump, nccopy etc.) and NCO and python are required to effectively work on mizuRoute input and output.


User's guide provides the following instructions on setup mizuRoute:

:doc:`1.2. Building model <Build_model>` provide stand-alone mizuRoute build instruction: external libraries requirement and step-by-step instruction on how to compile the model.

:doc:`1.3. Control file <Control_file>` provide information on how to setup a control file (configuration file), and various control variables.

:doc:`1.4. Input data <Input_data>` provide input data requirement for river routing, not including advanced features: lake, and water management.

:doc:`1.5. River model parameters <seg_hru_param>` provides information on river network topology and physical parameters.

:doc:`1.6. Lake model setup <lake>` section provide instructions on how to modify river network into river-lake network, lake model parameter specification, control variable related to lake model.
For scientific description is provided in :doc:`Lake water balance computation <../tech_note/lake_routing>`.

:doc:`1.7. Water management option <water_management>` section provides additional data inputs related to water managments, control variables related to water management.

:doc:`1.8. testCase data <testCase>` provide a link to testCase including model setup and data for river routing only.
