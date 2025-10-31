.. _River_routing_config:

River or In-Channel Routing Models
==================================

MizuRoute has five channel routing schemes:

1. Impulse response function (IRF)
2. Lagrangian kinematic wave (LKW)
3. Euler kinematic wave (KW)
4. Muskingum-Cunge (MC)
5. Diffusive Wave (DW)

Please see :ref:`River routing schemes <River_routing_schemes>` for more information on methods and numerical solutions of each method.
Additionally, mizuRoute has an option to compute accumulated runoff [m3/s] without any channel routing at any outlets of river reaches, if a user wants to get total runoff volume of upstream areas.

Currently, a user must provide the control key ``<routOpt>`` to select which routing option(s) to be used. ``<routOpt>`` is an integer that corresponds to the routing schemes (IRF=1, LKW=2, KW=3, MC=4, DW=5). Runoff accumulation option uses ``0``.
A user can run multiple routing schemes at the same time by setting ``<routOpt>`` to multiple IDs. For example, if a user wants to use runoff accumulation, KW, and DW, ``<routOpt>`` should be ``035``.


.. _channel_physical_parameters:

Channel physical parameters
---------------------------

IRF requires two parameters (wave celerity and diffusivity) along with **slope** and **length**, used to compute unit-hydrograph, while the other schemes uses a few channel physical parameters in addition to **slope** and **length**.
Regardless of the routing schemes, **slope** and **length** are the minimum channel physical parameters that must be included in a river input data See :ref:`River input data <River_network_data>`.
If users uses minimum set of channel physical parameters, mizuRoute will compute all the channel physical parameters using spatially constant parameter specified in the namelist.
To setup spatially distributed river channel parameters, a user can provide channel physical parameters at each river reach in river input data.

To use the channel parameters from the netCDF, make sure that the key ``<hydGeometryOption>`` in the control file **must be set to 0**. The default values is ``1`` (compute channel parameter internally).

As a default, a river channel does not have floodplains, meaning river water is always contained in a channel.
A user can add floodplain by adding the control key ``<floodplain>`` with ``T``. By adding floodplain, discharge tends to be attenuated due to greater water-riverbed contact area.
For active floodplain option, the channel bankfull depth needs to be computed as a default, or supplied in the river input netCDF.
Also, by adding floodplain, water storage over floodplain is computed in addition to total water storage, which may be used for furter flood mapping (outside mizuRoute) or feedback to the land model (for CESM coupled mode)
Note that currently floodplain is activate only for KW, MC, and DW routing schemes.

How the physical parameters control the shape of the channel cross-section of the channel without and with floodplain is depicted below:

.. image:: images/channel_xsection_no_fp.png
  :width: 600

.. image:: images/channel_xsection_fp.png
  :width: 600

Below is the control keys related to the option of channel geometry specification.

.. list-table:: Control key to method channel physical parameter specification
   :header-rows: 1
   :widths: 20 15 15 50
   :name: channel-parameter-specification

   * - Control key
     - Type
     - Default
     - Description
   * - ``<hydGeometryOption>``
     - int
     - ``1``
     - Options for channel physical parameters estimations

       * ``0`` → read from a river input data
       * ``1`` → compute internally
   * - ``<floodplain>``
     - logical
     - ``F``
     - Options to add a simple floodplain

       * ``F`` → no floodplain, channel is unlimited bank depth
       * ``T`` → add floodplain, and floodwater volume is computed


Then, a user needs to specify the variable name related channel physical properties in the control file.

.. list-table:: channel routing related control keys in the river input file
   :widths: 20 20 15 15 15 15 30
   :header-rows: 1
   :name: channel-parameter-variables

   * - Control key
     - Type
     - Variable type
     - Variable dimension
     - Variable unit
     - routing schemes
     - Description
   * - ``<varname_width>``
     - NetCDF variable name
     - real
     - seg
     - m
     - KWT, KW, MC, DW
     - channel bottom width
   * - ``<varname_man_n>``
     - NetCDF variable name
     - real
     - seg
     - \-
     - KWT, KW, MC, DW
     - manning n coefficient
   * - ``<varname_sideSlope>``
     - NetCDF variable name
     - real
     - seg
     - \-
     - KWT, KW, MC, DW
     - channel side slope. vertical:horisontal=1:sideSlope
   * - ``<varname_depth>``
     - NetCDF variable name
     - real
     - seg
     - m
     - KWT, KW, MC, DW
     - channel bankful depth for active floodplain


.. _Miscleneous_control_keys:

Miscleneous control keys
------------------------

There are few miscleneous control keys available for a channel routing. Note that these are not implemented to all the routing schemes.

.. list-table:: Miscleneous channel routing related control keys
   :widths: 20 15 15 15 50
   :header-rows: 1
   :name: Miscleneous channel routing control keys

   * - Control key
     - Type
     - Default
     - routing schemes
     - Description
   * - ``<min_length_route>``
     - real
     - 0.0
     - IRF, KW, MC, DW
     - minimum reach length [m] for routing to be performed. pass-through (outflow=infolw+local flow) is performed for length less than this threshold.
   * - ``<hw_drain_point>``
     - int
     - 2
     - IRF, KW, MC, DW
     - how to add local runoff in a reach of headwater HRUs.
       * ``1`` → top of reach
       * ``2`` → bottom of reach (default)


.. Full list of river parameters, both physical and topological ones, can be output in netCDF as river network augmentation mode.
.. Those augmented variables can be read in from augmented network netCDF and variable names need to be specified in :doc:`control file <control_file>`

.. To read additional augmented network parameters, <hydGeometryOption> and <topoNetworkOption> needs to be turned on (specified as 0) in :doc:`control file <control_file>`

.. Names of the river network variables (both network topology and physical parameters) can be also speficied in :doc:`control file <control_file>`,
.. if they are different than their default names. The format is

.. <varname_PARAMETER_DEFAULT_NAME>   NEW_NAME    !


.. Dimensions

.. +------------+-----------------------------------------------------------+
.. | Dimension  | Description                                               |
.. +============+===========================================================+
.. | seg        | river reach                                               |
.. +------------+-----------------------------------------------------------+
.. | hru        | river network catchment or hru (hydrologic response unit) |
.. +------------+-----------------------------------------------------------+
.. | upSeg      | immediate upstream reaches                                |
.. +------------+-----------------------------------------------------------+
.. | upHRU      | HRUs contributing to a reach                              |
.. +------------+-----------------------------------------------------------+
.. | upAll      | all the upstream reaches                                  |
.. +------------+-----------------------------------------------------------+

.. .. _physical_parameters:

.. physical parameters
.. *******************

.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | Variable      | Dimension  | Unit      | Type  | Description                                           |
.. +===============+============+===========+=======+=======================================================+
.. | width         | seg        | ``-``     | real  | channel width                                         |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | man_n         | seg        | ``-``     | real  | mannings n                                            |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | hruArea       | upHRU      | m2        | real  | area of each contributing HRU                         |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | weight        | upHRU      | ``-``     | real  | weight assigned to each HRU                           |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | basArea       | seg        | m2        | real  | total area of contributing HRUs                       |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | upsArea       | seg        | m2        | real  | area above the top of the reach. 0 if headwater       |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | totalArea     | seg        | m2        | real  | area above the bottom of the reach (bas + ups)        |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+
.. | timeDelayHist | uh         | sec       | real  | time delay histogram for each reach (only UH routing) |
.. +---------------+------------+-----------+-------+-------------------------------------------------------+

.. .. _Topology_parameters:

.. Topology parameters
.. *******************

.. Extra or augmented river reach and hru topology are typically computed internally. It is recommended to compute instead of generating outside mizuRoute

.. Variables

.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | Variable        | Dimension  | Unit      | Type  | Description                                                    |
.. +=================+============+===========+=======+================================================================+
.. | segIndex        | seg        | ``-``     | int   | reach Index                                                    |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | downSegId       | seg        | ``-``     | int   | downstream reach ID                                            |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | downSegIndex    | seg        | ``-``     | int   | downstream reach index                                         |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | upSegIds        | upSeg      | ``-``     | int   | Immediate upstream reach IDs for each reach                    |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | upSegIndices    | upSeg      | ``-``     | int   | immediate upstream reach indices for each reach                |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | allUpSegIndices | upAll      | ``-``     | int   | all the upstream reach indices for each reach                  |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | rchOrder        | seg        | ``-``     | int   | routing processing order                                       |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | goodBasin       | upSeg      | ``-``     | int   | flag to indicate immediate upstream HRUs are good HRU (area>0) |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | HRUindex        | hur        | ``-``     | int   | RN_HRU index                                                   |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | hruSegIndex     | hur        | ``-``     | int   | index of the reach below each HRU                              |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | hruContribIx    | upHRU      | ``-``     | int   | indices of HRUs contributing flow to each reach                |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+
.. | hruContribId    | upHRU      | ``-``     | int   | IDs of HRUs contributing flow to each reach                    |
.. +-----------------+------------+-----------+-------+----------------------------------------------------------------+



.. Impulse response function
.. --------------------------
..
.. Lagrangian kinematic wave
.. -------------------------
..
.. Euler kinematic wave
.. ---------------------
..
.. Muskingum–Cunge
.. ----------------
..
.. Diffusive wave
.. ---------------
