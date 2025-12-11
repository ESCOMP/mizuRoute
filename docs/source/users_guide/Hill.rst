.. _Hillslope_Routing:

Hillslope routing options
====================================

Runoff from hydrologic model can be areal mean depth over 1) hydrologic model HRU (or grid box) or 2) river network HRU.
The former case requires runoff remapping from hydrolgic model HRU to river network model HRU. See :ref:`remapping data <Runoff_mapping_data>` for how to configure remapping in the control file.

Once runoff depth is in river network HRU space, mizuRoute converts runoff depth to the volume by multiplying the depth by river network model HRU area.
The runoff volume at a river network HRU is input into the river reach that corresponds to the HRU.

However, there is a option to *delay* timing of runoff to account for the travel time of runoff to the river reach (hillslope routing).
MizuRoute uses unit-hydrograph formulated with gamma distribution with two parameter to accomplish this computation (though there are more methods that could be implemented)

To enable the unit-hydrograph based hillslope routing, a user is required to include:

.. list-table:: Control key to enable hillslope routing
   :header-rows: 1
   :widths: 20 15 15 50
   :name: channel-parameter-specification

   * - Control key
     - Type
     - Default
     - Description
   * - ``<doesBasinRoute>``
     - int
     - 1
     - hillslope routing options. 0-> no (already routed), 1->gamma distribution UH

If runoff from hydrologic model is already delayed with their hillslope routing, user should set it to 0 (no hillslope routing) to avoid duplicated computations.

There are two gamma distribution parameters (shape parameter: ``fshape`` and scale parameter: ``tscale``) needed, and currently they are provided as spatial constant parameter from namelist.
Please see :ref:`Spatially-constant parameter namelist <namelist_file>` for example setup

For more about gamma distribution based hillslop methods, please see :ref:`Hillslope routing method <Hillslope_routing_scheme>`.
