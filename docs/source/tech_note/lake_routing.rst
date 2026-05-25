.. _Lake_schemes:

Lake Routing Schemes
====================

The lake component in mizuRoute provides a flexible framework for representing the water balance of both natural lakes and managed reservoirs within a river routing network.


.. _Lake_Water_Balance_Representation_scheme:

Lake Water Balance Representation
---------------------------------

The lake module simulates the water balance by accounting for the primary fluxes entering and leaving the system.

**Input fluxes:**

- Upstream river discharge
- Direct precipitation or other areal fluxes onto the lake

**Output fluxes:**

- Evaporation or other areal fluxes from the lake
- Outflow discharge (for exorheic lakes)

Lake storage evolves dynamically as a function of these fluxes, allowing the model to resolve lake and reservoir behavior over time.

The water balance of a lake or reservoir in mizuRoute can be expressed as:

.. math::
   :label: lake_water_balance

   \frac{dS}{dt} = I - O + F_{\mathrm{eff}}\,A - F_{a,i}

where:

- :math:`S` [m³] is the lake or reservoir storage
- :math:`I` and :math:`O` [m³ s⁻¹] are inflow and outflow, respectively
- :math:`A` [m²] is the lake surface area
- :math:`F_{\mathrm{eff}}` [m s⁻¹] is the effective areal forcing over the lake
- :math:`F_{a,i}` [m³ s⁻¹] is the abstraction or injection flux provided as a time series

The flux :math:`F_{a,i}` is defined such that positive values represent abstraction (water removal), provided sufficient water is available in the river segment or lake, while negative values represent injection.

The effective forcing term can be decomposed as:

.. math::
   :label: effective_areal_forcing

   F_{\mathrm{eff}} = P - E + G

where:

- :math:`P` [m s⁻¹] is precipitation over the lake surface
- :math:`E` [m s⁻¹] is evaporation from the lake surface
- :math:`G` [m s⁻¹] represents groundwater exchange with the lake

In mizuRoute, precipitation and evaporation can be explicitly represented, depending on the model configuration. Users may choose to apply precipitation and evaporation separately, or use an effective precipitation approach in which these terms are already combined externally (e.g., by a land surface model providing runoff or net water flux to the lake).

Groundwater exchange is not explicitly represented in the current implementation. However, it may be incorporated externally either as part of an effective runoff or within the effective precipitation term, depending on how the forcing data are constructed.

For further details on the formulation, see
`Gharari et al., 2024 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR032400>`_.

.. _Lake_Types_scheme:

Lake Types
----------

Two types of lakes are supported:

- **Exorheic lakes**: lakes with an outlet, allowing water to flow downstream.
- **Endorheic lakes**: closed basins with no surface outflow, where water leaves the system only through evaporation or infiltration.

.. _Lake_parametric_models_scheme:

Lake Routing Formulations
-------------------------

One commonly used approach is **parametric lake models**, which relate lake inflow and storage to outflow through functional relationships. These formulations are computationally efficient and well suited for large-scale applications. mizuRoute includes multiple formulations to represent lake dynamics, each varying in complexity and data requirements; currently, three lake models are implemented (see :ref:`Lake_res_model` for details).

.. _Lake_Target_Volume_scheme:

Target Volume
-------------

The model also supports a simplified representation of reservoir operations through a **target volume** mechanism.

The target volume is provided as a time series in the water management input file for each lake or reservoir. This approach enables approximation of reservoir rule curves and their temporal variability, incorporation of in situ or satellite-based storage observations, and the ability to inform mizuRoute using outputs from more complex reservoir operation models during events such as floods or other regulated conditions; for more details, see :ref:`Lake_Target_Volume`.

.. _Lake_Abstraction_scheme:

Water Management Fluxes
-----------------------

In addition to physically based formulations, mizuRoute supports a data-driven representation of lake and reservoir management, where users can prescribe external input and output fluxes. These fluxes are provided as time series, similar to the treatment of regulated river segments, enabling the representation of human interventions in the hydrological system; furthermore, this framework allows interactions with other hydrological processes, such as exchanges with groundwater or other connected storage components. The magnitude of abstraction is constrained by the available water in the lake or reservoir, ensuring consistency with the water balance; if storage is insufficient (e.g., at or near zero), the prescribed abstraction is reduced accordingly, and may become zero when no water is available for extraction. For details on the required input format, see :ref:`Lake_Abstraction`.