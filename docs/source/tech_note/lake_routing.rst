.. _Lake_schemes:

Lake Routing Schemes
====================

The lake component in mizuRoute provides a flexible framework for representing
both natural lakes and managed reservoirs within a river routing network.

Lake Types
----------

Two types of lakes are supported:

- **Exorheic lakes**: Lakes with an outlet, allowing water to flow downstream.
- **Endorheic lakes**: Closed basins with no surface outflow, where water
  leaves the system only through evaporation or infiltration.

Water Balance Representation
----------------------------

The lake module simulates the water balance by accounting for the primary
fluxes entering and leaving the lake system.

**Input fluxes:**
- Upstream river discharge
- Direct precipitation onto the lake surface

**Output fluxes:**
- Evaporation from the lake surface
- Outflow discharge (for exorheic lakes)

Lake storage evolves dynamically as a function of these fluxes, enabling the
model to resolve lake behavior over time.

Lake Routing Formulations
-------------------------

mizuRoute includes multiple formulations to represent lake dynamics, each
varying in complexity and data requirements. Currently, three lake models are
implemented (see :ref:`Lake_res_model` for details).

One commonly used approach is **parametric lake models**, which relates
lake inflow and storage to outflow through functional
relationships. These formulations are computationally efficient and well suited
for large-scale applications.

Target Volume
-------------

The model also supports a simplified representation of reservoir operations
through a **target volume** mechanism.

- If the current storage is below the target volume, water is retained.
- If storage exceeds the target volume, excess water is released downstream.

The target volume is typically provided as a time series in the water management input file for each lake or reservoir. This approach enables a first-order approximation of reservoir rule curves without requiring explicit operational policies. Alternatively, it can be used by the user to incorporate more complex management information derived from expert knowledge or external water management models. For more details, see: :ref:Water management input file <WaterManagement_file>

.. note::

   This simplified operational scheme may introduce numerical uncertainties,
   particularly under rapidly varying flow conditions. Care should be taken
   when interpreting results for highly regulated systems.

Water Management Fluxes
-----------------------

In addition to physically based formulations, mizuRoute supports a
data-driven representation of lake and reservoir management.

Users can prescribe:

- Water withdrawals from lakes
- Artificial inflows (e.g., diversions)
- Managed releases

These fluxes are provided as time series, similar to the treatment of
regulated river segments. This functionality enables the representation of
human interventions in the hydrological system.

For details on the required input format, see:
:ref:`Water management input file <WaterManagement_file>`