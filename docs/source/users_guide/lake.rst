.. _Lake_model:

Lake Model
==========

To simulate lakes, the variable ``<is_lake_sim>`` in the control file **must be set to true**.
This flag is false by default and therefore if not specify otherwise mizuRoute does not simulate lakes.

.. list-table:: Global control keys for lake simulation
   :header-rows: 1
   :widths: 20 15 15 50

   * - Control key
     - Type
     - Default
     - Description
   * - ``<is_lake_sim>``
     - Logical (Global Flag)
     - ``F``
     - Indicates whether lakes are simulated.

       * ``F`` → lakes are not simulated.
       * ``T`` → lakes are included in the routing.
   * - ``<lakeRegulate>``
     - Logical (Global Flag)
     - ``T``
     - Controls whether all lakes are treated as natural or regulated.

       * ``F`` → all lakes are treated as natural (``lakeType = 1``) regardless of individual ``lakeModelType``.
       * ``T`` → regulation is applied (default).
   * - ``<LakeInputOption>``
     - Integer (Global Flag)
     - ``0``
     - Selects the type of fluxes used in lake simulation:

       * ``0`` → evaporation + precipitation (default, ignores runoff provided for lakes and reservoirs)
       * ``1`` → runoff (for cases that simulated runoff account for processes such as snow melt, ice formation, etc – resolving lakes and reservoirs water and energy balances in land surface model for example)
       * ``2`` → evaporation + precipitation + runoff




The following variables need to be specified in the network topology data for each element of network topology that is identified as lake.

.. list-table:: Lake-related control keys in the network topology file
   :widths: 20 20 15 15 15 15 30
   :header-rows: 1

   * - Control key
     - Type
     - Variable type
     - Variable dimension
     - Variable unit
     - Default
     - Description
   * - ``<varname_islake>``
     - NetCDF variable name
     - int
     - seg
     - flag (0/1)
     - ``-``
     - Flag to define whether the segment is a lake (``1`` = lake, ``0`` = reach).
   * - ``<varname_lakeModelType>``
     - NetCDF variable name
     - int
     - seg
     - categorical
     - ``-``
     - Defines the lake model type for the segment:
       ``0`` = Endorheic, ``1`` = Döll, ``2`` = Hanasaki, ``3`` = HYPE.
   * - ``<varname_LakeTargVol>``
     - NetCDF variable name
     - int
     - seg
     - flag (0/1)
     - ``-``
     - Flag to follow the provided target volume for the lake (``1`` = yes, ``0`` = no).



For further reading about the below formulation, please see
`Gharari et al., 2024 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR032400>`_.


.. _Lake_model_Doll:

Storage-based model (Döll)
--------------------------

The least complex lake model in *mizuRoute-Lake* is the Döll formulation
(based on Döll, 2003; Hanasaki, 2006).
The Döll formulation links the outflow from the lake to the ratio of
active storage to maximum active storage through a power function.

.. list-table::
   :widths: 20 15 10 10 45
   :header-rows: 1

   * - Variable
     - Dimension
     - Unit
     - Type
     - Description
   * - <varname_D03_MaxStorage>
     - seg
     - m³
     - real
     - Maximum active storage for Doll 2003 formulation
   * - <varname_D03_Coefficient>
     - seg
     - d⁻¹
     - real
     - Coefficient for Doll 2003 formulation (release coefficient)
   * - <varname_D03_Power>
     - seg
     - –
     - real
     - Power scaling for Doll 2003 formulation
   * - <varname_D03_S0>
     - seg
     - m³
     - real
     - Inactive storage for Doll 2003 formulation


.. _Lake_model_HYPE:

Elevation-based model (HYPE - Hydropower Reservoir Formulation)
---------------------------------------------------------------

The HYPE formulation describes the representation of a hydropower reservoir in *mizuRoute-Lake*.
This includes parameters for spillways, turbine operations, and reservoir management rules.

.. list-table::
   :widths: 20 15 10 10 45
   :header-rows: 1

   * - Variable
     - Dimension
     - Unit
     - Type
     - Description
   * - <varname_HYP_E_emr>
     - seg
     - m
     - real
     - Elevation of emergency spillway
   * - <varname_HYP_E_lim>
     - seg
     - m
     - real
     - Elevation below which primary spillway flow is restricted
   * - <varname_HYP_E_min>
     - seg
     - m
     - real
     - Elevation below which outflow is zero
   * - <varname_HYP_E_zero>
     - seg
     - m
     - real
     - Elevation at which lake/reservoir storage is zero
   * - <varname_HYP_Qrate_emr>
     - seg
     - m³ s⁻¹
     - real
     - Emergency rate of flow for each unit of elevation above HYP_E_emr
   * - <varname_HYP_Erate_emr>
     - seg
     - –
     - real
     - Power for the emergency spillway exponential flow curve
   * - <varname_HYP_Qrate_prim>
     - seg
     - m³ s⁻¹
     - real
     - Average yearly/long-term output from primary spillway
   * - <varname_HYP_Qrate_amp>
     - seg
     - –
     - real
     - Amplitude of the primary spillway outflow
   * - <varname_HYP_Qrate_phs>
     - seg
     - –
     - real
     - Phase of the primary spillway outflow (day of year; default = 100)
   * - <varname_HYP_prim_F>
     - seg
     - –
     - real
     - Reservoir primary spillway flag (1 if present, else 0)
   * - <varname_HYP_A_avg>
     - seg
     - m²
     - real
     - Average area of lake (unused if bathymetry is provided)
   * - <varname_HYP_Qsim_mode>
     - seg
     - –
     - real
     - Outflow calculation mode (1 = sum of emergency + primary spillway; else = maximum)


.. _Lake_model_Hanasaki:

Demand-based model (Hanasaki)
-----------------------------

The Hanasaki 2006 formulation represents reservoirs with explicit consideration of water demand.
It calculates target release based on storage, inflow, and demand, differentiating between “within-a-year”
and “multi-year” reservoirs.

.. list-table::
   :widths: 20 15 10 10 45
   :header-rows: 1

   * - Variable
     - Dimension
     - Unit
     - Type
     - Description
   * - <varname_H06_Smax>
     - seg
     - m³
     - real
     - Maximum reservoir storage
   * - <varname_H06_alpha>
     - seg
     - –
     - real
     - Fraction of active storage compared to total storage
   * - <varname_H06_envfact>
     - seg
     - –
     - real
     - Fraction of inflow that can be used to meet demand
   * - <varname_H06_S_ini>
     - seg
     - m³
     - real
     - Initial storage used for estimating release coefficient
   * - <varname_H06_c1>
     - seg
     - –
     - real
     - Coefficient 1 for target release for irrigation reservoir
   * - <varname_H06_c2>
     - seg
     - –
     - real
     - Coefficient 2 for target release for irrigation reservoir
   * - <varname_H06_exponent>
     - seg
     - –
     - real
     - Exponent for actual release for “within-a-year” reservoir
   * - <varname_H06_denominator>
     - seg
     - –
     - real
     - Denominator of actual release for “within-a-year” reservoir
   * - <varname_H06_c_compare>
     - seg
     - –
     - real
     - Criterion to distinguish “within-a-year” vs “multi-year” reservoir
   * - <varname_H06_frac_Sdead>
     - seg
     - –
     - real
     - Fraction of dead storage to maximum storage
   * - <varname_H06_E_rel_ini>
     - seg
     - –
     - real
     - Initial release coefficient
   * - <varname_H06_I_Jan> … <varname_H06_I_Dec>
     - seg
     - m³ s⁻¹
     - real
     - Average monthly inflow for each month
   * - <varname_H06_D_Jan> … <varname_H06_D_Dec>
     - seg
     - m³ s⁻¹
     - real
     - Average monthly demand for each month
   * - <varname_H06_purpose>
     - seg
     - –
     - real
     - Reservoir purpose flag (0=non-irrigation, 1=irrigation)
   * - <varname_H06_I_mem_F>
     - seg
     - –
     - real
     - Flag to transition to modelled inflow
   * - <varname_H06_D_mem_F>
     - seg
     - –
     - real
     - Flag to transition to modelled/provided demand
   * - <varname_H06_I_mem_L>
     - seg
     - year
     - real
     - Memory length in years for inflow
   * - <varname_H06_D_mem_L>
     - seg
     - year
     - real
     - Memory length in years for demand

