.. _Lake_model:

Lake Model
==========

To simulate lakes, the variable ``<is_lake_sim>`` in the control file **must always be set to true**. 
If it is false, lakes will not be simulated.

.. list-table:: General Control Variable in Control File
   :widths: 25 75
   :header-rows: 1

   * - Variable
     - Description
   * - <is_lake_sim>
     - Logical; indicates whether lakes are simulated
   * - <lakeRegulate>
     - Logical; F -> all lakes are treated as natural (lakeType=1) regardless of individual lakeModelType, T (default)
   * - <LakeInputOption>
     - Fluxes for lake simulation; 0 -> evaporation + precipitation (default), 1 -> runoff, 2 -> evaporation + precipitation + runoff



The following variables need to be specified in the network topology data for each segment that may contain a lake.

.. list-table:: Lake Variables that need to be specified in Network Topology
   :widths: 25 75
   :header-rows: 1

   * - Variable
     - Description
   * - <varname_islake>
     - Flag to define a lake (1 = lake, 0 = reach)
   * - <varname_lakeModelType>
     - Defines the lake model type (1 = Döll, 2 = Hanasaki, 3 = HYPE, 0 = Endorehic)
   * - <varname_LakeTargVol>
     - Flag to follow the provided target volume (1 = yes, 0 = no)

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

