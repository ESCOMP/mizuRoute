.. _Lake_res_model:

Lake and Reservoir Models
=========================

To simulate lakes, the variable ``<is_lake_sim>`` in the control file **must be set to true**.
This flag is false by default and therefore if not specify otherwise mizuRoute does not simulate lakes.

.. list-table:: Global control keys for lake simulation
   :header-rows: 1
   :widths: 20 15 15 50
   :name: lake-global-flags

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
   :name: lake-individual-flags

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


.. _Lake_model_Endorheic:

Endorheic (closed) model
------------------------

This lake type represents **endorheic or closed lakes**, where **outflow is set to zero**. Water enters the lake only from upstream river segments or direct precipitation and is lost through **evaporation** or **water abstraction**.

This model does not require any specific parameters.

To designate a lake as endorheic, set the corresponding flag in the network topology file to ``0`` (variable identified by ``<varname_lakeModelType>``).


.. _Lake_model_Doll:

Storage-based model (Döll)
--------------------------

The least complex lake model in *mizuRoute-Lake* is the Döll formulation
(based on Döll, 2003; Hanasaki, 2006).
The Döll formulation links the outflow from the lake to the ratio of
active storage to maximum active storage through a power function.

.. list-table:: Doll lake model control keys for variables in the network topology file
   :widths: 20 20 15 15 15 15 30
   :header-rows: 1
   :name: lake-doll-parameters

   * - Control key
     - Type
     - Variable type
     - Variable dimension
     - Variable unit
     - Default
     - Description
   * - ``<varname_D03_MaxStorage>``
     - NetCDF variable name
     - real
     - seg
     - m³
     - ``-``
     - Maximum active storage for Döll formulation.
   * - ``<varname_D03_Coefficient>``
     - NetCDF variable name
     - real
     - seg
     - d⁻¹
     - ``-``
     - Release coefficient for Döll formulation.
   * - ``<varname_D03_Power>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Power scaling exponent for Döll formulation.
   * - ``<varname_D03_S0>``
     - NetCDF variable name
     - real
     - seg
     - m³
     - ``-``
     - Inactive storage for Döll formulation.

For the case when the storage is larger than inactive storage, :math:`S_0 < S`, the outflow is calculated as:

.. math::
   :name: lake-doll-equation

   O = K \, (S - S_0) \left( \frac{S - S_0}{S_{\text{max}} - S_0} \right)^P

Outflow is set to zero, :math:`O = 0`, when the storage is equal to or smaller than inactive storage, :math:`S <= S0`.

Where:

- :math:`O` = outflow from the lake segment
- :math:`S` = storage of the lake segment
- :math:`S_0` = inactive storage of the lake segment (network topology variable name is identified by control key ``<varname_D03_S0>``)
- :math:`S_{\text{max}}` = maximum active storage of the lake segment (network topology variable name is identified by control key ``<varname_D03_MaxStorage>``)
- :math:`K` = release coefficient (network topology variable name is identified by control key ``<varname_D03_Coefficient>``)
- :math:`P` = power scaling exponent (network topology variable name is identified by control key ``<varname_D03_Power>``)



.. _Lake_model_Hanasaki:

Demand-based model (Hanasaki)
-----------------------------

The Hanasaki 2006 formulation represents reservoirs with explicit consideration of water demand.
It calculates target release based on storage, inflow, and demand, differentiating between “within-a-year”
and “multi-year” reservoirs.

.. list-table:: Hanasaki lake model control keys for variables in the network topology file
   :widths: 20 20 15 15 15 15 30
   :header-rows: 1
   :name: lake-hanasaki-parameters

   * - Control key
     - Type
     - Variable type
     - Variable dimension
     - Variable unit
     - Default
     - Description
   * - ``<varname_H06_Smax>``
     - NetCDF variable name
     - real
     - seg
     - m³
     - ``-``
     - Maximum reservoir storage
   * - ``<varname_H06_alpha>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Fraction of active storage compared to total storage
   * - ``<varname_H06_envfact>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Fraction of inflow that can be used to meet demand
   * - ``<varname_H06_S_ini>``
     - NetCDF variable name
     - real
     - seg
     - m³
     - ``-``
     - Initial storage used for estimating release coefficient
   * - ``<varname_H06_c1>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Coefficient 1 for target release for irrigation reservoir
   * - ``<varname_H06_c2>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Coefficient 2 for target release for irrigation reservoir
   * - ``<varname_H06_exponent>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Exponent for actual release for “within-a-year” reservoir
   * - ``<varname_H06_denominator>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Denominator of actual release for “within-a-year” reservoir
   * - ``<varname_H06_c_compare>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Criterion to distinguish “within-a-year” vs “multi-year” reservoir
   * - ``<varname_H06_frac_Sdead>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Fraction of dead storage to maximum storage
   * - ``<varname_H06_E_rel_ini>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Initial release coefficient
   * - ``<varname_H06_I_Jan>`` … ``<varname_H06_I_Dec>``
     - NetCDF variable name
     - real
     - seg
     - m³ s⁻¹
     - ``-``
     - Average monthly inflow for each month
   * - ``<varname_H06_D_Jan>`` … ``<varname_H06_D_Dec>``
     - NetCDF variable name
     - real
     - seg
     - m³ s⁻¹
     - ``-``
     - Average monthly demand for each month
   * - ``<varname_H06_purpose>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Reservoir purpose flag (0 = non-irrigation, 1 = irrigation)
   * - ``<varname_H06_I_mem_F>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Flag to transition to modelled inflow
   * - ``<varname_H06_D_mem_F>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Flag to transition to modelled/provided demand
   * - ``<varname_H06_I_mem_L>``
     - NetCDF variable name
     - int
     - seg
     - year
     - ``-``
     - Memory length in years for inflow
   * - ``<varname_H06_D_mem_L>``
     - NetCDF variable name
     - int
     - seg
     - year
     - ``-``
     - Memory length in years for demand



.. _Lake_model_HYPE:

Elevation-based model (Hydropower Reservoir Formulation from HYPE)
------------------------------------------------------------------

The HYPE formulation describes the representation of a hydropower reservoir in *mizuRoute-Lake*.
This includes parameters for spillways, turbine operations, and reservoir management rules.

.. list-table:: HYPE lake model control keys for variables in the network topology file
   :widths: 20 20 15 15 15 15 30
   :header-rows: 1
   :name: lake-hype-parameters

   * - Control key
     - Type
     - Variable type
     - Variable dimension
     - Variable unit
     - Default
     - Description
   * - ``<varname_HYP_E_emr>``
     - NetCDF variable name
     - real
     - seg
     - m
     - ``-``
     - Elevation of emergency spillway
   * - ``<varname_HYP_E_lim>``
     - NetCDF variable name
     - real
     - seg
     - m
     - ``-``
     - Elevation below which primary spillway flow is restricted
   * - ``<varname_HYP_E_min>``
     - NetCDF variable name
     - real
     - seg
     - m
     - ``-``
     - Elevation below which outflow is zero
   * - ``<varname_HYP_E_zero>``
     - NetCDF variable name
     - real
     - seg
     - m
     - ``-``
     - Elevation at which lake/reservoir storage is zero
   * - ``<varname_HYP_Qrate_emr>``
     - NetCDF variable name
     - real
     - seg
     - m³ s⁻¹
     - ``-``
     - Emergency rate of flow for each unit of elevation above HYP_E_emr
   * - ``<varname_HYP_Erate_emr>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Power for the emergency spillway exponential flow curve
   * - ``<varname_HYP_Qrate_prim>``
     - NetCDF variable name
     - real
     - seg
     - m³ s⁻¹
     - ``-``
     - Average yearly/long-term output from primary spillway
   * - ``<varname_HYP_Qrate_amp>``
     - NetCDF variable name
     - real
     - seg
     - –
     - ``-``
     - Amplitude of the primary spillway outflow
   * - ``<varname_HYP_Qrate_phs>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Phase of the primary spillway outflow (day of year; default = 100)
   * - ``<varname_HYP_prim_F>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Reservoir primary spillway flag (1 if present, else 0)
   * - ``<varname_HYP_A_avg>``
     - NetCDF variable name
     - real
     - seg
     - m²
     - ``-``
     - Average area of lake (unused if bathymetry is provided)
   * - ``<varname_HYP_Qsim_mode>``
     - NetCDF variable name
     - int
     - seg
     - –
     - ``-``
     - Outflow calculation mode (1 = sum of emergency + primary spillway; else = maximum of emergency or primary spillway)


For hydropower reservoirs, a sinusoidal function defines the target hydropower production outflow.
This function is shifted in time based on a day of the year, :math:`B_{\mathrm{phase}}`, as:

.. math::
   :label: HYPE_sin_equation

   F_{\mathrm{sin}} = \max \Big(0, 1 + A_{\mathrm{amp}} \sin\Big(\frac{2 \pi D_{\mathrm{julian}} + B_{\mathrm{phase}}}{365}\Big) \Big)

Next, the limiting factor is defined when the lake elevation is between :math:`E_{\mathrm{prim}}` and :math:`E_{\mathrm{lim}}`.
The linear scaling for restricted hydropower production is:

.. math::
   :label: HYPE_lim_equation

   F_{\mathrm{lim}} = \min \Big( \max \Big( \frac{E - E_{\mathrm{prim}}}{E_{\mathrm{lim}} - E_{\mathrm{prim}}}, 0 \Big), 1 \Big)

If the water level is below :math:`E_{\mathrm{prim}}`, :math:`F_{\mathrm{lim}} = 0`.
If the water level is above :math:`E_{\mathrm{lim}}`, :math:`F_{\mathrm{lim}} = 1`.

The production outflow for hydropower is then calculated as:

.. math::
   :label: HYPE_main_equation

   Q_{\mathrm{main}} = F_{\mathrm{sin}} \, F_{\mathrm{lim}} \, F_{\mathrm{managed}} \, Q_{\mathrm{avg,rate}}

If the reservoir elevation, :math:`E`, exceeds the emergency spillway elevation, :math:`E_{\mathrm{emg}}`, the emergency spillway is activated:

.. math::
   :label: HYPE_emg_equation

   Q_{\mathrm{emg}} = Q_{\mathrm{emg,rate}} (E - E_{\mathrm{emg}})^{P_{\mathrm{emg}}}

Finally, the outflow from the reservoir is either the maximum of :math:`Q_{\mathrm{emg}}` and :math:`Q_{\mathrm{main}}` or their summation (depending on mizuRoute settings):

.. math::
   :label: HYPE_outflow_equation

   O = \max(Q_{\mathrm{emg}}, Q_{\mathrm{main}})

Where the parameters are defined as:

- :math:`A_{\mathrm{amp}}` = amplitude of the sinusoidal function (network topology variable name is identified by control key ``<varname_HYP_Qrate_amp>``)
- :math:`B_{\mathrm{phase}}` = phase shift for the sinusoidal function (network topology variable name is identified by control key ``<varname_HYP_Qrate_phs>``)
- :math:`E_{\mathrm{prim}}` = primary spillway elevation (flow restricted below this) (network topology variable name is identified by control key ``<varname_HYP_E_prim>``)
- :math:`E_{\mathrm{lim}}` = elevation at which primary spillway flow is unrestricted (network topology variable name is identified by control key ``<varname_HYP_E_lim>``)
- :math:`F_{\mathrm{managed}}` = management factor (optional control) (network topology variable name is identified by control key ``<varname_HYP_prim_F>``)
- :math:`Q_{\mathrm{avg,rate}}` = average rated outflow of primary spillway/turbine (network topology variable name is identified by control key ``<varname_HYP_Qrate_prim>``)
- :math:`Q_{\mathrm{emg,rate}}` = emergency spillway flow coefficient (network topology variable name is identified by control key ``<varname_HYP_Qrate_emr>``)
- :math:`P_{\mathrm{emg}}` = emergency spillway exponent (network topology variable name is identified by control key ``<varname_HYP_Erate_emr>``)
- :math:`D_{\mathrm{julian}}` = Julian day of the year
- :math:`E` = reservoir elevation
- :math:`F_{\mathrm{sin}}` = sinusoidal target flow fraction
- :math:`F_{\mathrm{lim}}` = limiting factor due to reservoir elevation
- :math:`Q_{\mathrm{emg}}` = emergency spillway outflow
- :math:`Q_{\mathrm{main}}` = main hydropower production outflow
- :math:`O` = final outflow from the reservoir (m³/s) (network topology variable name is identified by control key ``<varname_HYP_Qsim_mode>``)