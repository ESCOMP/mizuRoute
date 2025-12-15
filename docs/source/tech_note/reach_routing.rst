
.. _Reach_routing_oveall_workflow:

Overall workflow
======================

Overall computation workflow is shown in :numref:`Figure_overall_comp_workflow`.

#. Compute HRU (Hydrologic Response Unit, or simply catchment) mean runoff [m/s], :math:`R_{lat}`, if runoff is given at hydrologic model HRU.

#. Convert :math:`R_{lat}` from depth unit to volume (:math:`R_{lat}` times HRU area) to get lateral runoff volume (:math:`q_{lat}`) [m3/s]

#. Perform overland routing to delay lateral runoff volume.

#. Route delayed lateral discharge at each river reach through river network.

The overland routing is optional, and currently simple unit hydrograph based on gamma distribution is used to delay instantaneous runoff.

.. _Figure_overall_comp_workflow:

.. figure:: images/overall_comp_workflow.png
 :width: 700
 :height: 400

 Overall routing procedures from runoff input into model to streamflow computation.

.. _Hillslope_routing_scheme:

Hillslope routing scheme
========================

MizuRoute uses a similar concept to a unit-hydrograph commonly used in engineering hydrology to account for travel time of instantaneous runoff to the river reach (delay and attenuate runoff).
A unit-hydrograph is defined by a hydrograph (time series of discharge) that is derived from a unit depth of excess rainfall on a drainage area within a specific time period.
Here, this concept is applied to directly runoff volume over a HRU, instead of rainfall excess over a drainage area.
Therefore, a unit-hydrograph represents a time series of lateral discharge into a river reach from a HRU derived from a unit volume of runoff.
Here, a probability density function (PDF) is used as a unit-hydrograph, so that cumulative sum of PDF is 1.0.
This means that runoff volume at a current time step is just distributed in the future, with the sum of future distributed runoff is equal to the current runoff volume (i.e., volume conserved)

To get actual delayed lateral flow series to the river reach, unit-hydrograph convolution is performed as below:

.. _Figure_uh_convolution:

.. figure:: images/uh_convolution.png
 :width: 700
 :height: 500

 Illustration of discrete Unit hydrograph convolution.

In mizuRoute, gamma distribution is used for PDF-based unit-hydrograph and written as:

.. math::
   :label: gamma_distribution

   f(t; a, \theta) = \frac{1}{\Gamma(a)\theta^{a}}t^{a - 1} e^{-\frac{t}{\theta}},
   \quad t > 0

where *t* is time [sec], *a* is a shape parameter [–] (a > 0), and :math:`\theta` is a timescale parameter [sec].
Both the shape and timescale parameters affect the peak time of (mode of the distribution: :math:`(a - 1)\theta` and flashiness (variance of the distribution: :math:`a\theta^2` of the unit-hydrograph (UH).
These UH should depend on the physical HRU characteristics. Continuous gamm distribution is descritized per time step before the convolution performed.
These shape and scale parameters are currently specified as a spatially constant parameter (see :ref:`Spatially-constant parameter namelist <namelist_file>`), though they could be provided as spatially distributed parameters from river data netCDF (potential future implementation)

Also, please see section 3.1 in :ref:`Mizukami et al. 2016 <Mizukami2016>` for hillslope routing theory.

.. _River_routing_schemes:

River routing schemes
======================

For river reach routing, mizuRoute include five different routing methods. The routing method(s) are applied to each river reach to compute outflow from the reach.
The methods are Impulse response function (routOpt=1 in mizuRoute), Lagrangian kinmatic wave (routOpt=2), Euler kinematic wave (routOpt=3), Muskingum Cunge (routOpt=4)), and Diffusive wave (routOpt=5).

This section describes each scheme including numerical implementation.
Impulse response function and lagrangian kinematic wave, implemented eariler, are also described in :ref:`Mizukami et al. (2016) <Mizukami2016>`.
The othter schemes were implemented afterwards, and described in appendices of :ref:`Cortes-Salazar et al. (2023) <Cortes2023>` (except for Euler kinematic wave).

A fundamental equations for river reach routing start with Saint-Venant equations that consists of two equations

.. math::
   :label: 0.1

   \frac{\partial Q }{\partial x} + \frac{\partial A }{\partial t} = q_{lat}

.. math::
   :label: 0.2

   \frac{\partial Q }{\partial t} + \frac{\partial }{\partial x}(\frac{Q^{2}}{A}) + gA\frac{\partial h }{\partial x} = gA(S_{0}-S_{f})

where *Q* is a discharge [m\ :sup:`3`\/s] at time t and a point of reach x,
*A* is a flow cross-sectional area [m\ :sup:`2`],
*h* is flow height [m],
:math:`S_{0}` is a slope of reach [m/m],
:math:`S_{f}` is a friction slope [m/m].
LHS of :eq:`0.2` consists of advection, inertia, and pressure gradient from the 1st to 3rd terms, while force temrs of RHS of :eq:`0.2` consists of gravity and frictional force from a river bed.

The frictional slope is written as:

.. math::
   :label: 0.3

   S_{f} = \frac{Q^2}{K^2}

.. math::
   :label: 0.4

   K = \frac{A}{n}R^{\frac{2}{3}}

where *K* (Eq. :eq:`0.4`) is a channel conveyance. *n* is manning coefficient [-] and *R* is hydraulic radius [m].

If advection and inertia terms are neglected (i.e., the 1st and 2nd terms in LHS of :eq:`0.2`) and without lateral flow in :eq:`0.1`, 1-D Saint-Venant equation :eq:`0.1` and :eq:`0.2` is reduced to

.. math::
   :label: 0.5

   \frac{\partial Q }{\partial t} + C\frac{\partial Q}{\partial x} = D\frac{\partial^2 Q}{\partial^2 x}

.. math::
   :label: 0.6

   C = \frac{1}{K}\frac{\partial K}{\partial A}

.. math::
   :label: 0.7

   D = \frac{K^2}{2QB}

where *C* (Eq. :eq:`0.6`) is a wave celerity [m/s] and *D* (Eq. :eq:`0.7`) is a diffusivity [m\ :sup:`2`\/s]. *B* is a top width of flow cross-sectional area [m].

If *D* is set to zero (i.e., diffusion is neglected), Eq. :eq:`0.5` becomes kinematic wave equation.
The other way to derive kinematic wave equation is to neglect pressure gradient term in addition to advection and inertia and pressure gradient terms (i.e., all the terms in LHS of Eq. :eq:`0.2`).

Eq. :eq:`0.5` serve a starting point of numerical implementation of
Impulse response function (section :numref:`Impulse_response_function`),
Euler Kinematic wave (section :numref:`Euler_kinematic_wave`),
Muskingum-Cunge (section :numref:`Muskingum-Cunge`),
Diffusive wave schemes (section :numref:`Diffusive_wave`).

The numerical implementations of Euler kinematic wave and Diffusive wave are essentially identical. Therefore, a user is referred to section :numref:`Diffusive_wave` for Euler kinematic wave numerical schemes.


.. _Impulse_response_function:

Impulse response function
--------------------------

The impulse response function (IRF) is derived from 1-D diffusive wave equation (Eq. :eq:`0.5`). Eq. :eq:`0.5` can be solved using convolution integrals:

.. math::
   :label: 0.8

   Q(x,t) = \int_{0}^t Q_{in}(t-s)h(x,s)ds

where :math:`Q_{in}` is a pulse or here input runoff volume in the river reach (e.g., inflow from HRU, or/and upstream reach), and :math:`h(x,s)` is the green function or *impulse response function*, and *x* is the distance from the inflow point in the reach to the reach location where discharge desires to be computed (i.e., HRU outlet).

.. math::
   :label: 0.9

   h(x,t) = \frac{x}{2t \sqrt{\pi t D}} \exp(-\frac{(Ct-x)^2}{4Dt})

where *C* (Eq. :eq:`0.6`) is a wave celerity [m/s] and *D* (Eq. :eq:`0.7`) is a diffusivity [m\ :sup:`2`\/s].

This is the same concept as river routing described in Lohmann et al., (1996). Also, please see section 3.2.2 in :ref:`Mizukami et al. 2016 <Mizukami2016>`.

Lohmann et al., (1996) and original mizuRoute (version 1) use the *source-to-sink* approach where unique IRFs are developed for each HRU outlet (sink) to its all upstream HRUs (sources), the IRFs are applied with lateral flows from the corresponding HRUs,
then routed flow from all the source HRUs are summed up to the river discharge at the target HRU outlet. This process is repeated for all the HRUs and each time step.

One modification made from Lohmann et al., 1996 and mizuRoute version 1 was to develop reach-specific IRF based on reach length, *C* and *D*, then apply it with inflow from the immediate upstream reach(es) to move water from upstream end to the outlet.
The computed river discharge at each river reach is used as inflow of the next downstream reach. Therefore, the IRF computation is performed from upstream HRUs to downstream in the order.
This approach is consistent with the other routing methods used in mizuRoute, and also requires less memory usage of IRF ordinates than source-to-sink approach, which requires storing IRFs of all the upstream HRUs for all the HRUs.

.. _Lagrangian_kinematic_wave:

Lagrangian kinmatic wave
--------------------------

Also, please see section 3.2.1 in Mizukami et al. 2016

.. _Euler_kinematic_wave:

Euler kinmatic wave
--------------------------

See section :numref:`Diffusive_wave` for details on numerical implementatin of Euler kinematic wave.


.. _Muskingum-Cunge:

Muskingum-Cunge
--------------------------


Muskingum-Cunge (M-C) routing formulation begins with a kinematic wave equation, :eq:`0.5` with *D* set to zero.
The kinematic wave equation can be discretized with weight factors X and Y to give:

.. math::
   :label: 4.1

   \frac{X(I_{t+1}-I_{t})+(1-X)(O_{t+1}-O_{t})}{\Delta t} + C \frac{Y(O_{t}-I_{t})+(1-Y)(O_{t+1}-I_{t})}{\Delta x}=0

where :math:`I_{t+1}` and :math:`I_{t}` are inflow to a reach segment (length is :math:`\Delta x`) at the end and beginning of the time step (time step is :math:`\Delta t` ) and :math:`O_{t+1}` and :math:`O_{t}` are outflow from a reach segment at the end and beginning of the time step.
The spatial weight factor Y is set to 0.5 and then Eq. :eq:`4.1` is rearranged, giving:

.. math::
   :label: 4.2

   O_{t+1} = \frac{-X+0.5 C_{n}}{1-X+0.5 C_{n}} I_{t+1} + \frac{X+0.5 C_{n}}{1-X+0.5 C_{n}} I_{t} + \frac{1-X-0.5 C_{n}}{1-X+0.5 C_{n}} O_{t}

where :math:`C_{n}` is Courant Number defined by :math:`C \frac{\Delta t}{\Delta x}`. Eq :eq:`4.2` is generally called Muskingum equation,
but Cunge (1969) found that the numerical diffusion in the explicit solution of Eq :eq:`4.2`, which can happen depending on weight factors, can match the physical diffusion by setting X (along with Y=0.5) to:

.. math::
   :label: 4.3

   X=0.5(1-\frac{Q}{BS_{0} C\Delta x})

where :math:`S_{0}` is the reach slope, *B* is a top widith of flow cross-section area. Here discharge *Q* and *B* can be estimated by 3-point Q values (:math:`I_{t+1}`, :math:`I_{t}`, and :math:`O_{t}`).
Note that *B* is a function of Q given channel cross-section assumption (see section x-x).
At every time step and reach, temporal weight factor X is update based on given 3-point discharge values. Since Muskingum-Cunge is explicitly solved, the solution can be unstable.
To stabilize the solution, the sub time step (:math:`\Delta t`) is determined at every simulation step so that the Courant number is less than unity


.. _Diffusive_wave:

Diffusive wave
--------------------------


To solve the diffusive wave equation for discharge Q, Eq. :eq:`0.5` is discretized using weighted averaged finite-difference approximations across two time steps in space
(Figure 1; i.e., second-order central difference in the RHS of :eq:`0.5` and first-order central difference for the second term of the LHS of :eq:`0.5`).

.. _Figure diffusive wave numerical discretization:

.. figure:: images/dw_discretization.png
 :width: 600

 Space and time discretization used for numerical solution of diffusive wave equation

The resulting discretized diffusive wave equation becomes:

.. math::
   :label: 5.5

   \frac{Q_{j}^{t+1} - Q_{j}^{t}}{\Delta t} + \frac{C}{2 \Delta x} \cdot ((1- \alpha )(Q_{j+1}^{t} - Q_{j-1}^{t})+ \alpha (Q_{j+1}^{t+1} - Q_{j-1}^{t+1})) = \\\\
   D \cdot (\frac{(1- \beta)(Q_{j+1}^{t} - 2Q_{j}^{t} + Q_{j-1}^{t})}{(\Delta x)^2} + \frac{\beta (Q_{j+1}^{t+1} - 2Q_{j}^{t+1} +Q_{j-1}^{t+1})}{(\Delta x)^2})

Rearranging Eq. :eq:`5.5` to:

.. math::
   :label: 5.6

   ( \alpha C_{a} - 2 \beta C_{d}) \cdot Q_{j+1}^{t+1} + (2+4 \beta C_{d}) \cdot Q_{j}^{t+1} - ( \alpha C_{a} + 2 \beta C_{d}) \cdot Q_{j-1}^{t+1} = \\\\
   -[(1- \alpha )C_{d} - 2(1- \beta )C_{d})] \cdot Q_{j+1}^{t} \\\\
   + [2-4(1- \beta )C_{d}] \cdot Q_{j}^{t} \\\\
   + [(1- \alpha )C_{a} + 2(1- \beta )C_{d})] \cdot Q_{j-1}^{t} \\\\

   C_{a} = \frac{C \Delta t}{ \Delta x}, C_{d} = \frac{D \Delta t}{( \Delta x)^{2}}

where :math:`\alpha` is the weight factor for the first-order space difference approximation of the second term of the LHS of :eq:`0.5`, and :math:`\beta` is a weight factor for the second-order space difference approximation in RHS of :eq:`0.5`.
If both weights are set to 1, the finite difference becomes a fully implicit scheme, while setting both weights to zero results in a fully explicit scheme. For default, mizuRoute uses a fully implicit finite-difference approximation (i.e., :math:`\alpha` = :math:`\beta` = 1).
Note that celerity (C) and diffusivity (D) include Q, which means the diffusive equation is actually non-linear. Here celerity (*C*) and diffusivity (*D*) are updated at every time step based on the discharges (Q) and flow area (A) at previous time step to liearize the diffusive equation.
Note that IRF routing is also based on diffusve equation. a major difference is that in IRF routing, celerity and diffusivity are provided as model parameters and constant in time, though they can be spatially distributed.

To apply the numerical solution of discretized diffusive wave equation for each reach, the internal nodes need to be defined within each reach.
The number of internal node is now hard-coded as 5 (in future, this will be made available as a control variable so that the number of the internal nodes can bespecified by a user via a control file.

:eq:`5.6` can be written as a system of linear equations that can be expressed in tridiagonal matrix form, :math:`A \cdot Q=b`, which can be solved with  with the Thomas' algorithm.

.. _Figure 4 internal nodes in a reach:

.. figure:: images/4_internal_nodes.png
 :width: 600

 An example of 4 internal nodes per reach.

For example, with 4 internal nodes as shown in, the matrix form of the equations are written as:

.. math::
   :label: 5.8

   \small A=
   \left[ \begin {array}{cccc}
   1&0&0&0&0\cr
   -(\alpha C_{d}+2\beta C_{d})&2+4\beta C_{d}&\alpha C_{a}-2\beta C_{d}&0&0\cr
   0&-(\alpha C_{d}+2\beta C_{d})&2+4\beta C_{d}&\alpha C_{a}-2\beta C_{d}&0\cr
   0&0&-(\alpha C_{d}+2\beta C_{d})&2+4\beta C_{d}&\alpha C_{a}-2\beta C_{d}\cr
   0&0&0&-1&1
   \end {array} \right]

.. math::
   :label: 5.7

   \small Q=
   \left[ \begin {array}{c}
   Q_{1}^{t+1} \cr
   Q_{2}^{t+1} \cr
   Q_{3}^{t+1} \cr
   Q_{4}^{t+1} \cr
   Q_{5}^{t+1}
   \end {array} \right]

.. math::
   :label: 5.9

   \small b=
   \left[ \begin {array}{c}
   Q_{1}^{t+1} \cr
   ((1-\alpha)C_{a} + 2(1-\beta)C_{d}) \cdot Q_{1}^{t} + (2-4(1-\beta)C_{d}) \cdot Q_{2}^{t} - ((1-\alpha)C_{a}-2(1-\beta)C_{d}) \cdot Q_{3}^{t} \cr
   ((1-\alpha)C_{a} + 2(1-\beta)C_{d}) \cdot Q_{2}^{t} + (2-4(1-\beta)C_{d}) \cdot Q_{3}^{t} - ((1-\alpha)C_{a}-2(1-\beta)C_{d}) \cdot Q_{4}^{t} \cr
   ((1-\alpha)C_{a} + 2(1-\beta)C_{d}) \cdot Q_{3}^{t} + (2-4(1-\beta)C_{d}) \cdot Q_{4}^{t} - ((1-\alpha)C_{a}-2(1-\beta)C_{d}) \cdot Q_{5}^{t} \cr
   a \cdot dx
   \end {array} \right]


The top row of the system of equations is upstream boundary conditions, which is inflow from upstream reaches (i.e., Dirichlet boundary condition).
The Bottom row of the system of equations is downstream boundary condition.
Here, Neumann boundary condition, which specifies the gradient of discharge between two adjacent nodes at the downstream end, is used.
Neumann boundary condition at the downstream end is written by:

.. math::
   :label: 5.10

   \frac{\partial Q}{\partial x}\Big{|}_{x=5}

which is discretized as :math:`Q_{5}^{t+1} - Q_{4}^{t+1} = a \cdot dx`. The gradient at downstream end :math:`a` is approximated by the Q computed at the nodes at previous time step.

What makes this numerical solution become **kinematic wave solution** is simply to set *D* to zero.

