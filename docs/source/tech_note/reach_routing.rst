River routing methods
======================

mizuRoute include five different routing methods. The routing method(s) are applied to each river reach to compute outflow from the reach. The methods are Impulse response function, Lagrangian kinmatic wave, Euler kinematic wave, Muskingum Cunge, and Diffusive wave.
This section describe each method including numerical implementation.

A fundamental equations for river reach routing start with Saint-Venant equations that consists of two equations

.. math::
   :label: 0.1

   \frac{\partial Q }{\partial x} + \frac{\partial A }{\partial t} = q_{lat}

.. math::
   :label: 0.2

   \frac{\partial Q }{\partial t} + \frac{\partial }{\partial x}(\frac{Q^{2}}{A}) + gA\frac{\partial h }{\partial x} = gA(S_{0}-S_{f})


.. _Impulse_response_function:

Impulse response function
--------------------------



.. _Lagrangian_kinematic_wave:

Lagrangian kinmatic wave
--------------------------



.. _Euler_kinematic_wave:

Euler kinmatic wave
--------------------------



.. _Muskingum-Cunge:

Muskingum-Cunge
--------------------------



.. _Diffusive_wave:

Diffusive wave
--------------------------


If advection and inertia terms are neglected,  1-D Saint-Venant equation :eq:`0.1` and :eq:`0.2` is reduced to

.. math::
   :label: 5.1

   \frac{\partial Q }{\partial t} + C\frac{\partial Q}{\partial x} = D\frac{\partial^2 Q}{\partial^2 x}

.. math::
   :label: 5.2

   C = \frac{1}{K}\frac{\partial K}{\partial A}

.. math::
   :label: 5.3

   D = \frac{K^2}{2Qw}

.. math::
   :label: 5.4

   K = \frac{A}{n}R^{\frac{2}{3}}

where :eq:`5.2` is wave cerlerity [m/s], :eq:`5.3` is a wave diffusivity [m\ :sup:`2`\/s], :eq:`5.4` is a channel conveyance.

To solve the diffusive wave equation for discharge Q, :eq:`5.1` is discretized using weighted averaged finite-difference approximations across two time steps in space
(i.e., second-order central difference in the first term of the R.H.S of :eq:`5.1` and first-order central difference for the second term of the L.H.S of :eq:`5.1`). The resulting discretized diffusive wave equation becomes

.. math::
   :label: 5.5

   ( \alpha C_{a} - 2 \beta C_{d}) \cdot Q_{j+1}^{t+1} + (2+4 \beta C_{d}) \cdot Q_{j}^{t+1} - ( \alpha C_{a} + 2 \beta C_{d}) \cdot Q_{j-1}^{t+1} = \\\\
   -[(1- \alpha )C_{d} - 2(1- \beta )C_{d})] \cdot Q_{j+1}^{t} \\\\
   + [2-4(1- \beta )C_{d}] \cdot Q_{j}^{t} \\\\
   + [(1- \alpha )C_{a} + 2(1- \beta )C_{d})] \cdot Q_{j-1}^{t} \\\\

   C_{a} = \frac{C \Delta t}{ \Delta x}, C_{d} = \frac{D \Delta t}{( \Delta x)^{2}}



where :math:`\alpha` is the weight factor for the first-order space difference approximation of the second term in :eq:`5.1`, and :math:`\beta` is a weight factor for the second-order space difference approximation of the first term :eq:`5.1`.
If both weights are set to 1, the finite difference becomes a fully implicit scheme, while setting both weights to zero results in a fully explicit scheme. If internal nodes are defined within each reach that can be specified by a user.
:eq:`5.5` can be written as a system of linear equations that can be expressed in tridiagonal matrix form and solved with the Thomas' algorithm.
For default, mizuRoute uses a fully implicit finite-difference approximation (i.e., :math:`\alpha` = :math:`\beta` = 1).
The solution of the implicit method requires downstream and upstream boundary conditions, being the latter inflow from upstream reaches. Downstream boundary condistion use the Neumann boundary condition, which specifies the gradient of discharge between the current and downstream reaches.
Note that in diffusive wave routing, celerity (C) and diffusivity (D) are updated at every time step based on the discharges (Q) and flow area (A) as opposed to IRF routing, in which celerity and diffusivity are provided as model parameters.
