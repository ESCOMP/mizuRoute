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


If advectiona and inertia terms are neglected,  1-D Saint-Venant equation :eq:`0.1` and :eq:`0.2` is reduced to

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

   K = \frac{A}{n}R^{frac{2}{3}}

where :eq:`5.2` is wave cerlerity [m/s], :eq:`5.3` is a wave diffusivity [m2/s], :eq:`5.4` is a channel conveyance
