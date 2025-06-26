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

blah blah

.. _Lagrangian_kinematic_wave:

Lagrangian kinmatic wave
--------------------------

blah blah

.. _Euler_kinematic_wave:

Euler kinmatic wave
--------------------------

blah blah

.. _Muskingum-Cunge:

Muskingum-Cunge
--------------------------

blah blah

.. _Diffusive_wave:

Diffusive wave
--------------------------

blah blah
