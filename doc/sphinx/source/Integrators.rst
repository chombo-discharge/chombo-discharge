.. _Chap:Integrators:

Time steppers
-------------

We have various implementation of :ref:`Chap:time_stepper` that allow different temporal integration of the equations of motion. Our favorite integrator is a second order Strang splitting with implicit diffusion in the middle. However, please refer to the :doxy:`Doxygen API <classtime__stepper>` to see all available integrators. 

Typically, time steppers are selected at compile time following the style in :ref:`Chap:NewSimulations`. However, the user may select time steppers at run-time by modifying his main file in the appropriate way. For each time stepper, there are various options available at run-time through an input script.

strang2
_______

The ``strang2`` time integrator is the one we favor the most. Currently, this class performs a Strang splitting of the convection-reaction-diffusion equations in the form: 

.. math::
   :nowrap:

   \begin{equation}
   n^{k+1} = \exp\left(\frac{\Delta t}{2}A\right)\exp\left(\Delta tD\right)\exp\left(\frac{\Delta t}{2}A\right)n^{k},
   \end{equation}

where :math:`A` is the advection-reaction operator and :math:`D` is the diffusion operator. Likewise, the surface charge updates occur in the evaluation of the advective updates since charge injection is a part of the advective scheme. The RTE and electric fields are updated after each Runge-Kutta stage, although there are debugging features that allow the user to turn this off.

In ``strang2``, the diffusion operator is a second order implicit Runge-Kutta method (TGA), whereas the advection-reaction solvers are strongly stability preserving explicit Runge Kutta methods of order 2, 3, or 4. The Runge-Kutta methods support stages up to 5 and can therefore be at efficiencies close to the explicit Euler method, but with significantly higher accuracy. Note that for a five-stage second order Runge-Kutta method, the CFL limit is 4, so that the Strang splitting has a CFL limitation of 8 (since :math:`A` is advanced twice with :math:`\Delta t/2`).

The following options are available for the ``strang2`` integrator; there is also support for adaptive integration, but this is a part of the development feature and is not stable for use. The user will, for the most part, only change the Runge-Kutta method and adjust the CFL accordingly. 
   

.. literalinclude:: links/strang2.options

