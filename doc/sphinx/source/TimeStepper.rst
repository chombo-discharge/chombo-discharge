.. _Chap:time_stepper:

time_stepper
------------

The :ref:`Chap:time_stepper` class handles the integration of the plasma equations. :ref:`Chap:time_stepper` is an abstract for which we have several implements (such as a second order Runge Kutta method). Writing new temporal integrators is an extensive task, but the base class :ref:`Chap:time_stepper` contains a lot of short hand functionality for updating equations, and also contains an interface to :ref:`Chap:plasma_kinetics`. Since :ref:`Chap:time_stepper` does not actually contain an advance method for the equations of motion, we recommend that the user refers to the API of the temporal integrator that he uses (see :ref:`Chap:Solvers`) for a full explanation of the various integrators. However, the following options for :ref:`Chap:time_stepper` are passed into *all* implementation classes:

.. literalinclude:: links/time_stepper.options

The input options above are, for the most part, self-explanatory. Mostly, they refer to the handling of the size of the time step, for example by passing the Courant-Friedrichs-Lewy number, or setting a minimum or maximum possible time step. However, all of these options *may* be handled differently by different integrators, since different schemes have different restrictions on stable time steps.

Finally, there is an option to allow radiative transport updates only at certain time steps by modifying (at his own peril) the ``fast_rte`` flag. Yet again, we remark that this flag may be handled differently by different solvers. 

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
