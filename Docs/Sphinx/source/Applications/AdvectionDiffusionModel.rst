.. _Chap:AdvectionDiffusionModel:

Advection-diffusion model
=========================

The advection-diffusion model simply advects a scalar quantity with EB and AMR.
The equation of motion is

.. math::

   \frac{\partial\phi}{\partial t} + \nabla\cdot\left(\mathbf{v}\phi - D\nabla\phi\right) = 0.

The full implementation for this model consists of the following classes:

* ``AdvectionDiffusionStepper``, which implements :ref:`Chap:TimeStepper`.
* ``AdvectionDiffusionSpecies``, which implements :ref:`Chap:CdrSpecies`, and thus sparses the initial condition into the problem.
* ``AdvectionDiffusionTagger``, which implements :ref:`Chap:CellTagger` and flags cells for refinement and coarsening.

This module only uses :ref:`Chap:CdrSolver`.

.. tip::
   
   The source code for this problem is found in :file:`$DISCHARGE_HOME/Physics/AdvectionDiffusion`.
   See `AdvectionDiffusionStepper <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhysics_1_1AdvectionDiffusion_1_1AdvectionDiffusionStepper.html>`_ for the C++ API for this time stepper.

Time stepping
-------------

Two time stepping algorithms are supported:

#. A second-order Runge-Kutta method (Heun's method).
#. An implicit-explicit method (IMEX) which uses corner-transport upwind (CTU, see :ref:`Chap:CdrCTU`) and a Crank-Nicholson diffusion solve.

   .. important::

      The Crank-Nicholson method is known to be marginally stable.

See :ref:`Chap:AdvectionDiffusionOptions` to see how to select between these two algorithms.

Heun's method
_____________

For Heun's method we perform the following steps:
Let :math:`\left[\nabla\cdot\mathbf{J}\left(\phi\right)\right]` be the finite-volume approximation to :math:`\nabla\cdot\left(\mathbf{v}\phi - D\nabla\phi\right)`, where :math:`\mathbf{J}\left(\phi\right) = \mathbf{v}\phi - D\nabla\phi`.
The Runge-Kutta advance is then

.. math::

   \phi^\prime &= \phi^k - \Delta t\left[\nabla\cdot\mathbf{J}\left(\phi^k\right)\right]\\
   \phi^{k+1} &= \phi^k - \frac{\Delta t}{2}\left(\left[\nabla\cdot\mathbf{J}\left(\phi^k\right)\right] + \left[\nabla\cdot\mathbf{J}\left(\phi^\prime\right)\right]\right)

.. warning::
   
   Note that when diffusion and advecting is coupled in this way, we do not include the transverse terms in the :ref:`Chap:CdrCTU` discretization and limit the time step by

   .. math::

      \Delta t \leq \frac{1}{\frac{\sum_{i=1}^d |v_i|}{\Delta x} + \frac{2dD}{\Delta x^2}},

   where :math:`d` is the dimensionality and :math:`D` is the diffusion coefficient. 

IMEX
____

For the implicit-explicit advance, we use the :ref:`Chap:CdrCTU` discretization to center the divergence at the half time step.
We seek the update

.. math::
   
   \frac{\phi^{k+1} - \phi^k}{\Delta t} - \frac{1}{2}\left[\nabla\left(D\nabla\phi^k\right) + \nabla\left(D\nabla\phi^{k+1}\right)\right] = -\nabla\cdot\mathbf{v}\phi^{k+1/2},

which is a Crank-Nicholson discretization of the diffusion equation with a source term centered on :math:`k+1/2`.
We use the :ref:`Chap:CdrCTU` for obtaining the edge states :math:`\phi^{k+1/2}` and then complete the update by solving the corresponding Helmholtz equation

.. math::

   \phi^{k+1} - \frac{\Delta t}{2}\nabla\left(D\nabla\phi^{k+1}\right) = \phi^k - \Delta t\nabla\cdot\mathbf{v}\phi^{k+1/2} + \frac{\Delta t}{2}\nabla\left(D\nabla\phi^k\right).

In this case the time step limitation is

.. math::

   \Delta t \leq \frac{\Delta x}{\sum_i\left|v_i\right|}.

.. warning::

   It is possible to use this module with any implementation of ``CdrSolver``, but the IMEX discretization only makes sense if the hyperbolic term can be centered on :math:`k+1/2`.

If the truncation order of :math:`\phi^{k+1/2}` is :math:`\mathcal{O}\left(\Delta t^2\right)`, the resulting IMEX discretization is globally second order accurate.
For the :ref:`Chap:CdrCTU` discretization the edge states are accurate to :math:`\mathcal{O}\left(\Delta t\Delta x\right)`, so the scheme is globally first order convergent (but with a small error constant).
		

Initial data
------------

Default behavior
________________

By default, the initial data for this problem is given by a super-Gaussian blob

.. math::

   \phi\left(\mathbf{x},t=0\right) = \phi_0\exp\left(-\frac{\left|\mathbf{x}-\mathbf{x}_0\right|^4}{2R^4}\right),

where :math:`\phi_0` is an amplitude, :math:`\mathbf{x}_0` is the blob center and :math:`R` is the blob radius.
These are set by the input options

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.options
   :language: text
   :lines: 11-13

Custom value
____________

For a more general way of specifying initial data, ``AdvectionDiffusionStepper`` has a public member function

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.H
   :language: c++
   :dedent: 6
   :lines: 209-214

Velocity field
--------------

Default behavior
________________

The default velocity field for this class is

.. math::

   v_x &= -r\omega\sin\theta, \\
   v_y &=  r\omega\cos\theta, \\
   v_z &= 0,

where :math:`r = \sqrt{x^2 + y^2}`, :math:`\tan\theta = \frac{x}{y}`.
I.e. the flow field is a circulation around the Cartesian grid origin.

To adjust the velocity field through :math:`\omega`, simply set the following option:

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.options
   :language: text
   :lines: 14

Custom value
____________

For a more general way of setting a user-specified velocity, ``AdvectionDiffusionStepper`` has a public member function

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.H
   :language: c++
   :dedent: 6
   :lines: 216-221


Diffusion coefficient
---------------------

Default behavior
________________

The default diffusion coefficient for this problem is set to a constant.
To adjust it,  :math:`\omega`, set

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.options
   :language: text
   :lines: 15

to a chosen value.

Custom value
____________

For a more general way of setting the diffusion coefficient, ``AdvectionDiffusionStepper`` has a public member function

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.H
   :language: c++
   :dedent: 6
   :lines: 223-228

Boundary conditions
-------------------

At the EB, this module uses a wall boundary condition (i.e. no flux into or out of the EB).
On domain edges (faces in 3D), the user can select between wall boundary conditions or outflow boundary conditions by selecting the corresponding input option for the solver.
See :ref:`Chap:CdrCTU`.

Cell refinement
----------------

The cell refinement is based on two criteria:

#. The amplitude of :math:`\phi`.
#. The local curvature :math:`\left|\nabla\phi\right|\Delta x/\phi`.

We refine if the curvature is above some threshold :math:`\epsilon_1` *and* the amplitude is above some threshold :math:`\epsilon_2`.
These can be adjusted through

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.options
   :language: text
   :lines: 23-27

By default, every cell that fails to meet the above two criteria is coarsened.

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/AdvectionDiffusion` is the simplest way.
A full description is available in the ``README.md`` contained in the folder:

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   python setup.py --help

.. _Chap:AdvectionDiffusionOptions:

Solver configuration
--------------------

The ``AdvectionDiffusionStepper`` and ``AdvectionDiffusionTagger`` classes come with user-configurable input options that can be adjusted at runtime.
The configuration options for ``AdvectionDiffusionStepper`` are given below:

.. literalinclude:: ../../../../Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.options
   :language: text

Example programs
----------------

Some example programs for this module are given in

* :file:`$DISCHARGE_HOME/Exec/Examples/AdvectionDiffusion/DiagonalFlowNoEB`
* :file:`$DISCHARGE_HOME/Exec/Examples/AdvectionDiffusion/PipeFlow`

Verification
------------

Spatial and temporal convergence tests for this module (and thus also the underlying solver implementation) are given in

* :file:`$DISCHARGE_HOME/Exec/Convergence/AdvectionDiffusion/C1`
* :file:`$DISCHARGE_HOME/Exec/Convergence/AdvectionDiffusion/C2`  

C1: Spatial convergence
_______________________

A spatial convergence test is given in :file:`$DISCHARGE_HOME/Exec/Convergence/AdvectionDiffusion/C1`.
The problem solves for an advected and diffused scalar in a rotational velocity in the presence of an EB:

.. _Fig:AdvectionDiffusionC1_1:
.. figure:: /_static/figures/AdvectionDiffusionC1_1.png
   :width: 40%
   :align: center

   Final state on a :math:`512^2` uniform grid. 

To compute the convergence rate we compute two solutions with grid spacings :math:`\Delta x` and :math:`\Delta x/2`, and estimate the solution error using the approach in :ref:`Chap:SpatialConvergence`. 
Figure :numref:`Fig:AdvectionDiffusionC1_2` shows the computed convergence rates with various choice of slope limiters.
We find 2nd order convergence in all three norms for sufficiently fine grid when using slope limiters, and first order convergence when limiters are turned off.
The reduced convergence rates at coarser grids occur due to 1) insufficient resolution of the initial density profile and 2) under-resolution of the geometry. 

.. _Fig:AdvectionDiffusionC1_2:
.. figure:: /_static/figures/AdvectionDiffusionC1_2.png
   :width: 70%
   :align: center

   Spatial convergence rates with various limiters. 

C2: Temporal convergence
________________________

A temporal convergence test is given in :file:`$DISCHARGE_HOME/Exec/Convergence/AdvectionDiffusion/C2`.
To compute the temporal convergence rate we compute two solutions using time steps :math:`\Delta t` and :math:`\Delta t/2`, and estimate the solution error using the approach in :ref:`Chap:TemporalConvergence`. 
Figure :numref:`Fig:AdvectionDiffusionC2` shows the computed convergence rates for the second order Runge-Kutta and the IMEX discretization. 
As expected, we find 2nd order convergence for Heun's method and first order convergence for the IMEX discretization. 

.. _Fig:AdvectionDiffusionC2:
.. figure:: /_static/figures/AdvectionDiffusionC2.png
   :width: 70%
   :align: center

   Temporal convergence rates. 
