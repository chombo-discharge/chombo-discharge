.. _Chap:BrownianWalkerModel:

Brownian walker
***************

The Brownian walker model runs a single microscropic drift-diffusion using the :ref:`Chap:ItoSolver`, solving

.. math::
   
   \Delta\mathbf{X} = \mathbf{V}\Delta t + \sqrt{2D\Delta t}\mathbf{W}

where :math:`\mathbf{X}` is the spatial position of a particle :math:`\mathbf{V}` is the drift velocity and :math:`D` is the diffusion coefficient *in the continuum limit*.   

.. tip::

   Source code is located in :file:`$DISCHARGE_HOME/Physics/BrownianWalker`.
   
The model consists of the following implementation files:

* :file:`CD_BrownianWalkerStepper.H/cpp` which implements the integration routines. 
* :file:`CD_BrownianWalkerSpecies.H/cpp` which implements initial conditions.
* :file:`CD_BrownianWalkerTagger.H/cpp` which implements mesh refinement and de-refinement criteria.

Solvers
-------

This application uses the following solvers:

* :ref:`Chap:ItoSolver`.




