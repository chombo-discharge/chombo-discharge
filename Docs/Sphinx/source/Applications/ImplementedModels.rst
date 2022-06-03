.. _Chap:ImplementedModels:

Implemented models
==================

Various models that use the ``chombo-discharge`` solvers and functionality have been implemented and are shipped with ``chombo-discharge``.
These models implement ``TimeStepper`` for advancing various types of equations, and have been set up both for unit testing and as building blocks for more complex applications. 
The models reside :file:`$DICHARGE_HOME/Physics` and are supplemented with Python tools for quickly setting up new applications that use the same code.
In general, users are encouraged to modify or copy these models and tailor them for their own applications.

The following models are currently supported:

#. :ref:`Chap:AdvectionDiffusionModel` Used for solving pure advection-diffusion problems in 2D/3D with complex geometries.
#. :ref:`Chap:BrownianWalkerModel` Used for solving advection-diffusion problems using computational particles.
   This is, essentially, a Particle-In-Cell code that uses classical particles (rather than kinetic ones).
#. :ref:`Chap:CdrPlasmaModel` Implements ``TimeStepper`` for solving low-temperature discharge plasma problems using fluids.
   Used e.g. for streamer simulations.
#. :ref:`Chap:ElectrostaticsModel` Implements ``TimeStepper`` for setting up a static calculation that solves the Poisson equation in 2D/3D.   
#. :ref:`Chap:GeometryModel` Implement ``TimeStepper`` with empty functionality. Often used when setting up new geometries/cases.
#. :ref:`Chap:RadiativeTransferModel` Implements ``TimeStepper`` for solving radiative transfer problems. Used e.g. for regression testing.

