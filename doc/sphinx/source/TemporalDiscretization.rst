.. _Chap:TemporalDiscretization:

Temporal discretization
======================

In this chapter we discuss the most commonly used temporal integrators for ``PlasmaC``, and discuss their input parameters. 

godunov
-------
The ``rk2_tga`` temporal integrator is our least sophisticated temporal integrator. ``rk2_tga`` uses a Godunov splitting between advection-reaction and diffusion, and advances the equations of motion as follows:

strang2
-------
The ``strang2`` temporal integrator uses a second order Strang splitting between advection-reaction and diffusion, and supports up to fifth order efficient Runge-Kutta schemes for the advection-reaction part. 

SISDC
-----
``SISDC`` is a semi-implicit spectral deferred correction method for the ``PlasmaC`` equation set and is an adaptive high-order discretization with implicit diffusion.

MISDC
-----
``MISDC`` is a semi-implicit spectral deferred correction method for the ``PlasmaC`` equation set and is an adaptive high-order discretization with implicit diffusion. 
