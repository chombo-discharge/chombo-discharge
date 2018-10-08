.. _Chap:Integrators:

Time steppers
-------------

We have various implementation of :ref:`Chap:time_stepper` that allow different temporal integration of the equations of motion. For the most part, we recommend using the second order Runge Kutta implementation with implicit diffusion. However, please refer to the :doxy:`Doxygen API <classtime__stepper>` to see all available integrators. 

Typically, time steppers are selected at compile time following the style in :ref:`Chap:NewSimulations`. However, the user may select time steppers at run-time by modifying his main file. For each time stepper, there are various options available at run-time through an input script.

Multirate SSP Runge-Kutta
_________________________

The multirate solver is the one that we use the most. This solver uses an adaptive-in-time method (without error control) which removes redundant elliptic solves at CFL bound time scales. This implies that the electric field and radiative transfer equations are updated only at the end of each time step. For SSP reasons, we only support first to third order Runge-Kutta stages. The tableus for these can be found in paper *On High Order Strong Stability Preserving Runge-Kutta and Multi Step Time Discretizations* by Sigal Gottlieb. The diffusion advance is done by using the TGA scheme. The following options are available:

.. literalinclude:: links/multirate_rkSSP.options

Second order Runge-Kutta
________________________

For the second order Runge-Kutta method, the following options are available:

.. literalinclude:: links/rk2_tga.options

In the above option, ``alpha`` contains the Butcher tableu. Note that only Heun's method with ``alpha`` equal to one, is monotone in time (a.k.a. strongly stability preserving) with a CFL of one. The above scheme solves the equations of motion by using a consistent second-order accurate Runge-Kutta scheme where the elliptic solves (RTE and Poisson) are synchronized at every stage. 

Other integrators
_________________

In this user guide, we do not discuss the other temporal integrators for two reasons. 1) These integrators are not sufficiently stable for full release and 2) We have not performed benchmarking or comparisons of these solvers. However, you may of course use these solvers nonetheless (at your own peril). 
