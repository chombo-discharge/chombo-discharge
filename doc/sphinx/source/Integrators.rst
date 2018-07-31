.. _Chap:Integrators:

Time steppers
-------------

We have various implementation of :ref:`Chap:time_stepper` that allow different temporal integration of the equations of motion. For the most part, we recommend using the second order Runge Kutta implementation, which we believe is bug-free. However, please refer to the :doxy:`Doxygen API <classtime__stepper>` to see all available integrators. 

Time steppers are selected at compile time following the style in :ref:`Chap:NewSimulations`. Typically, there are various options available at run-time through an input script.

Second order Runge-Kutta
________________________

For the second order Runge-Kutta method, the following options are available:

.. literalinclude:: links/rk2.options

In the above option, **alpha** contains the Butcher tableu. Note that only Heun's method with **alpha** equal to one, is monotone in time (a.k.a. strongly stability preserving). 

Other integrators
_________________

In this user guide, we do not discuss the other temporal integrators for two reasons. 1) These integrators are not sufficiently stable for full release and 2) We have not performed benchmarking or comparisons of these solvers. However, you may of course use these solvers nonetheless (at your own peril). 
