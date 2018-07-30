.. _Chap:time_stepper:

time_stepper
------------

The :ref:`Chap:time_stepper` class handles the integration of the plasma equations. :ref:`Chap:time_stepper` is an abstract for which we have several implements (such as a second order Runge Kutta method). Writing new temporal integrators is an extensive task, but the base class :ref:`Chap:time_stepper` contains a lot of short hand functionality for updating equations, and also contains an interface to :ref:`Chap:plasma_kinetics`. Since :ref:`Chap:time_stepper` does not actually contain an advance method for the equations of motion, we recommend that the user refers to the API of the temporal integrator that he uses (see :ref:`Chap:Solvers`) for a full explanation of the various integrators. However, the following options for :ref:`Chap:time_stepper` are passed into *all* implementation classes:

.. literalinclude:: links/time_stepper.options

The input options above are, for the most part, self-explanatory. Mostly, they refer to the handling of the size of the time step, for example by passing the Courant-Friedrichs-Lewy number, or setting a minimum or maximum possible time step. However, all of these options *may* be handled differently by different integrators, since different schemes have different restrictions on stable time steps.

Finally, there is an option to allow radiative transport updates only at certain time steps by modifying (at his own peril) the **fast_rte** flag. Yet again, we remark that this flag may be handled differently by different solvers. 
