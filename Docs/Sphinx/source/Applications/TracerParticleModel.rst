.. _Chap:TracerParticleModel:

Tracer particle model
=====================

The tracer particle model was written in order to illustrate the usage of :ref:`Chap:TracerParticleSolver`.
The module can trace particles passively in a fixed velocity fieldk.`).

.. tip::

   The source code is located in :file:`$DISCHARGE_HOME/Physics/TracerParticle`.

The model consists of the following implementation files:

* :file:`CD_TracerParticleStepper.H`, which implements :ref:`Chap:TimeStepper` as a subclass ``TracerParticleStepper``.

.. note::

   The current implementation does not incorporate solver-based adaptive mesh refinement, so refinement is restricted to refinement and coarsening through the :ref:`Chap:Driver` interface.

Basic problem
-------------

Currently, ``TracerParticleStepper`` simply instantiates a :ref:`Chap:TracerParticleSolver` the particles in a fixed velocity field.
There is currently no handling of particles that fall off the domain.

Initial conditions
__________________

The user can set the velocity field, which is either a velocity field that is diagonal to the Cartesian mesh, or a rotational flow. 
Likewise, the user can set the number of particles, which are sampled uniformly across the domain.

Integration algorithm
_____________________

``TracerParticleStepper`` implements the following algorithms for the particle advection: The explicit Euler method, Heun's method, and the classical fourth order Runge-Kutta method.

Solver configuration
--------------------

The ``TracerParticleStepper`` class comes with user-configurable input options that can be adjusted at runtime.
These configuration options are given below.

.. literalinclude:: ../../../../Physics/TracerParticle/CD_TracerParticleStepper.options
   :language: text

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/TracerParticle` is the simplest way.
A full description is available in the ``README.md`` contained in the folder:

.. literalinclude:: ../../../../Physics/TracerParticle/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   ./setup.py --help
