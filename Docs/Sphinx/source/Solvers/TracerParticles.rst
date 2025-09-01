.. _Chap:TracerParticles:

Tracer particles
================

Tracer particles are particles that move along a prescribed velocity field

.. math::

   \frac{\partial\mathbf{X}}{\partial t} = \mathbf{V}

where :math:`\mathbf{X}` is the particle position and :math:`\mathbf{V}` is the particle velocity.
The velocity is interpolated from a mesh-based field as

.. math::

   \mathbf{V} = \mathbf{v}\left(\mathbf{X}\right),

where :math:`\mathbf{v}` is a velocity field defined on the mesh. 
Such particles are useful, for example, for numerical integration along field lines.

.. tip::
   
   The ``chombo-discharge`` tracer particle functionality resides in :file:`$DISCHARGE_HOME/Source/TracerParticles`.

.. _Chap:TracerParticleSolver:

TracerParticleSolver
--------------------

The tracer particle solver is templated as

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 22-36
   :dedent: 0

where ``P`` is the particle type used for the solver.
The template constraints on ``P`` are

#. It *must* contain a function ``RealVect& position()``
#. It *must* contain a function ``const Real& weight() const``
#. It *must* contain a function ``RealVect& velocity()``.

Users are free to provide their own particle type provided that it meets these template constraints.
However, we also define a plug-and-play particle class that meets these requirements, see :ref:`Chap:TracerParticle`.

.. note::

   The ``TracerParticleSolver<P>`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTracerParticleSolver.html>`_.

.. _Chap:TracerParticle:

TracerParticle
--------------

The ``TracerParticle`` type inherits from :ref:`Chap:GenericParticle` particle class and is templated as

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticle.H
   :language: c++
   :lines: 25-32
   :dedent: 0

The class also defines two more members; the weight and a particle velocity.
These are accessible as

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticle.H
   :language: c++
   :lines: 51-77
   :dedent: 2

Note that, just as for ``GenericParticle``, the template arguments ``M`` and ``N`` indicates the number of scalars and vectors allocated to the particle, see :ref:`Chap:GenericParticle`.
These data fields can be used by applications for, e.g., storing integration variables (such as intermediate positions in a Runge-Kutta code).

Initialization
--------------

To initialize the solver, one can use the full constructor

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 56-61
   :dedent: 2

Getting the particles
---------------------

To obtain the solver particles, simply call

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 278-290
   :dedent: 2

This returns the ``ParticleContainer<P>`` holding the particles, see :ref:`Chap:ParticleContainer`. 

Setting :math:`\mathbf{v}`
--------------------------

To set the velocity field on the mesh, use

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 146-151
   :dedent: 2

This will associate the input velocity ``a_velocityField`` with :math:`\mathbf{v}`.

Interpolating velocities
------------------------

To compute :math:`\mathbf{V} = \mathbf{v}\left(\mathbf{X}\right)` for all particles that reside in the solver, use

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 198-202
   :dedent: 2

This will interpolate the velocities to the particle positions using the user-defined interpolation method (see :ref:`Chap:TracerInputOptions`).

If desirable, one can also interpolate a scalar field defined on the mesh onto the particle weight by calling

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 192-196
   :dedent: 2

The interpolation function is set by the user, see :ref:`Chap:TracerInputOptions`.
See :ref:`Chap:ParticleMesh` for further details.

Deposit particles
-----------------

To deposit the particles, call

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 185-190
   :dedent: 2

This will deposit the particle weights onto the input data holder.

The deposition function is set by the user, see :ref:`Chap:TracerInputOptions`.
Complete details regarding how the deposition functions work is available in :ref:`Chap:ParticleMesh`.

.. _Chap:TracerInputOptions:

Input options
-------------

Available input options for the tracer particle solver are given in the listing below.

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.options
   :caption: List on configuration options for ``TracerParticleSolver<P>``.
	     All options are run-time configurable.
