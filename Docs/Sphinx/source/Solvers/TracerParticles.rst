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
   :lines: 23-37
   :dedent: 0

where ``P`` is the payload type used for the solver.
The particles are stored Struct-of-Arrays in a ``ParticleContainerSoA``: the position and weight are
container-owned columns, while ``P`` declares the *extra* payload columns. The template constraints on
``P`` are

#. It *must* expose the velocity as ``SpaceDim`` scalar payload columns named ``D_DECL(vx, vy, vz)``.
#. It *must* have a ``ParticleTraits<P>`` specialization listing its columns (including ``vx/vy/vz``).

Users are free to provide their own payload provided that it meets these constraints.
However, we also define a plug-and-play payload for the tracer-particle stepper, see :ref:`Chap:TracerParticlePayload`.

.. note::

   The ``TracerParticleSolver<P>`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTracerParticleSolver.html>`_.

.. _Chap:TracerParticlePayload:

TracerParticlePayload
---------------------

``TracerParticlePayload`` is the plug-and-play payload used by the tracer-particle stepper. It declares the
interpolated velocity together with Runge-Kutta stage scratch, all as per-component SoA columns:

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticlePayload.H
   :language: c++
   :lines: 26-33
   :dedent: 0

The velocity columns ``D_DECL(vx, vy, vz)`` are the ones required by the solver; the remaining columns hold
intermediate Runge-Kutta integration state. Weight and position are container-owned and are therefore *not*
payload members.

Initialization
--------------

To initialize the solver, one can use the full constructor

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 58-63
   :dedent: 2

Getting the particles
---------------------

To obtain the solver particles, simply call

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 282-294
   :dedent: 2

This returns the ``ParticleContainerSoA<P>`` holding the particles. 

Setting :math:`\mathbf{v}`
--------------------------

To set the velocity field on the mesh, use

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 148-153
   :dedent: 2

This will associate the input velocity ``a_velocityField`` with :math:`\mathbf{v}`.

Interpolating velocities
------------------------

To compute :math:`\mathbf{V} = \mathbf{v}\left(\mathbf{X}\right)` for all particles that reside in the solver, use

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 201-205
   :dedent: 2

This will interpolate the velocities to the particle positions using the user-defined interpolation method (see :ref:`Chap:TracerInputOptions`).

If desirable, one can also interpolate a scalar field defined on the mesh onto the particle weight by calling

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 194-199
   :dedent: 2

The interpolation function is set by the user, see :ref:`Chap:TracerInputOptions`.
See :ref:`Chap:ParticleMesh` for further details.

Deposit particles
-----------------

To deposit the particles, call

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.H
   :language: c++
   :lines: 187-192
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
