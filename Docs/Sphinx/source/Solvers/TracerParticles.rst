.. _Chap:TracerParticles:

Tracer particles
================

Tracer particles are particles that move along a velocity field

.. math::

   \frac{\partial\mathbf{X}}{\partial t} = \mathbf{V}

where :math:`\mathbf{V}` is the particle velocity.
This is interpolated from a mesh-based field as

.. math::

   \mathbf{V} = \mathbf{v}\left(\mathbf{X}\right),

where :math:`\mathbf{v}` is a velocity field defined on the mesh. 
Such particles are useful, for example, for numerical integration along field lines.

``chombo-discharge`` defines AMR-ready tracer particle solvers in :file:`$DISCHARGE_HOME/Source/TracerParticles`.

.. _Chap:TracerParticleSolver:

TracerParticleSolver
--------------------

The tracer particle solver is templated as

.. code-block:: c++

   template <typename P>
   class TracerParticleSolver;

where ``P`` is the particle type used for the solver.
The template constraints on ``P`` are

#. It *must* contain a function ``RealVect& position()``
#. It *must* contain a function ``const Real& weight() const``
#. It *must* contain a function ``RealVect& velocity()``.

Users are free to provide their own particle type provided that it meets these template constraints.
However, we also define a plug-and-play particle class that meets these requirements, see :ref:`Chap:TracerParticle`.

.. note::

   The ``TracerParticleSolver`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTracerParticleSolver.html>`_.

.. _Chap:TracerParticle:

TracerParticle
--------------

The ``TracerParticle`` type inherits from the ``GenericParticle`` particle class and is templated as

.. code-block:: c++

   template <size_t M, size_t N>
   class TracerParticle<M,N> : public GenericParticle<M, N>

and also defines two more members: A particle weight and a particle velocity.
These are accesible as

.. code-block:: c++

   template <size_t M, size_t N>
   Real&
   TracerParticle<M, N>::weight();

   template <size_t M, size_t N>
   RealVect&
   TracerParticle<M, N>::velocity();

Note that, just as for ``GenericParticle``, the template arguments ``M`` and ``N`` indicates the number of scalars and vectors allocated to the particle.
See :ref:`Chap:GenericParticle`.

.. note::

   The ``TracerParticleSolver`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classTracerParticle.html>`_.

Initialization
--------------

To initialize the solver, one can use the full constructor

.. code-block:: c++

   template <typename P>
   TracerParticleSolver<P>::TracerParticleSolver(const RefCountedPtr<AmrMesh>& a_amr,
		                                 const RefCountedPtr<ComputationalGeometry> a_compGeom);

Getting the particles
---------------------

To obtain the solver particles simply call

.. code-block:: c++

   template <typename P>
   ParticleContainer<P>&
   TracerParticleSolver<P>::getParticles();

This returns the ``ParticleContainer`` holding the particles, see :ref:`Chap:ParticleContainer`. 

Setting :math:`\mathbf{v}`
--------------------------

To set the velocity use

.. code-block:: c++

   template <typename P>
   void
   TracerParticleSolver<P>::setVelocity(const EBAMRCellData& a_velocityField)

This will associate the input velocity ``a_velocityField`` with :math:`\mathbf{v}`.

Interpolating velocities
------------------------

To compute :math:`\mathbf{V} = \mathbf{v}\left(\mathbf{X}\right)` use

.. code-block:: c++

   template <typename P>
   void
   TracerParticleSolver<P>::interpolateVelocities();

This will interpolate the velocities to the particle positions using the user-defined interpolation method (see :ref:`Chap:TracerInputOptions`).

One can also interpolate a scalar field defined on the mesh onto the particle weight by calling

.. code-block:: c++

   template <typename P>
   void
   TracerParticleSolver<P>::interpolateWeight(const EBAMRCellData& a_scalar) noexcept;

Letting :math:`f` define the input field ``a_scalar``, this performs the operation :math:`w = f\left(\mathbf{X}\right)`.

Deposit particles
-----------------

To deposit the particles call

.. code-block:: c++

   template <typename P>
   void
   TracerParticleSolver<P>::deposit(EBAMRCellData& a_phi)

This will deposit the particle weights onto the input data holder. 

.. _Chap:TracerInputOptions:

Input options
-------------

Available input options for the tracer particle solver are

.. literalinclude:: ../../../../Source/TracerParticles/CD_TracerParticleSolver.options

The flags ``deposition`` and ``interpolation`` indicates which deposition and interpolation methods will be used.
Likewise, ``deposition_cf`` indicates the coarse-fine deposition strategy, see :ref:`Chap:Particles`.
The flags ``plot_weight`` and ``plot_velocity`` indicates whether or not to include the particle weights and velocities in plot files. 

