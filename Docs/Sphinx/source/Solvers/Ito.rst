.. _Chap:ItoDiffusion:

Îto diffusion
=============

The Îto diffusion model advances computational particles as drifting Brownian walkers

.. math::
   :label: ito_diffusion
	   
   \Delta\mathbf{X} = \mathbf{V}\Delta t + \sqrt{2D\Delta t}\mathbf{W}

where :math:`\mathbf{X}` is the spatial position of a particle, :math:`\mathbf{V}` the particle drift velocity, and :math:`D` is the diffusion coefficient *in the continuum limit*.
The vector term :math:`\mathbf{W}` indicates a random number sampled from a Gaussian distribution with mean value of 0 and standard deviation of 1.

.. tip::
   
   The code for Îto diffusion is given in :file:`/Source/ItoDiffusion`.

ItoParticle
-----------

The ``ItoParticle`` is used as the underlying particle type for running the Ito drift-diffusion solvers.
It derives from :ref:`Chap:GenericParticle` as follows:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoParticle.H
   :language: c++
   :lines: 39
   :dedent: 0

From the signature one can see that ``ItoParticle`` contains a number of extra class ``Real`` and ``RealVect`` class members.
These extra fields are used for storing the following information in the particle:

#. Particle weight, mobility, diffusion coefficient, energy (not currently used), and a holder for a scratch storage. 
#. The previous particle position, the velocity, and a holder for a ``RealVect`` scratch storage.

.. tip::
   
   Several member functions are available for obtaining the particle properties. See the full `ItoParticle C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoParticle.html>`_

.. _Chap:ItoSolver:

ItoSolver
---------

The ``ItoSolver`` class encapsulates the implementation of :eq:`ito_diffusion` in ``chombo-discharge``.
This class can advance a set of computational particles (see :ref:`Chap:ItoParticle`) with the following functionality:

#. Move particles the a microscopic drift-diffusion model.
#. Compute particle intersection with embedded boundaries and domain edges.
#. Deposit particles and other particle types on the mesh.
#. Interpolate velocities and diffusion coefficients to the particle positons.
#. Manage superparticle splitting and merging.

Internally, ``ItoSolver`` stores its particles in various ``ParticleContainer<ItoParticle>`` containers.
Although the particle velocities and diffusion coefficients can be manually assigned, they can also be interpolated from the mesh.
``ItoSolver`` stores the following properties on the mesh:

#. Mobility.
#. Diffusion coefficient.
#. Velocity function.

The reason for storing both the mobility and velocity function is to simply to improve flexibility when assigned the particle velocity :math:`\mathbf{V}`.
Note that the velocity function does *not* have to represent the particle velocity.
When using both the mobility and velocity function, one can compute the particle velocity as :math:`\mathbf{V} = \mu\mathbf{v}`, where :math:`\mathbf{v}` is a velocity field.
This is typically done for discharge simulations where for simplicity we assign :math:`\mathbf{v}` to be the electric field, and :math:`\mu` to the the field-dependent mobility.
Additional information is available in :ref:`Chap:ItoInterpolation`.

.. _Chap:ItoSpecies:

ItoSpecies
-----------

``ItoSpecies`` is a class for parsing solver information into ``ItoSolver``, e.g., whether or not the particle type is mobile or not.
The constructor for the ``ItoSpecies`` class is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 35-42
   :dedent: 2

Here, ``a_name`` indicates a variable name for the solver.
This variable will be used in, e.g., error messages and I/O functionality.
``a_chargeNumber`` indicates the charge number of the species and the two booleans ``a_mobile`` and ``a_diffusive`` indicates whether or not the solver is mobile or diffusive.

.. note::

   The C++ ``ItoSpecies`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoSpecies.html>`_.

Supplying initial data
______________________

Initial data for the ``ItoSolver`` is provided through ``ItoSpecies`` by providing it with the following:

#. Initial particles specified from a list (``List<ItoParticle>``) of particles.
#. Provide a density description from which initial particles are stochastically sampled within each grid cell.

In particular, there are two data members that must be populated:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 147-155
   :dedent: 2

These can either be populated during construction, or explicitly supplied via the following set functions:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 100-112
   :dedent: 2

When ``ItoSolver`` initializes the data in the solver, it will copy the particle list ``m_initialParticles`` from the species and into the solver.

.. tip::

   When using MPI, the user must ensure that each MPI rank does not provide duplicate particles.
   The ``ParticleOps`` class contains lots of supporting functionality for sampling particles with MPI, see the `ParticleOps C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classParticleOps.html>`_

When sampling particles from a mesh-based density, the solver will generate the particles so that the specified density is approximately reached within each grid cell.
If the density that is supplied does not lead to an integer number of particles in the grid cell (which is virtually always the case), the evaluation of the number of particles is stochastically evaluated.
E.g., if the density is :math:`\phi` and then grid cell volume is :math:`\Delta V`, and :math:`\phi\Delta V = 1.2`, then there is a 20% chance that there will be generated two particles within the grid cell, and 80% chance that only one particle will be generated.

.. warning::
   
   There is currently a hard limit that restricts the number of initial computational particles per cell to 32. 

Particle containers
-------------------

Internally, ``ItoSolver`` contains several ``ParticleContainer<ItoParticle>`` for storing various categories of particles.
These categories exist because the transport kernel will almost always lead to particles that leave the domain or intersect the EB.
Chemistry models that use ``ItoSolver`` for tracking particles might also require *new* particles to be added into the domain.

``ItoSolver`` defines an enum ``WhichContainer`` for classification of ``ParticleContainer<ItoParticle>`` data holders for holding particles that live on:

* Main particles (``WhichContainer::Bulk``). 
* The embedded boundary (``WhichContainer::EB``).
* On the domain edges/faces (``WhichContainer::Domain``).
* Representing ''source particles'' (``WhichContainer::Source``).
* Particles that live *inside* the EB (``WhichContainer::Covered``).

The particles are available from the solver through the function

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 635-640
   :dedent: 2

Usually, ``ItoSolver`` will perform a drift-diffusion advance and the user will then check if some of the particles crossed into the EB.
The solver can then automatically fill the boundary particles containers, see :ref:`Chap:ParticleIntersection`.

Computing the particle velocity
-------------------------------

For the ``ItoSolver`` the particle velocity is computed as

.. math::

   \mathbf{V} = \mu\left(\mathbf{X}\right)\mathbf{v}\left(\mathbf{X}\right)

where :math:`\mu` is a particle mobility and :math:`\mathbf{v}` is a velocity field defined on the mesh.
Note that both :math:`\mu` and :math:`\mathbf{v}` are defined on the mesh.
The solver can, alternatively, also compute the velocity as

.. math::

   \mathbf{V} = \left(\mu\mathbf{v}\right)\left(\mathbf{X}\right),

i.e. through interpolation of :math:`\mu\mathbf{v}` to the particle position.
Regardless of which method is chosen (see :ref:`Chap:ItoInterpolation`), both :math:`\mu` and :math:`\mathbf{v}` exist on the mesh (stored as ``EBAMRCellData``).



Transport kernel
----------------

The transport kernels for the ``ItoSolver`` will simply consist of particle updates of the following type:

.. code-block:: c++

   Real a_dt;
   List<ItoParticle> particles;
   
   for (ListIterator<ItoParticle>& lit(particles); lit.ok(); ++lit) {
      ItoParticle& p = lit();

      p.oldPosition() = p.position();
      p.position()   += p.velocity() * a_dt + sqrt(2.0*p.diffusion()*a_dt) * this->randomGaussian();
   }

The function ``randomGaussian`` implements a diffusion hopping and returns a 2D/3D dimensional vector with values drawn from a normal distribution with standard width of one and mean value of zero.
The implementation uses the random number generators in :ref:`Chap:Random`.
The user can choose to truncate the normal distribution, see :ref:`Chap:ItoInput`.

Remapping particles
-------------------

To remap, call the ``ItoSolver`` remapping functions as

.. code-block:: c++

   void
   ItoSolver::remap();

This will remap the particles to the correct grid patches.

To remap the other ``ParticleContainer`` data holders (holding e.g. particles that intersected the EB), there's an alternative function

.. code-block:: c++

   void
   ItoSolver::remap(const WhichContainer a_container);

where ``a_container`` is one of the particle containers.

Deposition
----------

To deposit the particle weights onto the grid one can use

.. code-block:: c++

   void
   ItoSolver::depositParticles();

The particles are deposited into the data holder ``m_phi``.
The data can then be fetched with

.. code-block:: c++

   EBAMRCellData&
   ItoSolver::getPhi();
   
To deposit a different particle data holder into ``m_phi`` one can use

.. code-block:: c++

   void
   ItoSolver::depositParticles(const WhichContainer a_container);

This can be used, for example, to deposit the EB particles on the mesh.
More general methods also exist, see the ``ItoSolver`` C++ API `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoSolver.html>`_.

One can also deposit the following quantities on the mesh:

* Conductivity, which deposits :math:`\mu W`.
* Diffusivity, which deposits :math:`D W`.
* Energy, which deposits :math:`\epsilon W`.    

Here, :math:`W` is the particle weight, :math:`\mu` is the particle mobility, :math:`D` is the particle diffusion coefficient and :math:`\epsilon` is the particle energy.
These functions exist as

.. code-block:: c++

   void
   ItoSolver::depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

   void
   ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

   void
   ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

.. important::

   The ``ItoSolver`` deposition method is specified in the input script, see :ref:`Chap:ItoInput`.

.. _Chap:ItoInterpolation:

Velocity interpolation
----------------------

Computing the particle velocity is done by first computing the particle mobility and then computing the particle velocity.
For interpolating the mobility to the particle position one will call

.. code-block:: c++

   void
   ItoSolver::interpolateMobilities();

which will compute :math:`\mu\left(\mathbf{X}\right)` using the user-specified deposition/interpolation method for computing the mobility.
The solver can switch between two ways of computing the mobility.
The first is to compute :math:`\mu\left(\mathbf{X}\right)` directly.
The other is to compute the mobility as

.. math::

   \mu = \frac{\left(\mu\left|\mathbf{v}\right|\right)\left(\mathbf{X}\right)}{\left|\mathbf{v}\left(\mathbf{X}\right)\right|}.

When computing the particle velocity as :math:`\mathbf{V} = \mu\left(\mathbf{X}\right)\mathbf{v}\left(\mathbf{X}\right)`, the latter method ensures that :math:`\mathbf{V} = \left(\mu\mathbf{v}\right)\left(\mathbf{X}\right)`.

.. note::

   The user selects between the two mobility interpolation methods in the input script.
   See :ref:`Chap:ItoInput`.

After the mobility has been appropriately set, the velocity can be interpolated from

.. code-block:: c++

   void
   ItoSolver::interpolateVelocities();

The above will compute :math:`v\left(\mathbf{X}\right)` and set the velocity as :math:`\mathbf{V} = \mu\left(\mathbf{X}\right)\mathbf{v}\left(\mathbf{X}\right)`.

.. important::

   The ``ItoSolver`` interpolation method is specified in the input script, see :ref:`Chap:ItoInput`.

Particle intersections
----------------------

It will happen that particles occasionally hit the embedded boundary or leave through the domain sides.
In this case one might want to keep the particles in separate data holders rather than discard them.
``ItoSolver`` supplies the following routine for transferring the particles to the containers that hold the EB and domain particles:

.. code-block:: c++

   void
   ItoSolver::intersectParticles(const EbIntersection a_ebIntersection, const bool a_deleteParticles);

Here, ``EbIntersection`` is a just an enum for putting logic into how the intersection is computed.
Valid options are ``EbIntersection::Bisection`` and ``EbIntersection::Raycast``.
These algorithms are discussed in :ref:`Chap:ParticleEB`.
The flag ``a_deleteParticles`` specifies if the original particles should be deleted when populating the other particle containers.

After calling ``intersectParticles``, the particles that crossed the EB or domain walls are available through the ``getParticles`` routine, see :ref:`Chap:ItoSolver`. 
   

Computing time steps
--------------------

The drift time step routines are implemented such that one restricts the time step such that the fastest particle does not move more than a specified number of grid cells.
This routine is implemented as

.. code-block::

   Real
   ItoSolver::computeAdvectiveDt() const;

which returns a CFL-like condition :math:`\Delta x/\textrm{max}(v_x, v_y, v_z)` on the the various AMR levels and patches.

The signatures for the diffusion time step are similar to the ones for drift:

.. code-block:: c++

   Real
   ItoSolver::computeDiffusiveDt() const;

and this returns another CFL-like condition :math:`\Delta x^2 / (2dD)` for all the particles, where :math:`d` is the spatial dimension.
Note that there is no fundamental limitation to how far the particles can move, unless the user explicitly makes this restriction in the input options, see :ref:`Chap:ItoInput`.

A combination of the advection and diffusion time step routines also exists as

.. code-block::

   Real
   ItoSolver::computeDt() const;

This routine computes the time step

.. math::

   \Delta t = \frac{1}{\frac{\Delta x}{\textrm{max}(v_x, v_y, v_z)} + \frac{\Delta x^2}{2dD}},

Superparticles
--------------

``ItoSolver`` currently handles superparticles through kD-trees, see :ref:`Chap:SuperParticles`, re-initialization, or user-based criteria. 
The function for splitting and merging the particles is in all cases

.. code-block:: c++

   void
   ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell);

Calling this function will merge/split the particles.

The default behavior in ``ItoSolver`` is to not merge the particles, but the user can set the merging algorithm through the input script, or supply one externally.
In order to specify the merging algorithm the user must set the ``ItoSolver.merge_algorithm`` to one of the following:

* ``none`` - No particle merging/splitting is performed.
* ``equal_weight_kd`` Use a kD-tree with bounding volume hierarchies to partition and split/merge the particles. This conserves the particle center-of-mass.
* ``reinitialize`` Re-initialize the particles in each grid cell, ensuring that weights are as uniform as possible.
* ``reinitialize_bvh`` Re-initialize the particles in each node of a kD tree. Weights are as uniform as possible. 
* ``external`` Use an externally injected particle merging algorithm. In order to use this feature the user must supply one through

  .. code-block:: c++

     // Set a particle merging algorithm
     virtual void
     setParticleMerger(const std::function<void(List<ItoParticle>& a_particles, const CelInfo& a_cellinfo, const int a_numParticles)>);

  where the input function is a function which merges the input particles, possibly also taking into account geometric information in the cell.

.. tip::
   
   ``ItoSolver`` uses the kD-node implementation from :ref:`Chap:SuperParticles` and partitioners for splitting the particles into two subsets with equal weights.

.. _Chap:ItoIO:

I/O
---

.. _Chap:ItoPlot:

Plot files
__________

``ItoSolver`` can output the following variables to plot files:

* :math:`\phi`, i.e. the deposited particle weights (``ItoSolver.plt_vars = phi``)
* :math:`\mathbf{v}`, the advection field (``ItoSolver.plt_vars = vel``).
* :math:`D`, the diffusion coefficient  (``ItoSolver.plt_vars = dco``).

It can also plot the corresponding particle data holders:

* Ito particles (``ItoSolver.plt_vars = part``).
* EB particles  (``ItoSolver.plt_vars = eb_part``).
* Domain particles  (``ItoSolver.plt_vars = domain_part``).
* Source particles  (``ItoSolver.plt_vars = source_part``).

.. _Chap:ItoCheck:  

Checkpoint files
________________

When writing checkpoint files, ``ItoSolver`` can either

* Add the particles to the HDF5 file,
* Checkpoint the corresponding fluid data.

The user specifies this through the input script variable ``ItoSolver.checkpointing``, see :ref:`Chap:ItoInput`.
If checkpointing fluid data then a subsequent restart will generate a new set of particles.

.. warning::

   If writing particle checkpoint files, simulation restarts must also *read* as if the checkpoint file contains particles. 

.. _Chap:ItoInput:

Input options
-------------

I/O
___

Plot variables are specified using ``ItoSolver.plt_vars``, see :ref:`Chap:ItoPlot`).
If adding the various particle container data holders to the plot file, the deposition method for those is specified using ``ItoSolver.plot_deposition``.

If using fluid checkpointing for simulation restarts, the flag ``ItoSolver.ppc_restart`` determines the maximum number of particles that will initialized in each grid cell during a restart. 


Particle-mesh
_____________

To specify the mobility interpolation, use ``ItoSolver.mobility_interp``.
Valid options are ``direct`` and ``velocity``, see :ref:`Chap:ItoInterpolation`.

Deposition and coarse-fine deposition (see :ref:`Chap:ParticleMesh`) is controlled using the flags

* ``ItoSolver.deposition`` for the base deposition scheme.
  Valid options are ``ngp``, ``cic``, and ``tsc``.
* ``ItoSolver.deposition_cf`` for the coarse-fine deposition strategy.
  Valid options are ``interp``, ``halo``, or ``halo_ngp``.

To modify the deposition scheme in cut-cells, one can enforce NGP interpolation and deposition through

* ``ItoSolver.irr_ngp_deposition`` for enforcing NGP deposition. Valid options are ``true`` or ``false``.
* ``ItoSolver.irr_ngp_interp`` for enforcing NGP interpolation. Valid options are ``true`` or ``false``.  

Checkpoint-restart
__________________

Available input options for the ``ItoSolver`` are listed below:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.options

Example application
-------------------

An example application of usage of the ``ItoSolver`` is found in

* :file:`$DISCHARGE_HOME/Physics/BrownianWalker`, see :ref:`Chap:BrownianWalkerModel`.
