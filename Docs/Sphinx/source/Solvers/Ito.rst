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

.. _Chap:ItoParticle:

ItoParticle
-----------

The ``ItoParticle`` is used as the underlying particle type for running the Ito drift-diffusion solvers.
It is a Struct-of-Arrays payload (see :ref:`Chap:ParticleSoA`) whose columns are

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoParticle.H
   :language: c++
   :lines: 30-54

In addition to the container-owned position and weight, ``ItoParticle`` stores the payload columns above.
These extra fields are used for storing the following information in the particle:

#. Mobility, diffusion coefficient, energy (not currently used), and a holder for a scratch scalar storage.
#. The previous particle position, the velocity, and a holder for a ``RealVect`` scratch storage.

.. tip::
   
   Several member functions are available for obtaining the particle properties. See the full `ItoParticle C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoParticle.html>`_

.. _Chap:ItoSolver:

ItoSolver
---------

The ``ItoSolver`` class encapsulates the implementation of :eq:`ito_diffusion` in ``chombo-discharge``.
This class can advance a set of computational particles (see :ref:`Chap:ItoParticle`) with the following functionality:

#. Move particles with a microscopic drift-diffusion model.
#. Compute particle intersection with embedded boundaries and domain edges.
#. Deposit particles and other particle types on the mesh.
#. Interpolate velocities and diffusion coefficients to the particle positions.
#. Manage superparticle splitting and merging.

Internally, ``ItoSolver`` stores its particles in various ``ParticleContainer<ItoParticle>`` containers.
Although the particle velocities and diffusion coefficients can be manually assigned, they can also be interpolated from the mesh.
``ItoSolver`` stores the following properties on the mesh:

#. Mobility.
#. Diffusion coefficient.
#. Velocity function.

The reason for storing both the mobility and velocity function is simply to improve flexibility when assigning the particle velocity :math:`\mathbf{V}`.
Note that the velocity function does *not* have to represent the particle velocity.
When using both the mobility and velocity function, one can compute the particle velocity as :math:`\mathbf{V} = \mu\mathbf{v}`, where :math:`\mathbf{v}` is a velocity field.
This is typically done for discharge simulations where for simplicity we assign :math:`\mathbf{v}` to be the electric field, and :math:`\mu` to the field-dependent mobility.
Additional information is available in :ref:`Chap:ItoInterpolation`.

.. _Chap:ItoSpecies:

ItoSpecies
-----------

``ItoSpecies`` is a class for parsing solver information into ``ItoSolver``, e.g., whether or not the particle type is mobile or not.
The constructor for the ``ItoSpecies`` class is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 37-44
   :dedent: 2

Here, ``a_name`` indicates a variable name for the solver.
This variable will be used in, e.g., error messages and I/O functionality.
``a_chargeNumber`` indicates the charge number of the species and the two booleans ``a_mobile`` and ``a_diffusive`` indicate whether or not the solver is mobile or diffusive.

.. note::

   The C++ ``ItoSpecies`` API is available at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoSpecies.html>`_.

Supplying initial data
______________________

Initial data for the ``ItoSolver`` is provided through ``ItoSpecies`` by providing it with the following:

#. Initial particles specified from a container (``ParticleSoA<ItoParticle>``) of particles.
#. Provide a density description from which initial particles are stochastically sampled within each grid cell.

In particular, there are two data members that must be populated:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 155-163
   :dedent: 2

These can either be populated during construction, or explicitly supplied via the following set functions:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSpecies.H
   :language: c++
   :lines: 108-120
   :dedent: 2

When ``ItoSolver`` initializes the data in the solver, it will copy the particle list ``m_initialParticles`` from the species and into the solver.

.. tip::

   When using MPI, the user must ensure that each MPI rank does not provide duplicate particles.
   The ``ParticleOps`` class contains lots of supporting functionality for sampling particles with MPI, see the `ParticleOps C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classParticleOps.html>`_

When sampling particles from a mesh-based density, the solver will generate the particles so that the specified density is approximately reached within each grid cell.
If the density that is supplied does not lead to an integer number of particles in the grid cell (which is virtually always the case), the evaluation of the number of particles is stochastically evaluated.
E.g., if the density is :math:`\phi` and the grid cell volume is :math:`\Delta V`, and :math:`\phi\Delta V = 1.2`, then there is a 20% chance that there will be generated two particles within the grid cell, and 80% chance that only one particle will be generated.

.. tip::
   
   The number of initially sampled particles is set through ``ItoSolver.ppc_restart``.

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
   :lines: 691-697
   :dedent: 2

Usually, ``ItoSolver`` will perform a drift-diffusion advance and the user will then check if some of the particles crossed into the EB.
The solver can then automatically fill the boundary particles containers, see :ref:`Chap:ParticleIntersection`.

Remapping particles
-------------------

``ItoSolver`` has two functions for remapping particles:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 956-967
   :dedent: 2

The bottom function lets the user remap any ``ParticleContainer<ItoParticle>`` that lives in the solver.
Here, ``a_container`` indicates which particle container to remap.

Particle deposition
-------------------

``ItoSolver`` contains several member functions for depositing various particle properties onto the mesh.
The most general version is given below:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 386-403
   :dedent: 2

This version permits the user to deposit an arbitrary per-particle quantity from a particle container ``a_particles`` onto some pre-allocated mesh storage ``a_phi``.
The quantity to be deposited is supplied through the ``a_gather`` callable, which returns a ``Real`` value for each particle in the SoA container.

.. important::

   The ``ItoSolver`` deposition methods are specified in the input script, see :ref:`Chap:ItoInput`.
   Both the base deposition scheme (e.g., NGP or CIC) must be specified, as well as the handling near refinement boundaries. 

A simpler version that deposits the bulk particles as a density on the mesh is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 317-325
   :dedent: 2

The particles are deposited into the class member ``m_phi``, which stores the particle density on the mesh. 
This data can then be fetched with

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 714-719
   :dedent: 2
   
For the full list of available deposition functions, see the ``ItoSolver`` C++ API `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoSolver.html>`_.

AMR synchronization after a deposit
___________________________________

The deposition functions put the particles on the mesh and, if redistribution is enabled, redistribute the cut-cell mass.
They do *not* average the result down onto the coarser levels, and they do not fill its ghost cells.
Whether that is wanted depends on what happens next, which only the caller knows: a deposit that is the final state of a field needs it, whereas a deposit that is merely one term of a larger sum does not, because only the assembled sum needs synchronizing.
Callers that need it should therefore call

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 353-362
   :dedent: 2

after depositing.
Note that ``ItoSolver::depositParticles`` already does this internally, so ``m_phi`` (i.e. the data returned by ``getPhi()``) is always synchronized.

Deposition of other quantities
______________________________

One can also deposit the following quantities on the mesh:

* Conductivity, which deposits :math:`\mu W`.
* Diffusivity, which deposits :math:`D W`.

Here, :math:`W` is the particle weight, :math:`\mu` is the particle mobility, :math:`D` is the particle diffusion coefficient.
It is up to the user to first interpolate or directly set the particle mobilities and diffusion coefficients before depositing the conductivity onto the mesh.

Functionality for the above deposited quantities exist as the following functions:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 219-231,243-255
   :dedent: 2

.. _Chap:ItoInterpolation:

Particle interpolation
----------------------

Interpolating particle velocities for ``ItoSolver`` is done by interpolating the mobility and particle velocities to the mesh,

.. math::

   \mathbf{V} = \mu\left(\mathbf{X}\right) \mathbf{v}\left(\mathbf{X}\right).

There is, however, some freedom in choosing how the mobility coefficient is calculated, which is discussed below.
In either case, there is some interpolation from a mesh-based variable onto the particle position :math:`\mathbf{X}`.
This interpolation method is always parsed from an options file, and is usually an NGP or CIC scheme.

.. important::

   When interpolating particle properties from the mesh, the user must first ensure that ghost cells are properly updated.


The separation into a mobility function and a velocity field is motivated by the introduction of an electric conductivity that permits a rather simple velocity relation as :math:`\mathbf{v} = \mu\mathbf{E}`, where :math:`\mathbf{E}` is the electric field.
Complete interpolation of the particle velocity consists of calling two functions:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 809-814,791-797
   :dedent: 2

Here, the calling sequence is such that the mobilities must be interpolated first, and then the velocity fields. 

Mobility coefficient interpolation
__________________________________

The mobility coefficient of a particle is usually interpolated directly, i.e.,

.. math::

   \mu = \mu\left(\mathbf{X}\right).

The other option is to compute the mobility as

.. math::

   \mu = \frac{\left(\mu\left|\mathbf{v}\right|\right)\left(\mathbf{X}\right)}{\left|\mathbf{v}\left(\mathbf{X}\right)\right|}.

This method ensures that the particle velocity becomes :math:`\mathbf{V} = \left(\mu\mathbf{v}\right)\left(\mathbf{X}\right)`.

.. tip::

   One can switch between the two interpolation methods in the ``ItoSolver`` run-time input options.

Diffusion coefficient interpolation
___________________________________

Interpolation of the diffusion coefficient is always done using an interpolation method

.. math::

   D = D\left(\mathbf{X}\right).

The function signature is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 824-829
   :dedent: 2

Particle intersections
----------------------

It will happen that particles occasionally hit the embedded boundary or leave through the domain sides.
In this case one might want to keep the particles in separate data holders rather than discard them.
``ItoSolver`` supplies several functions for transferring the particles to separate data containers when they intersect the EB or domain.
The most relevant function is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 482-489
   :dedent: 2

Here, ``EBIntersection`` is just an enum for putting logic into how the intersection is computed.
Valid options are ``EBIntersection::Bisection`` and ``EBIntersection::Raycast``.
These algorithms are discussed in :ref:`Chap:ParticleEB`.
The flag ``a_deleteParticles`` specifies if the original particles should be deleted when populating the other particle containers (again, see :ref:`Chap:ParticleEB`).

After calling ``intersectParticles``, the particles that crossed the EB or domain walls are available through the ``getParticles`` routine, see :ref:`Chap:ItoSolver` and can then be parsed separately by user code. 

Computing time steps
--------------------

While ``ItoSolver`` has no fundamental requirement on the time steps that can be used, several functions are available for computing various types of drift and diffusion related time steps.

.. important::

   All time step calculations below are imposed on the particles and not on the mesh variables.

Advective time step
___________________

The drift time step routines are implemented such that one restricts the time step such that the fastest particle does not move more than a specified number of grid cells.
This routine is implemented as

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 1100-1112
   :dedent: 2

which returns a CFL-like condition

.. math::

   \Delta t = \frac{\Delta x}{\textrm{max}(\left|v_x\right|, \left|v_y\right|, \left|v_z\right|)}.

Diffusive time step
___________________

The signatures for the diffusion time step are similar to the ones for drift:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 1133-1145
   :dedent: 2

which returns a CFL-like condition

.. math::

   \Delta t = \frac{\Delta x^2}{2dD},

where :math:`d` is the spatial dimension and :math:`D` is the particle diffusion coefficient.

Advective-diffusive time step
_____________________________

A combination of the advection and diffusion time step routines also exists as

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 993-1010
   :dedent: 2

This time step limitation is inspired by fully explicit and non-split fluid models, and is calculated as

.. math::

   \Delta t = \frac{1}{\frac{\Delta x}{\left|v_x\right| +  \left|v_y\right| +  \left|v_z\right|} + \frac{\Delta x^2}{2dD}}.

Superparticle management
------------------------

It can occasionally be necessary to merge or split computational particles.
This occurs in, e.g., plasma simulations where chemical reactions lead to exponential growth of particles. 
``ItoSolver`` handles superparticles via a configurable merger functor selected at parse time through ``ItoSolver.merge_algorithm``; the user can also supply a custom functor through ``setParticleCellMerger``.
The entry point for splitting and merging is in all cases

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 901-907
   :dedent: 2

Calling this function will merge/split the particles.

.. important::

   Most merging algorithms are performed within each grid cell, and particles must therefore be sorted by their cell index (``organizeParticlesByCell``) before calling the merging routine.
   The exceptions are ``nn_pair_tree``, ``nn_pair_onecell``, ``nn_pair_hash``, ``kd_carve``, ``kd_patch``, and ``kd_skin_nn``, which are distributed AMR-level merges dispatched over the whole container rather than cell by cell.
   All of these except ``kd_patch`` match particles across patch and rank boundaries and therefore require that a particle ghost halo has been filled; ``kd_patch`` is patch-local and instead requires that no ghost halo is present.

In order to specify the merging algorithm the user must set the ``ItoSolver.merge_algorithm`` to one of the following:

* ``none`` - No particle merging/splitting is performed.
* ``equal_weight_kd`` Use a kd-tree with bounding volume hierarchies to partition and split/merge the particles. This conserves the particle center-of-mass.
* ``reinitialize`` Re-initialize the particles in each grid cell, ensuring that weights are as uniform as possible.
* ``reinitialize_bvh`` Re-initialize the particles in each node of a kd-tree. Weights are as uniform as possible.
* ``nn_sfc`` Reach the target particle count by space-filling-curve nearest-neighbour clustering: when there are more particles than the target the nearest neighbours (along a Hilbert curve) are merged until the target count is reached, and when there are fewer the highest-weight particles are split. This gives spatially tight groups but does not equalize the weights.
* ``nn_pair_tree`` A distributed, MPI-safe nearest-neighbour *pair* merge that reaches the target particle count over the whole AMR hierarchy, searching for candidates via one whole-patch PointCloudBVH per patch. Over-full cells are drained by matching each over-crowded particle with its true nearest neighbour across patch and rank boundaries (a propose/judge/verdict protocol over a particle ghost halo) and merging the pair to its weighted centroid; because a single round merges pairs, the round is repeated until every cell reaches the target, or until ``ItoSolver.nn_pair_max_rounds`` rounds have run. Under-full cells are then brought up to the target by splitting the heaviest particle into two co-located daughters (floor/ceil weights, so integer weights stay integer). Tunable through ``ItoSolver.nn_pair_iterate``, ``ItoSolver.nn_pair_fallback``, ``ItoSolver.nn_pair_max_rounds`` and ``ItoSolver.nn_pair_max_cell_dist``.
* ``nn_pair_onecell`` The same distributed nearest-neighbour pair merge and drain/split protocol as ``nn_pair_tree``, but candidates are found via one PointCloudBVH per occupied grid cell instead of one per patch: a query only ever searches its own cell and its Moore-adjacent neighbours, so the merge distance is structurally fixed at Chebyshev cell distance 1 and ``ItoSolver.nn_pair_max_cell_dist`` does not apply. Tunable through ``ItoSolver.nn_pair_iterate``, ``ItoSolver.nn_pair_fallback`` and ``ItoSolver.nn_pair_max_rounds``.
* ``nn_pair_hash`` The same distributed nearest-neighbour pair merge and drain/split protocol as ``nn_pair_tree``, but candidates are found via one PointCloudHashGrid (a uniform spatial hash grid) per patch instead of a PointCloudBVH. Identical tunables and behaviour to ``nn_pair_tree`` (``ItoSolver.nn_pair_iterate``, ``ItoSolver.nn_pair_fallback``, ``ItoSolver.nn_pair_max_rounds``, ``ItoSolver.nn_pair_max_cell_dist``); only the per-patch spatial-index backend differs.
* ``kd_carve`` A distributed, MPI-safe whole-patch merge built around a spatial partition ("kd tree") rather than a nearest-neighbour graph: each patch splits its local-plus-ghost particles purely by position, bisecting the longest axis and never snapping to the grid, so particles near a cell face can merge across it -- at the count median while a node is still larger than ``ItoSolver.kd_split_weight_leaf_dx`` cell widths, and at the weight median once it is smaller, so that the resulting super-particles come out with comparable weights. How far splitting goes is set by a live per-cell quota: every group becomes exactly one super-particle at its weighted centroid, so the number of groups centred in a cell is that cell's post-merge population, and any split that would push a cell past the target count is refused. A group entirely clear of any patch/rank boundary merges immediately with no communication; a group that touches one is resolved by a single, non-iterative "z-buffer carve" -- competing groups from neighbouring patches are ranked by a deterministic key and the tightest one wins each contested particle. Unlike the ``nn_pair_*`` family there is no drain loop and no ``nn_pair_iterate``/``nn_pair_fallback``/``nn_pair_max_cell_dist`` equivalent to tune; the ghost width is fixed at 1. The maximum per-axis extent of a mergeable group is fixed at one cell width -- a physical safety bound applied as a filter on the finished partition, not a tunable, and matching the hardcoded width-1 ghost halo. Tunable through ``ItoSolver.kd_split_weight_leaf_dx`` (count-median/weight-median crossover, in cell widths; 0 disables the weight median).
* ``kd_patch`` The same whole-patch kd-tree build and per-cell quota as ``kd_carve``, but with no boundary tier: no particle ghost halo is filled, every group merges regardless of boundary exposure, and no particle is ever contested, so each patch reduces only the particles it owns. Cheaper and communication-free, at the cost of no coordination across patch boundaries.
* ``kd_skin_nn`` The same whole-patch kd-tree build and per-cell quota as ``kd_carve``, but the boundary tier is the nearest-neighbour pair merge instead of the carve arbitration. The particles are split across two containers: every group holding no ghost particle is committed into one of them (holding no ghost is the whole safety condition -- such a group's members are all resident in this patch, and any neighbouring patch that draws one of them into a group of its own sees it as a ghost and is disqualified by the same test), and the contested remainder -- the *skin* -- is drained from the other by the same algorithm as ``nn_pair_onecell``, whose Chebyshev-1 search radius matches the width-1 ghost halo exactly. Unlike ``kd_carve`` boundary exposure is never consulted, so an exposed but uncontested group merges locally, which is what makes the skin small. The skin tier drains each cell to what is *left* of its target after the already-merged interior results are counted against it, so the two tiers share one budget. Tunable through ``ItoSolver.kd_split_weight_leaf_dx`` and, for the skin tier, ``ItoSolver.nn_pair_iterate``, ``ItoSolver.nn_pair_fallback`` and ``ItoSolver.nn_pair_max_rounds``; ``ItoSolver.nn_pair_max_cell_dist`` does not apply.
* ``external`` Use an externally injected particle merging algorithm. In order to use this feature the user must supply one through ``setParticleCellMerger``.

The user can set the merging algorithm through the input script (see :ref:`Chap:ItoInput`), or supply one externally by setting the merge algorithm to ``external``.
In addition, the user must first supply a particle merging function:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 93-103
   :dedent: 2

In the code above, ``ParticleManagement::ParticleMerger<P>`` is an alias:

.. literalinclude:: ../../../../Source/Particle/CD_ParticleManagement.H
   :language: c++
   :lines: 82-84

.. tip::
   
   ``ItoSolver`` uses the kd-tree implementation from :ref:`Chap:SuperParticles` and partitioners for splitting the particles into two subsets with equal weights.

Example transport kernel
------------------------

Transport kernels for the particles within ``ItoSolver`` will typically be imposed externally by the user through a ``TimeStepper`` subclass that advances the particles.
For completeness, we here include a simple transport kernel for the ``ItoSolver`` which simply consists of a drift-diffusion kick.
``ItoParticle`` is a Struct-of-Arrays payload (see :ref:`Chap:Particles`), so the kernel operates on a ``ParticleSoA<ItoParticle>`` leaf and addresses each particle by index; the position is a container-owned column accessed through ``position(i)``/``setPosition(i, ...)``, while the interpolated velocity, diffusion coefficient, and old position are payload columns accessed through ``get<...>(i)``.
The loop below is the per-patch inner kernel and is run inside the usual level/patch iteration (see :ref:`Chap:Particles`):

.. code-block:: c++

   // One grid patch. The velocity columns (vx/vy/vz) have already been filled
   // by ItoSolver::interpolateVelocities().
   ParticleSoA<ItoParticle>& leaf = particles[lvl][dit()];

   for (std::size_t i = 0; i < leaf.size(); i++) {
      const RealVect      x = leaf.position(i);
      const RealVect      v = RealVect(D_DECL(leaf.get<&ItoParticle::vx>(i),
                                              leaf.get<&ItoParticle::vy>(i),
                                              leaf.get<&ItoParticle::vz>(i)));
      const ParticleReal& D = leaf.get<&ItoParticle::diffusion>(i);

      // Store the old position in the payload's old-position columns.
      D_TERM(leaf.get<&ItoParticle::old_x>(i) = x[0];,
             leaf.get<&ItoParticle::old_y>(i) = x[1];,
             leaf.get<&ItoParticle::old_z>(i) = x[2];);

      // Drift-diffusion kick.
      leaf.setPosition(i, x + v * a_dt + sqrt(2.0 * D * a_dt) * this->randomGaussian());
   }

The function ``randomGaussian`` implements a diffusion hopping and returns a 2D/3D dimensional vector with values drawn from a normal distribution with standard width of one and mean value of zero.
The implementation uses the random number generators in :ref:`Chap:Random`.

.. _Chap:ItoIO:

I/O
---

.. _Chap:ItoPlot:

Plot files
__________



For a complete list of available plot variables, see :ref:`Chap:ItoInput`. 

..
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

Several input options are available for configuring the run-time configuration of ``ItoSolver``, which are listed in :numref:`ItoInputOptions`.

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.options
   :name: ItoInputOptions
   :caption: Input options for the ``ItoSolver`` class.
	     Most are run-time configurable; ``seed`` and ``diffusion_grad_drift`` are read once at setup.
	     

Plot file variables
___________________

Plot variables are specified using ``ItoSolver.plt_vars``, see :ref:`Chap:ItoPlot`.
To add a variable to HDF5 output files, one can modify the ``ItoSolver.plt_vars`` input variable to include, e.g., the following variables:

* :math:`\phi`, i.e. the deposited particle weights (``ItoSolver.plt_vars = phi``)
* :math:`\mathbf{v}`, the advection field (``ItoSolver.plt_vars = vel``).
* :math:`D`, the diffusion coefficient  (``ItoSolver.plt_vars = dco``).

Particle-mesh configuration
___________________________

To specify the mobility interpolation, use ``ItoSolver.mobility_interp``.
Valid options are ``direct`` and ``velocity``, see :ref:`Chap:ItoInterpolation`.

Deposition and coarse-fine deposition (see :ref:`Chap:ParticleMesh`) are controlled using the flags

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

Example application(s)
----------------------

Example applications that use ``ItoSolver`` are found in

* :file:`$DISCHARGE_HOME/Physics/BrownianWalker`, see :ref:`Chap:BrownianWalkerModel`.
* :file:`$DISCHARGE_HOME/Physics/ItoKMC`, see :ref:`Chap:ItoKMC`.  
