.. _Chap:ItoDiffusion:

Îto diffusion
=============

The Îto diffusion model advances computational particles as drifting Brownian walkers

.. math::
   :label: ito_diffusion
	   
   \Delta\mathbf{X} = \left(\mathbf{V} + \nabla D\right)\Delta t + \sqrt{2D\Delta t}\mathbf{W}

where :math:`\mathbf{X}` is the spatial position of a particle, :math:`\mathbf{V}` the particle drift velocity, and :math:`D` is the diffusion coefficient *in the continuum limit*.
The vector term :math:`\mathbf{W}` indicates a random number sampled from a Gaussian distribution with mean value of 0 and standard deviation of 1.

.. _Chap:ItoDiffusionGradient:

The :math:`\nabla D` drift correction
-------------------------------------

:eq:`ito_diffusion` is an Îto-interpretation update, and the :math:`\nabla D` term in its drift is what makes it
transport the same diffusive flux as the fluid solvers.
Without that term the update is :math:`\Delta\mathbf{X} = \mathbf{V}\Delta t + \sqrt{2D\Delta t}\mathbf{W}`,
whose Fokker-Planck equation is

.. math::

   \frac{\partial n}{\partial t} = -\nabla\cdot\left(\mathbf{V}n\right) + \nabla^2\left(Dn\right),

i.e. it transports the flux :math:`\mathbf{V}n - n\nabla D - D\nabla n`.
:ref:`Chap:CdrSolver` integrates

.. math::

   \frac{\partial n}{\partial t} = -\nabla\cdot\left(\mathbf{V}n\right) + \nabla\cdot\left(D\nabla n\right),

i.e. the flux :math:`\mathbf{V}n - D\nabla n`.
The two differ by :math:`\nabla\cdot\left(n\nabla D\right)` and coincide exactly where :math:`D` is constant in
space, which is why the discrepancy is invisible in any test with a constant diffusion coefficient.
Adding :math:`\nabla D` to the drift, as :eq:`ito_diffusion` does, recovers the fluid form; equivalently, it
integrates the stochastic equation in the Hänggi-Klimontovich (anti-Îto) interpretation.

The correction is on by default and can be turned off with

.. code-block:: text

   ItoSolver.diffusion_grad_drift = true

in which case the solver reverts to the plain Îto update and transports :math:`\mathbf{V}n - \nabla\left(Dn\right)`.
``ItoSolver`` computes :math:`\nabla D` from the mesh diffusion coefficient and interpolates it to the particle
positions; it is the time stepper that adds :math:`\Delta t\nabla D` to the particle displacement.
Setting ``ItoSolver.plt_vars`` to include ``grad_dco`` writes :math:`\nabla D` to the plot files.

.. note::

   The correction is defined for a scalar, isotropic :math:`D`.
   It is also only meaningful when :math:`D` is a mesh field: a solver that assigns particle diffusion
   coefficients parametrically from the particle energy has no mesh field whose gradient means anything.
   Finally, the time step estimates in :ref:`Chap:ItoTimeStep` do not account for this term -- see the
   ``computeDt`` documentation for why.

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
   :lines: 730-736
   :dedent: 2

Usually, ``ItoSolver`` will perform a drift-diffusion advance and the user will then check if some of the particles crossed into the EB.
The solver can then automatically fill the boundary particles containers, see :ref:`Chap:ParticleIntersection`.

Remapping particles
-------------------

``ItoSolver`` has two functions for remapping particles:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 995-1006
   :dedent: 2

The bottom function lets the user remap any ``ParticleContainer<ItoParticle>`` that lives in the solver.
Here, ``a_container`` indicates which particle container to remap.

Particle deposition
-------------------

``ItoSolver`` contains several member functions for depositing various particle properties onto the mesh.
The most general version is given below:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 396-413
   :dedent: 2

This version permits the user to deposit an arbitrary per-particle quantity from a particle container ``a_particles`` onto some pre-allocated mesh storage ``a_phi``.
The quantity to be deposited is supplied through the ``a_gather`` callable, which returns a ``Real`` value for each particle in the SoA container.

.. important::

   The ``ItoSolver`` deposition methods are specified in the input script, see :ref:`Chap:ItoInput`.
   Both the base deposition scheme (e.g., NGP or CIC) must be specified, as well as the handling near refinement boundaries. 

A simpler version that deposits the bulk particles as a density on the mesh is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 327-335
   :dedent: 2

The particles are deposited into the class member ``m_phi``, which stores the particle density on the mesh. 
This data can then be fetched with

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 753-758
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
   :lines: 363-372
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
   :lines: 229-241,253-265
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
   :lines: 848-853,830-836
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
   :lines: 863-868
   :dedent: 2

Particle intersections
----------------------

It will happen that particles occasionally hit the embedded boundary or leave through the domain sides.
In this case one might want to keep the particles in separate data holders rather than discard them.
``ItoSolver`` supplies several functions for transferring the particles to separate data containers when they intersect the EB or domain.
The most relevant function is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 521-528
   :dedent: 2

Here, ``EBIntersection`` is just an enum for putting logic into how the intersection is computed.
Valid options are ``EBIntersection::Bisection`` and ``EBIntersection::Raycast``.
These algorithms are discussed in :ref:`Chap:ParticleEB`.
The flag ``a_deleteParticles`` specifies if the original particles should be deleted when populating the other particle containers (again, see :ref:`Chap:ParticleEB`).

After calling ``intersectParticles``, the particles that crossed the EB or domain walls are available through the ``getParticles`` routine, see :ref:`Chap:ItoSolver` and can then be parsed separately by user code. 

.. _Chap:ItoTimeStep:

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
   :lines: 1139-1151
   :dedent: 2

which returns a CFL-like condition

.. math::

   \Delta t = \frac{\Delta x}{\textrm{max}(\left|v_x\right|, \left|v_y\right|, \left|v_z\right|)}.

Diffusive time step
___________________

The signatures for the diffusion time step are similar to the ones for drift:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 1172-1184
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
   :lines: 1032-1049
   :dedent: 2

This time step limitation is inspired by fully explicit and non-split fluid models, and is calculated as

.. math::

   \Delta t = \frac{1}{\frac{\Delta x}{\left|v_x\right| +  \left|v_y\right| +  \left|v_z\right|} + \frac{\Delta x^2}{2dD}}.

Superparticle management
------------------------

It can occasionally be necessary to merge or split computational particles.
This occurs in, e.g., plasma simulations where chemical reactions lead to exponential growth of particles. 
``ItoSolver`` handles superparticles via a configurable merger functor selected at parse time through ``ItoSolver.merge_method``; the user can also supply a custom functor through ``setParticleCellMerger``.
The entry point for splitting and merging is in all cases

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 940-946
   :dedent: 2

Calling this function will merge/split the particles.

.. important::

   Most merging algorithms are performed within each grid cell, and particles must therefore be sorted by their cell index (``organizeParticlesByCell``) before calling the merging routine.
   Both ``kd_cell`` and ``reinitialize`` are of this kind.
   The exceptions are ``kd_patch``, ``kd_amr`` and ``nn_amr``, which are distributed merges dispatched over the whole container rather than cell by cell.
   ``kd_amr`` and ``nn_amr`` match particles across patch and rank boundaries and therefore require that a particle ghost halo has been filled; ``kd_patch`` is patch-local and instead requires that no ghost halo is present.

The merging algorithm is selected in two steps: ``ItoSolver.merge_method`` names the *scope* the merge groups particles over, and a small set of ``kd_*``/``nn_*`` specifiers say how it groups within that scope.
Splitting the two apart is what keeps the selector list short: ``kd_patch`` and ``kd_amr`` run the same tree build over a different scope, and the three nearest-neighbour backends differ only in the spatial index they search.

``ItoSolver.merge_method`` takes one of:

* ``none`` - No particle merging/splitting is performed.
* ``kd_cell`` Partition each grid cell's particles with a kd tree and reduce every leaf to one particle. Cell-scoped, so a leaf never reaches past the cell it belongs to and the per-cell density is preserved exactly. Reads ``ItoSolver.kd_partition`` and ``ItoSolver.kd_placement``.
* ``kd_patch`` The same tree build, but one tree per *patch* rather than per cell: nodes are bisected on the longest axis and are never snapped to the grid, so a leaf can straddle a cell face and merge particles across it. How far splitting goes is set by a live per-cell quota -- the number of leaves centred in a cell is that cell's post-merge population, and a split that would push a cell past ``ItoSolver.particles_per_cell`` is refused. No particle ghost halo is filled and no particle is ever contested, so each patch reduces only the particles it owns: communication-free, at the cost of no coordination across patch boundaries.
* ``kd_amr`` The same whole-patch tree build and per-cell quota as ``kd_patch``, but resolved across patch and rank boundaries, so a leaf may draw members owned by a neighbour. ``ItoSolver.kd_amr_boundary`` picks how that is done:

  * ``carve`` -- a single, non-iterative "z-buffer carve": competing leaves from neighbouring patches are ranked by a deterministic key and the tightest one wins each contested particle. There is no drain loop and nothing to tune beyond the tree itself; the ghost width is fixed at 1, and the maximum per-axis extent of a mergeable leaf is fixed at one cell width to match it.
  * ``nn`` -- every leaf holding no ghost particle is committed immediately (holding no ghost is the whole safety condition: such a leaf's members are all resident in this patch, and any neighbouring patch that draws one of them into a leaf of its own sees it as a ghost and is disqualified by the same test), and the contested remainder -- the *skin* -- is drained by the nearest-neighbour pair merge, whose Chebyshev-1 search radius matches the width-1 ghost halo exactly. Boundary exposure is never consulted, so an exposed but uncontested leaf merges locally, which is what keeps the skin small. The two tiers share one per-cell budget: the skin drains each cell to what is *left* of its target after the interior results are counted against it. Also reads ``ItoSolver.nn_iterate``, ``ItoSolver.nn_fallback`` and ``ItoSolver.nn_max_rounds``; ``ItoSolver.nn_max_cell_dist`` does not apply.
* ``nn_amr`` A distributed, MPI-safe nearest-neighbour *pair* merge that reaches the target particle count over the whole AMR hierarchy. Over-full cells are drained by matching each over-crowded particle with its true nearest neighbour across patch and rank boundaries (a propose/judge/verdict protocol over a particle ghost halo) and merging the pair to its weighted centroid; because a single round merges pairs, the round is repeated until every cell reaches the target, or until ``ItoSolver.nn_max_rounds`` rounds have run. Under-full cells are then brought up to the target by splitting the heaviest particle into two co-located daughters (floor/ceil weights, so integer weights stay integer). ``ItoSolver.nn_search`` picks the spatial index the candidate search is backed by:

  * ``tree`` -- one whole-patch ``PointCloudBVH`` per patch.
  * ``hash`` -- one whole-patch ``PointCloudHashGrid`` (a uniform spatial hash grid) per patch. Identical behaviour to ``tree``; only the backend differs.
  * ``onecell`` -- one ``PointCloudBVH`` per occupied grid cell. A query only ever searches its own cell and its Moore-adjacent neighbours, so the merge distance is structurally fixed at Chebyshev cell distance 1 and ``ItoSolver.nn_max_cell_dist`` does not apply.

  Tunable through ``ItoSolver.nn_iterate``, ``ItoSolver.nn_fallback``, ``ItoSolver.nn_max_rounds`` and (for ``tree``/``hash``) ``ItoSolver.nn_max_cell_dist``.
* ``reinitialize`` Discard each cell's spatial information entirely and rebuild it: the cell's physical particles are divided into near-equal integer weights and placed at uniformly random positions in the cell. Cell-scoped. Requires integer weights.
* ``external`` Use an externally injected particle merging algorithm. In order to use this feature the user must supply one through ``setParticleCellMerger``.

Every ``kd_*`` method reads the same two specifiers, because every one of them builds a kd tree and reduces each leaf to one particle -- only the scope of the tree differs.

* ``ItoSolver.kd_partition`` -- how a node is divided into two children. ``weight`` divides it so the two halves carry as nearly equal a *weight* as possible, which keeps the merged weights equal; the per-cell build divides a particle across the split when it must, so it can create particles during the build, while the patch and AMR builds never do and so cannot reach exact equality. ``count`` divides on particle *count* and never divides a particle, so it creates none and can never overshoot the target mid-build, but it leaves the merged weights uneven -- measurably so under repeated merging, where the effective particle count per cell falls well below the target. ``hybrid`` is the two combined by node size: the count median while a node is wider than ``ItoSolver.kd_hybrid_leaf_dx`` cell widths, the weight median once it is narrower. That is the rule the patch and AMR scopes want, since a split of a node spanning several cells apportions leaves *between* cells, which is a question of counts, while a narrow node is where the weight median earns its keep. ``weight`` and ``count`` are the two limits of the same rule -- a crossover above every node's size and one below it -- which is why one setting serves all three scopes. ``weight_capped`` is ``weight`` with a pre-pass that divides every particle heavier than a leaf's target weight into pieces before the tree is built; a particle heavier than that target cannot be equalized away by any partition, since it alone sets a floor on whichever leaf it lands in. The pre-pass matters most for ``kd_patch``/``kd_amr``, whose build may not divide a particle mid-tree at all and so cannot otherwise reach equal weights.
* ``ItoSolver.split_placement`` -- where the pieces of a *split* particle go. A split divides one particle's weight among several pieces; ``center`` puts every piece at the parent's position, ``jitter`` draws each from the local mean interparticle spacing around the parent, and ``cell`` draws each uniformly in the owning cell. All three keep every piece inside the parent's own cell, so the merge still moves no mass across a cell face, and all three fall back to ``center`` in a cut cell.

  ``center`` perturbs nothing spatially, but it makes the pieces indistinguishable, and a placement that inherits member positions (``sample``) then reproduces them. Over repeated merges the number of *distinct* positions a cell holds falls well below its particle count -- measured at a median of 6 distinct positions per 16 particles after 500 merges. That loss is invisible in the weight statistics, which stay perfectly uniform, while the cell's effective spatial sampling degrades by the same factor. ``jitter`` restores it, at the cost of moving mass by less than the distance between neighbouring particles.

.. note::

   The jitter kernel is *truncated* to the cell, not clamped and not reflected. Clamping stacks every out-of-range draw onto the wall. Reflecting folds the overhanging part of the kernel back onto the part that stayed, doubling the density within half a bandwidth of each wall -- and because the fold happens independently per direction, a parent near a corner is folded once per direction and the corner receives :math:`2^D` times the density it should. Truncation is the correct conditional distribution and has no directional structure.

* ``ItoSolver.kd_placement`` -- where the particle a leaf reduces to is placed. ``centroid`` puts it at the leaf's weighted centroid, which conserves the leaf's centre of mass exactly. ``sample`` puts it at one of the leaf's own particles, drawn with probability proportional to weight. ``random`` puts it at a uniformly random point in the leaf's bounding box (falling back to the centroid in cut cells, to keep the particle inside the embedded boundary).

  The centroid keeps a leaf's mean position and discards its spread, so it contracts the sub-cell distribution a little at every merge and, since the leaves of a locally uniform population are near enough equal-volume sub-boxes, drives the particles onto a regular sub-cell lattice that is the same in every grid cell and therefore coherent across the whole domain. ``sample`` is distribution-preserving -- a weight-proportional draw from a partition of a sample is itself a sample of the same distribution -- so neither the lattice nor the contraction occurs, at the cost of conserving the per-leaf centre of mass in expectation rather than exactly. ``random`` also avoids the lattice, but draws from the leaf's bounding box rather than the leaf itself, which is not quite the same distribution. Prefer ``sample`` wherever sub-cell position matters -- near a refinement boundary, or near an embedded boundary, where a centroid is not a position any particle occupied.

.. note::

   For ``kd_cell`` the leaf lies inside one cell, so every placement keeps the leaf's weight in the cell it came from and the per-cell density is preserved exactly.
   For ``kd_patch`` and ``kd_amr`` the leaves are deliberately not snapped to the grid, so a leaf straddling a cell face may place its merged particle on either side whichever placement is chosen -- that is what lets those scopes merge across a face in the first place, and it is why their per-cell particle count is a target rather than a guarantee.
   In all three, a placement that would land inside an embedded boundary falls back to the leaf's centroid.


The user can set the merging algorithm through the input script (see :ref:`Chap:ItoInput`), or supply one externally by setting the merge algorithm to ``external``.
In addition, the user must first supply a particle merging function:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 94-104
   :dedent: 2

In the code above, ``ParticleManagement::ParticleMerger<P>`` is an alias:

.. literalinclude:: ../../../../Source/Particle/CD_ParticleManagement.H
   :language: c++
   :lines: 195-197

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

How the cut cells are treated when depositing is selected with a single flag,

* ``ItoSolver.deposition_irreg``, with valid options

  * ``native`` -- deposit as-is, with no cut-cell treatment. A cut cell then holds :math:`\kappa n`.
  * ``ngp`` -- put the particle's entire cloud in its own cell when that cell is cut.
  * ``mirror`` -- even extension of the density about the embedded boundary. **Under this option a cut
    cell holds** :math:`n`, **not** :math:`\kappa n`, unlike every other option here; consumers of
    :math:`\phi` must opt into that change of meaning knowingly. Incompatible with
    ``ItoSolver.deposition = ngp``.
  * ``redistribute`` -- hybrid divergence, with the cut-cell mass excess redistributed to the neighbours.
  * ``redistribute_blended`` -- as ``redistribute``, but blended with the non-conservative divergence.

These are mutually exclusive answers to one question, which is why they are one selector rather than
independent flags: ``mirror`` combined with redistribution would apply a :math:`1/\kappa` correction to
a field that no longer needs one.

.. note::

   ``mirror`` is **experimental**. It is verified to add mass only where it should and never to remove
   any, but its acceptance suite is not finished, and the change of meaning above -- a cut cell holding
   :math:`n` rather than :math:`\kappa n` -- has not been audited against every consumer of
   :math:`\phi`. Use it for deliberate experiments, not production runs.

Cut-cell *interpolation* is separate and remains a boolean:

* ``ItoSolver.irr_ngp_interp`` for enforcing NGP interpolation. Valid options are ``true`` or ``false``.

Checkpoint-restart
__________________

Available input options for the ``ItoSolver`` are listed below:

Example application(s)
----------------------

Example applications that use ``ItoSolver`` are found in

* :file:`$DISCHARGE_HOME/Physics/BrownianWalker`, see :ref:`Chap:BrownianWalkerModel`.
* :file:`$DISCHARGE_HOME/Physics/ItoKMC`, see :ref:`Chap:ItoKMC`.  
