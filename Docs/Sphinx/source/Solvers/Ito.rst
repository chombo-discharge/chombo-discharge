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
   :lines: 608-613
   :dedent: 2

Usually, ``ItoSolver`` will perform a drift-diffusion advance and the user will then check if some of the particles crossed into the EB.
The solver can then automatically fill the boundary particles containers, see :ref:`Chap:ParticleIntersection`.

Remapping particles
-------------------

``ItoSolver`` has two functions for remapping particles:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 848-859
   :dedent: 2

The bottom function lets the user remap any ``ParticleContainer<ItoParticle>`` that lives in the solver.
Here, ``a_container`` indicates which particle container to remap.

Particle deposition
-------------------

``ItoSolver`` contains several member functions for depositing various particle properties onto the mesh.
The most general version is given below:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 308-322
   :dedent: 2

This version permits the user to select any particle container ``a_particles`` and deposit them onto some pre-allocated mesh storage ``a_phi``.
Note that the template type ``P`` does not need to be ``ItoParticle``, although this is the most common use case.

.. important::

   The ``ItoSolver`` deposition methods are specified in the input script, see :ref:`Chap:ItoInput`.
   Both the base deposition scheme (e.g., NGP or CIC) must be specified, as well as the handling near refinement boundaries. 

A simpler version that deposits the bulk particles as a density on the mesh is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 268-274
   :dedent: 2

The particles are deposited into the class member ``m_phi``, which stores the particle density on the mesh. 
This data can then be fetched with

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 629-634
   :dedent: 2
   
For the full list of available deposition functions, see the ``ItoSolver`` C++ API `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classItoSolver.html>`_.

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
   :lines: 170-177, 192-201
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


The separation into a mobility function and a velocity field is motivated by the introduction of an electric conductivity that permits a rather simple velocity velocity relation as :math:`\mathbf{v} = \mu\mathbf{E}`, where :math:`\mathbf{E}` is the electric field.
Complete interpolation of the particle velocity consists of calling two functions:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 724-729, 705-712
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

The function signatures is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 738-743
   :dedent: 2

Particle intersections
----------------------

It will happen that particles occasionally hit the embedded boundary or leave through the domain sides.
In this case one might want to keep the particles in separate data holders rather than discard them.
``ItoSolver`` supplies several functions for transferring the particles to separate data containers when they intersect the EB or domain.
The most relevant function is

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 389-404
   :dedent: 2

Here, ``EbIntersection`` is a just an enum for putting logic into how the intersection is computed.
Valid options are ``EbIntersection::Bisection`` and ``EbIntersection::Raycast``.
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
   :lines: 962-966
   :dedent: 2

which returns a CFL-like condition

.. math::

   \Delta t = \frac{\Delta x}{\textrm{max}(\left|v_x\right|, \left|v_y\right|, \left|v_z\right|)}.

Diffusive time step
___________________

The signatures for the diffusion time step are similar to the ones for drift:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 983-987
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
   :lines: 881-890
   :dedent: 2

This time step limitation is inspired by fully explicit and non-split fluid models, and is calculated as

.. math::

   \Delta t = \frac{1}{\frac{\Delta x}{\left|v_x\right| +  \left|v_y\right| +  \left|v_z\right|} + \frac{\Delta x^2}{2dD}}.

Superparticle management
------------------------

It can occasionally be necessary to merge or split computational particles.
This occurs in, e.g., plasma simulations where chemical reactions lead to exponential growth of particles. 
``ItoSolver`` can currently handle superparticles through several internal functions, and is also equipped with an interface in which the user can inject an external particle-handling routine.  
The function for splitting and merging the particles is in all cases

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 766-772
   :dedent: 2

Calling this function will merge/split the particles.

.. important::
   
   Particle merging is currently performed within each grid cell, and particles must therefore be sorted by their cell index before calling the merging routine.


	    
In order to specify the merging algorithm the user must set the ``ItoSolver.merge_algorithm`` to one of the following:

* ``none`` - No particle merging/splitting is performed.
* ``equal_weight_kd`` Use a kD-tree with bounding volume hierarchies to partition and split/merge the particles. This conserves the particle center-of-mass.
* ``reinitialize`` Re-initialize the particles in each grid cell, ensuring that weights are as uniform as possible.
* ``reinitialize_bvh`` Re-initialize the particles in each node of a kD tree. Weights are as uniform as possible. 
* ``external`` Use an externally injected particle merging algorithm. In order to use this feature the user must supply one through

The user can set the merging algorithm through the input script (see :ref:`Chap:ItoInput`), or supply one externally by setting the merge algorithm to ``external``.
In addition, the user must first supply a particle merging function:

.. literalinclude:: ../../../../Source/ItoDiffusion/CD_ItoSolver.H
   :language: c++
   :lines: 68-73
   :dedent: 2

In the code above, ``ParticleManagement::ParticleMerger<P>`` is an alias:

.. literalinclude:: ../../../../Source/Particle/CD_ParticleManagement.H
   :language: c++
   :lines: 33-41
   :dedent: 2  

.. tip::
   
   ``ItoSolver`` uses the kD-node implementation from :ref:`Chap:SuperParticles` and partitioners for splitting the particles into two subsets with equal weights.

Example transport kernel
------------------------

Transport kernels for the particles within ``ItoSolver`` will typically be imposed externally by the user through a ``TimeStepper`` subclass that advances the particles.
For completeness, we here include a simple transport kernel for the ``ItoSolver`` which simply consists of a drift-diffusion kick:

.. code-block:: c++

   List<ItoParticle> particles;
   
   for (ListIterator<ItoParticle>& lit(particles); lit.ok(); ++lit) {
      ItoParticle& p = lit();

      p.oldPosition() = p.position();
      p.position()   += p.velocity() * a_dt + sqrt(2.0*p.diffusion()*a_dt) * this->randomGaussian();
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
	     All options are run-time configurable.
	     

Plot file variables
___________________

Plot variables are specified using ``ItoSolver.plt_vars``, see :ref:`Chap:ItoPlot`).
To add a variable to HDF5 output files, one can modify the ``ItoSolver.plt_vars`` input variable to include, e.g., the following variables:

* :math:`\phi`, i.e. the deposited particle weights (``ItoSolver.plt_vars = phi``)
* :math:`\mathbf{v}`, the advection field (``ItoSolver.plt_vars = vel``).
* :math:`D`, the diffusion coefficient  (``ItoSolver.plt_vars = dco``).

Particle-mesh configuration
___________________________

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

Example application(s)
----------------------

Example applications that use ``ItoSolver`` are found in

* :file:`$DISCHARGE_HOME/Physics/BrownianWalker`, see :ref:`Chap:BrownianWalkerModel`.
* :file:`$DISCHARGE_HOME/Physics/ItoKMC`, see :ref:`Chap:ItoKMC`.  
