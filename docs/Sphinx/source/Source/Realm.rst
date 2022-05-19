.. _Chap:Realm:

Realm
=====

``Realm`` is a class for centralizing EBAMR-related grids and operators for a specific AMR grid. 
For example, a ``Realm`` consists of a set of grids (i.e. a ``Vector<DisjointBoxLayout>``) as well as *operators*, e.g. functionality for filling ghost cells or averaging down a solution from a fine level to a coarse level.
One may think of a ``Realm`` as a fully-fledged AMR hierarchy with associated multilevel operators, i.e. how one would usually do AMR.



Dual grid
---------

The reason why ``Realm`` exists at all is due to individual load balancing of algorithmic components. 
The terminology *dual grid* is used when more than one ``Realm`` is used in a simulation, and in this case the user/developer has chosen to solve the equations of motion over a different set of ``DisjointBoxLayout`` on each level.
This approach is very useful when using computational particles since users can quickly generate separate Eulerian sets of grids for fluids and particles, and the grids can then be load balanced separately.
Note that every ``Realm`` consists of the same boxes, i.e. the physical domain and computational grids are the same for all realms. 
The difference lies primarily in the assignment of MPI ranks to grids; i.e. the load-balancing and domain decomposition.

.. _Fig:DualMesh:
.. figure:: /_static/figures/DualMesh.png
   :width: 360px
   :align: center

   Sketch of dual grid approach.
   Each rectangle represents a grid patch. a) Load balancing with the number of grid cells. b) Load balancing with the number of particles.

:numref:`Fig:DualMesh` shows an example of a dual-grid approach.
In  this figure we have a set of grid patches on a particular grid level.
In the top panel the grid patches are load-balanced using the grid patch volume as a proxy for the computational load.
The numbers in each grid patch indicates the MPI rank ownership of the patches.
In the bottom panel we have introduced computational particles in some of the patches.
For particles, the computational load is better defined by the number of computational particles assigned to the patch, and so using the number of particles as a proxy for the load yields different rank ownership over the grid patches.

Realm registration
------------------

To register a ``Realm``, users will have ``TimeStepper`` allocate the desired number of realms in the pure routine ``registerRealms()``, as follows:

.. code-block:: c++

   void myTimeStepper::registerRealms(){
      m_amr->registerRealm(Realm::Primal);
      m_amr->registerRealm("particleRealm");
      m_amr->registerRealm("otherParticleRealm");
   }

Since at least one realm is required, ``Driver`` will *always* register the realm ``"Primal"``.

During regrid, all realms are initially load balanced with the grid patch volume as the load proxy.
However, users can change load balancing individually for each realm through the load balancing routines in :ref:`Chap:TimeStepper`.    


Operator registration
---------------------

Internally, an instantiation of ``Realm`` contains the grids and the geometric information (e.g. ``EBISLayout``), as well as any operators that the user has seen fit to *register*.
Various operators are available for e.g. gradient stencils, conservative coarsening, ghost cell interpolation, filling a patch with interpolation data, redistribution, and so on.
Since operators always incur overhead and not all applications require *all* operators, they must be *registered*. 
If a solver needs an operator for, say, piecewise linear ghost cell interpolation, the solver needs to *register* that operator through the ``AmrMesh`` public interface:

.. code-block:: c++

   m_amr->registerOperator(s_eb_pwl_interp,   m_realm, m_phase);

Once an operator has been registered, ``Realm`` will define those operators during initialization e.g. regrids.
Run-time aborts with error messages are issued if an AMR operator is used, but has not been registered.

More commonly, ``chombo-discharge`` solvers will contain a routine that registers the operators that the solver needs.
A valid ``TimeStepper`` implementation *must* register all required operators in the function ``registerOperators()``, which is mostly as simple as:

.. code-block:: c++
		
   FieldSolver myPoissonSolver;
   CdrSolver myCdrSolver;

   void myTimeStepper::registerOperators(){
      myPoissonSolver->registerOperators();
      myCdrSolver->registerOperators();
   }

This will register the operators needed for ``FieldSolver`` and ``CdrSolver`` solver classes.
Note that if the solvers register the same operators, these operators are only defined once in :ref:`Chap:AmrMesh`.

Available operators
-------------------

The current operators are currently available:

#. Gradient ``s_eb_gradient``.
#. Irregular cell centroid interpolation, ``s_eb_irreg_interp``.
#. Coarse grid conservative coarsening, ``s_eb_coar_ave``.
#. Piecewise linear interpolation (with slope limiters), ``s_eb_fill_patch``.
#. Linear ghost cell interpolation, ``s_eb_pwl_interp``.
#. Flux registers, ``s_eb_flux_reg``.
#. Redistribution registers, ``s_eb_redist``.
#. Non-conservative divergence stencils, ``s_eb_noncons_div``.
#. Multigrid interpolators, ``s_eb_multigrid`` (used for multigrid).     
#. Signed distance function defined on grid, ``s_levelset``.
#. Particle-mesh support, ``s_eb_particle_mesh``.   

Solvers will typically allocate a subset of these operators, but for multiphysics code that use both fluid and particles, most of these will be in use.


Interacting with realms
-----------------------

Users will not interact with ``Realm`` directly.
Every ``Realm`` is owned by ``AmrMesh``, and the user will only interact with realms through the public ``AmrMesh`` interface, for example by fetching operators for performing AMR operations. 
In addition, data that is defined on one realm can be copied to another; ``EBAMRData<T>`` takes care of this.
You will simply call a copier:

.. code-block:: c++

   EBAMRCellData realmOneData;
   EBAMRCellData realmTwoData;

   realmOneData.copy(realmTwoData);

The rest of the functionality uses the public interface of :ref:`Chap:AmrMesh`.
For example for coarsening of multifluid data:

.. code-block:: c++

   std::string multifluidRealm;
   MFAMRCellData multifluidData;
   AmrMesh amrMesh;

   amrMesh.averageDown(multifluidData, multifluidRealm);
