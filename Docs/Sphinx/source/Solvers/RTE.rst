.. _Chap:RadiativeTransfer:

Radiative transfer
==================

.. tip::

   The source code for the radiative transfer solvers reside in :file:`Source/RadiativeTransfer`

.. _Chap:RtSolver:   

RtSolver
--------

Radiative transfer solvers are supported in the form of

* Diffusion solvers, i.e. first order Eddington solvers, which takes the form of a Helmholtz equation.
* Using Monte Carlo sampling of discrete photons.
  
The solvers share a parent class ``RtSolver``, and code that uses only the ``RtSolver`` interface will should be able to switch between the two implementations.
Note, however, that the radiative transfer equation is inherently deterministic while Monte Carlo photon transport is inherently stochastic. 
The diffusion approximation relies on solving an elliptic equation in the stationary case and a parabolic equation in the time-dependent case, while the Monte-Carlo approach solves solves for fully transient or ''stationary'' transport.

.. tip::
   
   The source code for the solver is located in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer` and it is a fairly lightweight abstract class.
   As with other solvers, ``RtSolver`` can use a specified :ref:`Chap:Realm`.

To use the ``RtSolver`` interface the user must cast from one of the inherited classes (see :ref:`Chap:DiffusionRTE` or :ref:`Chap:MonteCarloRTE`).
Since most of the ``RtSolver`` is an interface which is implemented by other radiative transfer solvers, documentation of boundary conditions, kernels and so on are found in the implementation classes.

.. _Chap:RtSpecies:

RtSpecies
---------

The class ``RtSpecies`` is an abstract base class for parsing necessary information into radiative transfer solvers.
When creating a radiative transfer solver one will need to pass in a reference to an ``RtSpecies`` instantiation such that the solvers can look up the required infromation.
Currently, ``RtSpecies`` is a lightweight class where the user needs to implement the function

.. code-block:: c++

  virtual Real RtSpecies::getAbsorptionCoefficient(const RealVect a_pos) const = 0;

The absorption coefficient is used in the diffusion (see :ref:`Chap:DiffusionRTE`) and Monte Carlo (see :ref:`Chap:MonteCarloRTE`) solvers. 

One can also assign a name to the species through the member variable ``RtSpecies::m_name``.

.. _Chap:DiffusionRTE:

Diffusion approximation
-----------------------

EddingtonSP1
____________

The first-order diffusion approximation to the radiative transfer equation is encapsulated by the ``EddingtonSP1`` class which implements a first order Eddington approximation of the radiative transfer equation.
``EddingtonSP1`` implements ``RtSolver`` using both stationary and transient advance methods (e.g. for stationary or time-dependent radiative transport).
The source code is located in :file:`$DISCHARGE_HOME/RadiativeTransfer`. 

Equation of motion
__________________

In the diffusion approximation, the radiative transport equation is

.. math::
   :label: TransientDiffusionRTE

   \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

where :math:`\kappa` is the absorption coefficient (i.e., inverse absorption length).
Note that in the context below, :math:`\kappa` is *not* the volume fraction of a grid cell but the absorption coefficient.
This is called the Eddington approximation, and the radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`.

In the stationary case this yields a Helmholtz equation

.. math::
   :label: StationaryDiffusionRTE

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

Implementation
______________

``EddingtonSP1`` uses multigrid methods for solving :eq:`TransientDiffusionRTE` and :eq:`StationaryDiffusionRTE`, see :ref:`Chap:LinearSolvers`.
The class implements ``RtSolver::advance()``, which can switch between :eq:`TransientDiffusionRTE` and :eq:`StationaryDiffusionRTE`.
Note that for both the stationary and time-dependent cases the absorption coefficient :math:`\kappa` in :eq:`TransientDiffusionRTE` and :eq:`StationaryDiffusionRTE` are filled using the ``RtSpecies`` implementation provided to the solver.
Also note that the absorption coefficient does not need to be constant in space. 


Stationary kernel
^^^^^^^^^^^^^^^^^

For the stationary kernel we solve :eq:`StationaryDiffusionRTE` directly, using a single multigrid solve.
See :ref:`Chap:LinearSolvers` for discretization details. 

Transient kernel
^^^^^^^^^^^^^^^^

For solving :eq:`TransientDiffusionRTE`, ``EddingtonSP1`` implements the backward Euler method, while explicit discretizations are not currently available. 
The Euler discretization is

.. math::

   \left(1+ \kappa \Delta t\right)\Psi^{k+1} - \Delta t \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi^{k+1}\right) = \Psi^{k} + \frac{\Delta t\eta^{k+1}}{c},

Again, this is a Helmholtz equation for :math:`\Psi^{k+1}` which is solved using geometric multigrid.

.. _Chap:EddingtonSP1BC:
   
Boundary conditions
___________________

Simplified domain boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``EddingtonSP1`` solver supports the following boundary conditions on domain faces and EBs.
The domain boundary condition *type*, which is either Dirichlet, Neumann, or Larsen (a special type of Robin boundary condition) is always passed in through the input file.
If the user passes in a value, say ``neumann 0.0``, for a particular domain side/face, then the class will use a homogeneous Neumann boundary for the entire domain edge/face. 

Custom domain boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to use more complex boundary conditions by passing in ``dirichlet_custom``, ``neumann_custom``, or ``larsen_custom`` options.
In this case the ``EddingtonSP1`` solver will use a specified function at the domain edge/face.
To specify that function, ``EddingtonSP1`` has a member function

.. code-block:: c++

   void setDomainSideBcFunction(const int a_dir,
                                const Side::LoHiSide a_side,
				const std::function<Real(const RealVect a_pos, const Real a_time)> a_function);

which species a boundary condition value for one of the edges (faces in 3D).
Note that the boundary condition *type* is still Dirichlet, Neumann, or Larsen (depending on whether or not ``dirichlet_custom``, ``neumann_custom``, or ``larsen_custom`` was passed in). 
For example, to set the boundary condition on the left :math:`x` face in the domain, one can create a ``EddingtonSP1DomainBc::BcFunction`` object as follows:

.. code-block:: c++

   // Assume this has been instantiated. 
   RefCountedPtr<EddingtonSP1> eddingtonSolver;

   // Make a lambda which we can bind to std::function. 
   auto myValue = [](const RealVect a_pos, const Real a_time) -> Real {
      return a_pos[0] * a_time;
   }

   // Set the domain bc function in the solver. 
   eddingtonSolver.setDomainSideBcFunction(0, Side::Lo, myValue);

.. note::

   If the user specifies one of the custom boundary conditions but does not set the function, it will issue a run-time error.
   
Embedded boundaries
^^^^^^^^^^^^^^^^^^^

On the EB, we currently only support constant-value boundary conditions.
In the input script, the user can specify

* ``dirichlet <value>`` For setting a constant Dirichlet boundary condition everywhere. 
* ``neumann <value>`` For setting a constant Neumann boundary condition everywhere. 
* ``larsen <value>`` For setting a constant Larsen boundary condition everywhere. 

Boundary condition types
^^^^^^^^^^^^^^^^^^^^^^^^

#. **Dirichlet**.
   For Dirichlet boundary conditions we specify the value of :math:`\Psi` on the boundary.
   Note that this involves reconstructing the gradient :math:`\partial_n\Psi` on domain faces and edges, see :ref:`Chap:LinearSolverDirichletBC`. 
#. **Neumann**.
   For Neumann boundary conditions we specify the value of :math:`\partial_n\Psi` on the boundary.
   Note that the linear solver interface also supports setting :math:`B\partial_n\Psi` on the boundary (where :math:`B` is the Helmholtz equation :math:`B` coefficient).
   However, the ``EddingtonSP1`` solver does not use this functionality.
#. **Larsen**.
   The Larsen boundary condition is an absorbing boundary condition, taking the form of a Robin boundary as follows:

   .. math::

      \kappa\partial_n\Psi + \frac{3\kappa^2}{2}\frac{1-3r_2}{1-2r_1}\Psi = g,

   where :math:`r_1` and :math:`r_2` are reflection coefficients and :math:`g` is a surface source, see :cite:`Larsen2002` for details.
   Note that when the user specifies the boundary condition value (e.g. by setting the BC function), he is setting the surface sourge :math:`g`.
   In the majority of cases, however, we will have :math:`r_1 = r_2 = g = 0` and the BC becomes

   .. math::

      \partial_n\Psi + \frac{3\kappa}{2}\Psi = 0.

Solver configuration
____________________

The ``EddingtonSP1`` implementation has a number of configurable options for running the solver, and these are given below:

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_EddingtonSP1.options

Basic options
^^^^^^^^^^^^^

Basic input options to ``EddingtonSP1`` are as follows:

* ``EddingtonSP1.verbosity`` for controlling solver verbosity.
* ``EddingtonSP1.stationary`` for setting whether or not the solver is stationary.
* ``EddingtonSP1.reflectivity`` for controlling the reflectivity in the Larsen boundary conditions.
  Only relevant if ``EddingtonSP1.stationary = false``.
* ``EddingtonSP1.kappa_scale`` Switch for multiplying the source with with the volume fraction or not.
  Note that the multigrid Helmholtz solvers require a diagonal weighting of the operator, including the right-hand side.
  If ``EddingtonSP1.kappa_scale = false`` then the solver will assume that this weighting of the source term has already been made.
* ``EddingtonSP1.plt_vars`` For setting which solver plot variables are included in plot files.
* ``EddingtonSP1.use_regrid_slopes`` For setting turning on/off slopes when regridding the solution.

Setting boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Boundary conditions are parsed through the flags

* ``EddingtonSP1.ebbc`` Which sets the boundary conditions on the EBs.
* ``EddingtonSP1.bc.dim.side`` Which sets the boundary conditions on the domain sides, see :ref:`Chap:EddingtonSP1BC` for details. 

Multigrid settings
^^^^^^^^^^^^^^^^^^

All parameters that begin with the form ``EddingtonSP1.gmg_`` indicate a tuning parameter for geometric multigrid.

* ``EddingtonSP1.gmg_verbosity``.
  Controls the multigrid verbosity.
  Setting it to a number :math:`> 0` will print multigrid convergence information.
* ``EddingtonSP1.gmg_pre_smooth``.
  Controls the number of relaxations on each level during multigrid downsweeps.
* ``EddingtonSP1.gmg_post_smooth``.
  Controls the number of relaxations on each level during multigrid upsweeps.
* ``EddingtonSP1.gmg_bott_smooth``.
  Controls the number of relaxations before entering the bottom solve. 
* ``EddingtonSP1.gmg_min_iter``.
  Sets the minimum number of iterations that multigrid will perform. 
* ``EddingtonSP1.gmg_max_iter``.
  Sets the maximum number of iterations that multigrid will perform. 
* ``EddingtonSP1.gmg_exit_tol``.
  Sets the exit tolerance for multigrid.
  Multigrid will exit the iterations if :math:`r < \lambda r_0` where :math:`\lambda` is the specified tolerance, :math:`r = |L\Phi -\rho|` is the residual and :math:`r_0` is the residual for :math:`\Phi = 0`.  
* ``EddingtonSP1.gmg_exit_hang``.
  Sets the minimum permitted reduction in the convergence rate before exiting multigrid.
  Letting :math:`r^k` be the residual after :math:`k` multigrid cycles, multigrid will abort if the residual between levels is not reduce by at least a factor of :math:`r^{k+1} < (1-h)r^k`, where :math:`h` is the "hang" factor.
* ``EddingtonSP1.gmg_min_cells``.
  Sets the minimum amount of cells along any coordinate direction for coarsened levels.
  Note that this will control how far multigrid will coarsen. Setting a number ``gmg_min_cells = 16`` will terminate multigrid coarsening when the domain has 16 cells in any of the coordinate direction. 
* ``EddingtonSP1.gmg_bottom_solver``.
  Sets the bottom solver type. 
* ``EddingtonSP1.gmg_cycle``.
  Sets the multigrid method.
  Currently, only V-cycles are supported.
* ``EddingtonSP1.gmg_ebbc_order``.
  Sets the stencil order on EBs when using Dirichlet boundary conditions. 
  Note that this is also the stencil radius.
  See :ref:`Chap:LinearSolvers` for details. 
* ``EddingtonSP1.gmg_ebbc_weight``.
  Sets the least squares stencil weighting factor for least squares gradient reconstruction on EBs when using Robin or Dirichlet boundary conditions.
  See :ref:`Chap:LeastSquares` for details.   
* ``EddingtonSP1.gmg_smoother``.
  Sets the multigrid smoother.

Runtime parameters
^^^^^^^^^^^^^^^^^^

The following parameters for ``EddingtonSP1`` are run-time configurable:

* All multigrid tuning parameters, i.e. parameters starting with ``EddingtonSP1.gmg_``.
* Plot variables, i.e. ``EddingtonSP1.plt_vars``.
* Kappa scaling (for algorithmic adjustments), i.e. ``EddingtonSP1.kappa_scale``. 

.. _Chap:MonteCarloRTE:

Monte Carlo sampling
--------------------

``McPhoto`` defines a class which can solve radiative transfer problems using discrete photons that travel "instantenously" or transiently.
The class derives from ``RtSolver`` and can thus be used by problems that only require the ``RtSolver`` interface.
``McPhoto`` can provide a rather complex interaction with boundaries, such as computing the intersection between a photon path and a geometry, and thus it can capture e.g. shadows.

The Monte Carlo sampling is a particle-based radiative transfer solver, and particle-mesh operations (see :ref:`Chap:ParticleMesh`) are required in order to deposit the photons on a mesh when computing densities.

.. tip::

   The ``McPhoto`` class is defined in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer/CD_McPhoto.H`.
   See the `McPhoto C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classMcPhoto.html>`_ for further details.

The solver has multiple data holders for systemizing photons, which is especially useful during transport kernels where some of the photons might strike a boundary:

* In-flight photons
* Bulk-absorbed photons, i.e. photons absorbed on the mesh.
* EB-absorbed photons, i.e. photons that struck the EB during a transport step.
* Domain-absorbed photons, i.e. photons that struck the domain edge/face during a transport step.
* Source photons, for letting the user pass in externally generated photons into the solver.

Photon particle
_______________

The ``Photon`` particle is a simple encapsulation of a computational particle and is used by ``McPhoto``.
It derives from ``GenericParticle<2,1>`` and stores (in addition to the particle position):

* The particle weight.
* The particle mean absorption coefficient.
* The particle velocity/direction.

.. tip::

   The ``Photon`` class is defined in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer/CD_Photon.H`

When defining the ``McPhoto`` class, the particle's absorption coefficient is computed from the implementation of the absorption function method in :ref:`Chap:RtSpecies`.

Generating photons
__________________

There are several ways users can generate computational photons that are to be transported by the solver.

#. Fetch the *source photons* by calling ``McPhoto::getSourcePhotons()`` and fill the returned data holder.
   The photons can then be added to the ``McPhoto`` instantiation and one of the transport kernels can be called.

#. If the source term :math:`\eta` has been filled, the user can call ``McPhoto::advance`` to have the solver generate the computational photons and than transport them.

   .. important::

      The ``advance`` function is *only* meant to be used together with a mesh-based source term that the user has filled prior to calling the method.

      When using the ``advance``, the number of photons that are generated are limit to a user-specified number (see :ref:`Chap:McPhotoOptions` for further details).

Transport modes
_______________

``McPhoto`` can be run as a fully transient (in which photons are tracked in time) or as an instantaneous solver (where photons are absorbed immediately on the mesh).
These two differ in the way the transport problem over a time step :math:`\Delta t` is approach, but both methods include intersection tests with geometries and domain edges/faces, 

Instantaneous transport
^^^^^^^^^^^^^^^^^^^^^^^

When using instantaneous transport, any photon generated in a time step is immediately absorbed on the boundary through the following steps:

#. Optionally, have the solver generate photons to be transport (or add them externally).
#. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`.
   Here, :math:`\kappa` is computed by calling the underlying :ref:`Chap:RtSpecies` absorption function.
   The absorbed position of the photon is set to :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

   .. warning::

      In instantaneous mode photons might travel infinitely long, i.e. there is no guarantee that :math:`c\Delta t \leq r`.
#. Deposit the photons on the mesh.

Transient transport
^^^^^^^^^^^^^^^^^^^

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time.
For each photon, do the following:

#. Compute an absorption length :math:`r` by sampling the absorption function at the current photon position.
   
#. Each photon is advanced over the time step :math:`\Delta t` such that the position is

   .. math::

      \mathbf{x} = \mathbf{x}_0 + \mathbf{c}\Delta t.

#. Check if :math:`\left|\mathbf{x}-\mathbf{x}_0\right| < r` and if it is, absorb the photon on the mesh.
   
Other transport kernels
^^^^^^^^^^^^^^^^^^^^^^^

In addition to the above two methods, the solver interface permits users to add e.g. source photons externally and add them to the solvers' transport kernel. 

.. _Chap:McPhotoOptions:

Solver configuration
___________________

.. literalinclude:: ../../../../Source/RadiativeTransfer/CD_McPhoto.options

* ``McPhoto.verbosity`` for controlling the solver verbosity.
* ``McPhoto.instantaneous`` for setting the transport mode.
* ``McPhoto.max_photons_per_cell`` for restricting the number of photons generated per cell when having the solver generate the computational photons. This is only relevant when calling the ``advance`` method.
* ``McPhoto.num_sampling_packets`` for using sub-sampling when generating and transport photons in instantaneous mode through the ``advance`` function.
  This permits the ``McPhoto.max_photons_per_cell`` to partition the photon transport into packets where a fewer number of photons are generated during each step. Note that this will deposit the photons on the mesh for each packet, and the absorbed photons are only available as a density (i.e., the computational photons that were absorbed are lost).
  This can reduce memory for certain types of applications when using many computational photons.
* ``McPhoto.blend_conservation`` is a dead option marked for future removal (it blends a non-conservative divergence when depositing in cut-cells).
* ``McPhoto.transparent_eb`` for turning on/off transparent boundaries. Mostly used for debugging.
* ``McPhoto.plt_vars`` for setting plot variables. 
* ``McPhoto.intersection_alg`` sets the intersection algorithm when computing collisions with EBs.
  Ray-casting and bisection methods are supported.
* ``McPhoto.bisect_step`` sets bisection step (physical length) when calculation intersection tests using the bisection algorithm (i.e., this parameter is irrelevant if ``McPhoto.intersection_alg = raycast``).
* ``McPhoto.deposition`` for setting the deposition method.
  Currently, NGP and CIC methods are supported (see :ref:`Chap:ParticleMesh`).
* ``McPhoto.deposition_cf`` for setting the deposition strategy near coarse-fine boundaries.
  Currently, *interp* and *halo* are supported, see :ref:`Chap:ParticleMesh`.
* ``McPhoto_bc_<coord>_<low/high>`` sets the boundary condition on domain edges/faces.
* ``McPhoto.photon_generation`` for setting the photon generation method (details are given below).
* ``McPhoto.source_type`` for setting the photon generation method (details are given below).

Clarifications
^^^^^^^^^^^^^^

When computational photons are generated through the solver, users might have filled the source term differently depending on the application.
For example, users might have filled the source term with the number of photons generated per unit volume and time, or the *physical* number of photons to be generated. 
The two input options ``McPhoto.photon_generation`` and ``McPhoto.source_type`` contain the necessary specifications for ensuring that the user-filled source term can be translated properly for ensuring that the correct number of physical photons are generated.
Firstly, ``McPhoto.source_type`` contains the specification of what the source term contains, e.g.

* ``number`` if the source term contains the physical number of photons.
* ``volume`` if the source terms contains the physical number of photons generated per unit volume.
* ``volume_rate`` if the source terms contains the physical number of photons generated per unit volume and time.
* ``rate`` if the source terms contains the physical number of photons generated per unit time.

When ``McPhoto`` calculates the number of physical photons in a cell, it will automatically determine from ``McPhoto.source_type``, :math:`\Delta V` and :math:`\Delta t` how many physical photons are to be generated in each grid cell.

``McPhoto.photon_generation`` permits the user to turn on/off Poisson sampling when determining how many photons will be generated.
If this is set to *stochastic*, the solver will first compute the number of physical photons :math:`\overline{N}_\gamma^{\text{phys}}` following the procedure above, and then run a Poisson sampling such that the final number of physical photons is

.. math::

   N_{\gamma}^{\text{phys}} = P\left(\overline{N}_{\gamma}^{\text{phys}}\right).

Otherwise, if ``McPhoto.photon_generation`` is set to *deterministic* then the solver will generate

.. math::

   N_{\gamma}^{\text{phys}} = \overline{N}_{\gamma}^{\text{phys}}

photons.
Again, these elements are important because users might choose to run such stochastic samplings outside of ``McPhoto``.

.. important::

   All of the above procedures are done *per-cell*.
  

Example application
-------------------

Example applications that use ``RtSolver`` are found in:

* :ref:`Chap:RadiativeTransferModel`.
* :ref:`Chap:CdrPlasmaModel`.
* :ref:`Chap:ItoKMC`.    
