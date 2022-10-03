.. _Chap:RadiativeTransfer:

Radiative transfer
==================

Radiative transfer is supported in the diffusion (i.e. Eddington or Helmholtz) approximation and with Monte Carlo sampling of discrete photons.
The solvers share a common interface (a parent class), but note that the radiative transfer equation is inherently deterministic while Monte Carlo photon transport is inherently stochastic. 
The diffusion approximation relies on solving an elliptic equation in the stationary case and a parabolic equation in the time-dependent case, while the Monte-Carlo approach solves solves for fully transient or ''stationary'' transport.

.. note::

   The source code for the radiative transfer solvers reside in :file:`Source/RadiativeTransfer`

.. _Chap:RtSpecies:

RtSpecies
---------

The class ``RtSpecies`` is an abstract base class for parsing necessary information into radiative transfer solvers.
When creating a radiative transfer solver one will need to pass in a pointer to ``RtSpecies`` such that the solvers can look up the required infromation.
Currently, ``RtSpecies`` is a lightweight class where the user needs to implement the function

.. code-block:: c++

  virtual Real RtSpecies::getAbsorptionCoefficient(const RealVect a_pos) const = 0;

The absorption coefficient is used in the diffusion (see :ref:`Chap:DiffusionRTE`) and Monte Carlo (see :ref:`Chap:MonteCarloRTE`) solvers. 

One can also assign a name to the species through the member variable ``RtSpecies::m_name``.

.. _Chap:RtSolver:

RtSolver
--------

``RtSolver`` is the base class for encapsulating a radiative transfer solver.
The source code for the solver is located in :file:`$DISCHARGE_HOME/Source/RadiativeTransfer` and it is a fairly lightweight abstract class.
As with other solvers, ``RtSolver`` can use a specified :ref:`Chap:Realm`.

To use the ``RtSolver`` interface the user must cast from one of the inherited classes (see :ref:`Chap:DiffusionRTE` or :ref:`Chap:MonteCarloRTE`).
Since most of the ``RtSolver`` is an interface which is implemented by other radiative transfer solvers, documentation of boundary conditions, kernels and so on are found in the implementation classes.

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
Note that in the context below, :math:`\kappa` is *not* the volume fraction of a grid cell. 
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

For solving :eq:`TransientDiffusionRTE`, ``EddingtonSP1`` implements both the backward Euler method and the Twizell-Gumel-Arigu (TGA) scheme.
Explicit discretizations are not available. 
The Euler discretization is

.. math::

   \left(1+ \kappa \Delta t\right)\Psi^{k+1} - \Delta t \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi^{k+1}\right) = \Psi^{k} + \frac{\Delta t\eta^{k+1}}{c},

Again, this is a Helmholtz equation for `\Psi^{k+1}` which is solved using geometric multigrid. 
Expressions for the TGA scheme are found in :cite:`Twizell1996`, but note that the TGA scheme requires a solution to two elliptic equations (thus it has approximately twice the cost). 

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

.. code-block:: text
	     
   # ====================================================================================================
   # EddingtonSP1 class options
   # ====================================================================================================
   EddingtonSP1.stationary          = true         # Stationary solver
   EddingtonSP1.reflectivity        = 0.           # Reflectivity
   EddingtonSP1.use_tga             = false        # Use TGA for integration
   EddingtonSP1.kappa_scale         = true         # Kappa scale source or not (depends on algorithm)
   EddingtonSP1.plt_vars            = phi src      # Plot variables. Available are 'phi' and 'src'
   
   EddingtonSP1.ebbc                = larsen 0.0   # Bc on embedded boundaries
   EddingtonSP1.bc.x.lo             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.x.hi             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.y.lo             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.y.hi             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.z.lo             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.z.hi             = larsen 0.0   # Bc on domain side. 'dirichlet', 'neuman', or 'larsen'
   EddingtonSP1.bc.z.hi             = larsen 0.0   # Boundary on domain. 'neumann' or 'larsen'
   
   EddingtonSP1.gmg_verbosity       = -1           # GMG verbosity
   EddingtonSP1.gmg_pre_smooth      = 8            # Number of relaxations in downsweep
   EddingtonSP1.gmg_post_smooth     = 8            # Number of relaxations in upsweep
   EddingtonSP1.gmg_bott_smooth     = 8            # NUmber of relaxations before dropping to bottom solver
   EddingtonSP1.gmg_min_iter        = 5            # Minimum number of iterations
   EddingtonSP1.gmg_max_iter        = 32           # Maximum number of iterations
   EddingtonSP1.gmg_exit_tol        = 1.E-6        # Residue tolerance
   EddingtonSP1.gmg_exit_hang       = 0.2          # Solver hang
   EddingtonSP1.gmg_min_cells       = 16           # Bottom drop
   EddingtonSP1.gmg_bottom_solver   = bicgstab     # Bottom solver type. Valid options are 'simple <number>' and 'bicgstab'
   EddingtonSP1.gmg_cycle           = vcycle       # Cycle type. Only 'vcycle' supported for now
   EddingtonSP1.gmg_ebbc_weight     = 2            # EBBC weight (only for Dirichlet)
   EddingtonSP1.gmg_ebbc_order      = 2            # EBBC order (only for Dirichlet)
   EddingtonSP1.gmg_smoother        = red_black    # Relaxation type. 'jacobi', 'red_black', or 'multi_color'

Basic options
^^^^^^^^^^^^^

Basic input options to ``EddingtonSP1`` are as follows:

* ``EddingtonSP1.stationary`` for setting whether or not the solver is stationary.
* ``EddingtonSP1.reflectivity`` for controlling the reflectivity in the Larsen boundary conditions.
* ``EddingtonSP1.use_tga`` for switching between backward Euler and TGA time discretizations.
  Only relevant if ``EddingtonSP1.stationary = false``.
* ``EddingtonSP1.kappa_scale`` Switch for multiplying the source with with the volume fraction or not.
  Note that the multigrid Helmholtz solvers require a diagonal weighting of the operator.
  If ``EddingtonSP1.kappa_scale = false`` then the solver will assume that this weighting of the source term has already been made.
* ``EddingtonSP1.plt_vars`` For setting which solver plot variables are included in plot files. 

Setting boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Boundary conditions are parsed through the flags

* ``EddingtonSP1.ebbc`` Which sets the boundary conditions on the EBs.
* ``EddingtonSP1.bc.dim.side`` Which sets the boundary conditions on the domain sides, see :ref:`Chap:EddingtonSP1BC` for details. 

Tuning multigrid performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
  Sets the least squares stencil weighting factor for least squares gradient reconstruction on EBs when using Dirichlet boundary conditions. 
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

Monte Carlo methods
-------------------

All types of moment-closed radiative transfer equations contain nonphysical artifacts (which may or may not be acceptable).
For example, in the diffusion approximation the radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`, implying that photons can leak around boundaries.
I.e. the diffusion approximation does not correctly describe shadows.
It is possible to go beyond the diffusion approximation by also solving for higher-order moments like the radiative flux.
While such methods can describe shadows, they do, contain other nonphysical features.

Both ''stationary'' and transient Monte Carlo methods are offered as an alternative to the diffusion approximation. 

photon particle
_______________

The Îto particle is a computational particle class in `chombo-discharge` which can be used together with the particle tools in `Chombo`.
The following data fields are implemented in the particle:

.. code-block:: c++
   
   RealVect m_position;
   RealVect m_velocity;
   Real m_mass;
   Real m_kappa;

To obtain the fields, the user will call

.. code-block:: c++

   RealVect& position();
   RealVect& velocity();
   Real& mass();
   Real& diffusion();


All functions also have ``const`` versions.
Note that the field ``m_mass`` is the same as the *weight* of the computational particle.
The following functions are used to set the various properties:

.. code-block:: c++

   setPosition(const RealVect a_pos);
   setVelocity(const RealVect a_vel);
   setMass(const Real a_mass);
   setDiffusion(const Real a_diffusion;

Interaction with boundaries
___________________________


Stationary Monte Carlo
______________________

The stationary Monte Carlo method proceeds as follows.

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` where :math:`\mathcal{P}` is a Poisson distribution. The user may also choose to use pseudophotons rather than physical photons by modifying photon weights. Each photon is generated in the cell centroid :math:`\mathbf{x}_0` and given a random propagation direction :math:`\mathbf{n}`.

2. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`. The absorbed position of the photon is :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}` or reflect it off symmetry boundaries. 

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).
      
The Monte Carlo methods use computational particles for advancing the photons in exactly the same way a Particle-In-Cell method would use them for advancing electrons. Although a computational photon would normally live on the finest grid level that overlaps its position, this is not practical for all particle deposition kernels. For example, for cloud-in-cell deposition schemes it is useful to have the restrict the interpolation kernels to the grid level where the particle lives. In Chombo-speak, we therefore use a buffer region that extends some cells from a refinement boundary where the photons are not allowed to live. Instead, photons in that buffer region are transferred to a coarser level, and their deposition clouds are first interpolated to the fine level before deposition on the fine level happens. Selecting a deposition scheme and adjusting the buffer region is done through an input script associated with the solver. 
   
Transient Monte Carlo
_____________________

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time. The transient method is implemented as follows:

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` as above, and append these to the already existing photons. Each photon is given a uniformly distributed random creation time within :math:`\Delta t`. 
   
2. Each photon is advanced over the time step :math:`\Delta t` by a sequence of :math:`N` substeps (:math:`N` may be different for each photon).

   a. We compute :math:`N` such that we sample :math:`N\Delta \tau = \Delta t` with :math:`c\kappa\Delta\tau < 1`.

   b. A photon at position :math:`\mathbf{x}_0` is moved a distance :math:`\Delta \mathbf{x} = c\mathbf{n}\Delta\tau`. For each step we compute the absorption probability :math:`p = \kappa\left|\Delta\mathbf{x}\right|` where :math:`p\in[0,1]` is a uniform random number. If the photon is absorbed on this interval, draw a new uniform random number :math:`r \in [0,1]` and absorb the photon at the position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`. If the photon is not absorbed, it is moved to position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}`.

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).

Limitations
-----------

Example application
-------------------

An example application of usage of the ``RtSolver`` is found in :ref:`Chap:RadiativeTransferModel`. 
