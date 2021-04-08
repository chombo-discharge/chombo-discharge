.. _Chap:PoissonSolver:
   
Poisson
=======

The code for the Poisson solver is given in :file:`/src/poisson` *and* :file:`/src/elliptic`.
The first folder contains the interface and code for only the Poisson solver, whereas the :file:`/src/elliptic` folder also contains generic code for elliptic problems.

The Poisson solver is named ``poisson_solver`` and currently only has one specific implementation: ``poisson_multifluid_gmg`` which uses a multifluid embedded boundary formulation together with geometric multigrid.

Setting up the solver
---------------------

In order to set up the solver, one must provide the

Creating a solver is often done with smart pointer casts like so:

.. code-block:: c++

   RefCountedPtr<poisson_solver> poisson = RefCountedPtr<poisson_solver> (new poisson_multifluid_gmg());

In addition, one must parse run-time options to the class, provide the ``amr_mesh`` and ``computational_geometry`` instances, and set the initial conditions.
This is done as follows:

.. code-block:: c++

   poisson->parse_options();              // Parse class options
   poisson->set_amr(amr);                 // Set amr - we assume that `amr` is an object
   poisson->set_computational_geometry(); // Set the computational geometry
   poisson->allocate_internals();         // Allocate storage for potential etc.
   poisson->set_potential(potential);     // Set the potential

The argument in the function ``set_potential(...)`` is a function pointer of the type

.. code-block:: c++

   Real potential(const Real a_time)

and allows setting a time-dependent potential for the solver. 

The `PlasmaC` field solver has a lot of supporting functionality, but essentially relies on only one critical function:
Solving for the potential.
This is done by calling a class-specific function

.. code-block:: c++

   bool solve(MFAMRCellData& phi, const MFAMRCellData& rho, const EBAMRIVData& sigma);

where ``phi`` is the resulting potential that was computing with the space charge density ``rho`` and surface charge density ``sigma``.

Currently, only one field solver is implemented and this solver uses a geometric multigrid method for solving for the potential.
The solver supports three phases: electrodes, gas, and dielectric.

Boundary conditions
-------------------

Domain boundary conditions for the solver must be set by the user through an input script, whereas the boundary conditions on internal surfaces are Dirichlet by default.
Note that on multifluid-boundaries the Dirichlet boundary condition is enforced by the conventional matching boundary condition that follows from Gauss` law.

Tuning multigrid performance
----------------------------

The Poisson equation is currently solved with a geometric multigrid method (GMG). 
Various switches are enabled that adjust the performance of GMG, and these are listed below

.. code-block:: bash

   # ====================================================================================================
   # POISSON_MULTIFLUID_GMG_GMG CLASS OPTIONS (MULTIFLUID GMG SOLVER SETTINGS)
   # ====================================================================================================
   poisson_multifluid_gmg.bc_x_low  = neumann           # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.bc_x_high = neumann           # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.bc_y_low  = dirichlet_ground  # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.bc_y_high = dirichlet_live    # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.bc_z_low  = neumann           # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.bc_z_high = neumann           # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
   poisson_multifluid_gmg.plt_vars  = phi rho E res     # Plot variables. Possible vars are 'phi', 'rho', 'E', 'res'

   poisson_multifluid_gmg.auto_tune         = false     # Do some auto-tuning
   poisson_multifluid_gmg.gmg_verbosity     = -1        # GMG verbosity
   poisson_multifluid_gmg.gmg_pre_smooth    = 12        # Number of relaxations in downsweep
   poisson_multifluid_gmg.gmg_post_smooth   = 12        # Number of relaxations in upsweep
   poisson_multifluid_gmg.gmg_bott_smooth   = 12        # NUmber of relaxations before dropping to bottom solver
   poisson_multifluid_gmg.gmg_min_iter      = 5         # Minimum number of iterations
   poisson_multifluid_gmg.gmg_max_iter      = 32        # Maximum number of iterations
   poisson_multifluid_gmg.gmg_tolerance     = 1.E-10    # Residue tolerance
   poisson_multifluid_gmg.gmg_hang          = 0.2       # Solver hang
   poisson_multifluid_gmg.gmg_bottom_drop   = 4         # Bottom drop
   poisson_multifluid_gmg.gmg_bc_order      = 2         # Boundary condition order for multigrid
   poisson_multifluid_gmg.gmg_bottom_solver = bicgstab  # Bottom solver type. 'simple', 'bicgstab', or 'gmres'
   poisson_multifluid_gmg.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)
   poisson_multifluid_gmg.gmg_cycle         = vcycle    # Cycle type. Only 'vcycle' supported for now
   poisson_multifluid_gmg.gmg_relax_type    = gsrb      # Relaxation type. 'jacobi', 'gauss_seidel', or 'gsrb'

Adjusting output
----------------

The user may plot the potential, the space charge, the electric, and the GMG residue as follows:

.. code-block:: bash

   poisson_multifluid_gmg.plt_vars  = phi rho E res     # Plot variables. Possible vars are 'phi', 'rho', 'E', 'res'

   
