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
Domain boundary conditions for the solver must be set by the user through an input script, whereas the boundary conditions on internal surfaces are Dirichlet by default.

Setting up the solver
---------------------

In order to set up the solver, one must provide the
