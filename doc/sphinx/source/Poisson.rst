.. _Chap:PoissonSolver:
   
Poisson
=======

The `PlasmaC` field solver has a lot of supporting functionality, but essentially relies on only one critical function: Solving for the potential.
This is done by calling a class-specific function

.. code-block:: c++

   bool solve(MFAMRCellData& phi, const MFAMRCellData& rho, const EBAMRIVData& sigma);

where ``phi`` is the resulting potential that was computing with the space charge density ``rho`` and surface charge density ``sigma``.

Currently, only one field solver is implemented and this solver uses a geometric multigrid method for solving for the potential.
The solver supports three phases: electrodes, gas, and dielectric.
Domain boundary conditions for the solver must be set by the user through an input script, whereas the boundary conditions on internal surfaces are Dirichlet by default.
