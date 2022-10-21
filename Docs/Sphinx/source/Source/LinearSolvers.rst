.. _Chap:LinearSolvers:

Linear solvers
==============

.. _Chap:Helmholtz:

Helmholtz equation
------------------

The Helmholtz equation is represented by

.. math::

   \alpha a\left(\mathbf{x}\right)\Phi + \beta\nabla\cdot\left[b\left(\mathbf{x}\right)\nabla\Phi\right] = \rho
   
where :math:`\alpha` and :math:`\beta` are constants and :math:`a\left(\mathbf{x}\right)` and :math:`b\left(\mathbf{x}\right)` are spatially dependent and piecewise smooth.

To solve the Helmholtz equation, it is solved in the form

.. math::

   \kappa L\Phi = \kappa\rho,

where :math:`L` is the Helmholtz operator above.
The preconditioning by the volume fraction :math:`\kappa` is done in order to avoid the small-cell problem encountered in finite-volume discretizations on EB grids.

Discretization and fluxes
_________________________

The Helmholtz equation is solved by assuming that :math:`\Phi` lies on the cell-center.
The :math:`b\left(\mathbf{x}\right)`-coefficient lies on face centers and EB faces. 
In the general case the cell center might lie inside the embedded boundary, and the cell-centered discretization relies on the concept of an extended state.
Thus, :math:`\Phi` does not satisfy a discrete maximum principle.

.. _Fig:HelmholtzFluxes:
.. figure:: /_static/figures/CutCell.png
   :width: 40%
   :align: center

   Location of fluxes for finite volume discretization.

The finite volume update require fluxes on the face centroids rather than the centers.
These are constructed by first computing the fluxes to second order on the face centers, and then interpolating them to the face centroids.
For example, the flux :math:`F_3` in the figure above is

.. math::

   F_3 = \beta b_{i,j+1/2}\frac{\Phi_{i,j+1} - \Phi_{i,j}}{\Delta x}.

Boundary conditions
___________________

The finite volume discretization of the Helmholtz equation require fluxes through the EB and domain faces. 
Below, we discuss how these are implemented.

.. note::

   ``chombo-discharge`` supports spatially dependent boundary conditions

Neumann
^^^^^^^
Neumann boundary conditions are straightforward since the flux through the EB or domain faces are specified directly.

From the above figure, the fluxes :math:`F_{\textrm{EB}}` and :math:`F_{\textrm{D}}` are specified.

.. _Chap:LinearSolverDirichletBC:

Dirichlet
^^^^^^^^^

Dirichlet boundary conditions are more involved since only the value at the boundary is prescribed, but the finite volume discretization requires a flux. 
On the domain boundaries the fluxes are face-centered and we therefore use finite differencing for obtaining a second order accurate approximation to the flux at the boundary.

On the embedded boundaries the flux is more complicated to compute, and requires us to compute an approximation to the normal gradient :math:`\partial_n\Phi` at the boundary.
Our approach is to approximate this flux by expanding the solution as a polynomial around a specified number of grid cells.
By using more grid cells than there are unknown in the Taylor series, we formulate an over-determined system of equations up to some specified order.
As a first approximation we include only those cells in the quadrant or half-space defined by the normal vector, see :numref:`Fig:GradientReconstruction`. 
If we can not find enough equations, the strategy is to 1) drop order and 2) include all cells around the cut-cell.

.. _Fig:GradientReconstruction:
.. figure:: /_static/figures/GradientReconstruction.png
   :width: 40%
   :align: center

   Examples of neighborhoods (quadrant and half-space) used for gradient reconstruction on the EB. 

Once the cells used for the gradient reconstruction have been obtained, we use weighted least squares to compute the approximation to the derivative to specified order (for details, see :ref:`Chap:LeastSquares`). 
The result of the least squares computation is represented as a stencil:

.. _Eq:DirichletElliptic:
.. math::

   \frac{\partial\Phi}{\partial n} = w_{\textrm{B}}\Phi_{\textrm{B}} + \sum_{\mathbf{i}} w_{\mathbf{i}}\Phi_{\mathbf{i}},

where :math:`\Phi_{\textrm{B}}` is the value on the boundary, the :math:`w` are weights for grid points :math:`\mathbf{i}`, and the sum runs over cells in the domain.

Note that the gradient reconstruction can end up requiring more than one ghost cell layer near the embedded boundaries.
For example, :numref:`Fig:StencilRegion` shows a typical stencil region which is built when using second order gradient reconstruction on the EB.
In this case the gradient reconstruction requires a stencil with a radius of 2, but as the cut-cell lies on the refinement boundary the stencil reaches into two layers of ghost cells.
For the same reason, gradient reconstruction near the cut-cells might require interpolation of corner ghost cells on refinement boundaries. 

.. _Fig:StencilRegion:
.. figure:: /_static/figures/StencilRegion.png
   :width: 40%
   :align: center

   Example of the region of a second order stencil for the Laplacian operator with second order gradient reconstruction on the embedded boundary. 

Robin
^^^^^

Robin boundary conditions are in the form

.. math::

   A\partial_n\Phi + B\Phi = C,

where :math:`A`, :math:`B`, and :math:`C` are constants.
This boundary conditions is enforced through the flux

.. math::

   \partial_n\Phi = \frac{1}{A}\left(C - B\Phi\right),

which requires an evaluation of :math:`\Phi` on the domain boundaries and the EB.

For domain boundaries we extrapolate the cell-centered solution to the domain edge, using standard first order finite differencing.

On the embedded boundary, we approximate :math:`\Phi\left(\mathbf{x}_{\text{EB}}\right)` by linearly interpolating the solution with a least squares fit, using cells which can be reached with a monotone path of radius one around the EB face (see :ref:`Chap:LeastSquares` for details).
The Robin boundary condition takes the form

.. math::

   \partial_n\Phi = \frac{C}{A} - \frac{B}{A}\sum_{\mathbf{i}} w_{\mathbf{i}}\Phi_{\mathbf{i}}.
   
Currently, we include the data in the cut-cell itself in the interpolation, and thus also use unweighted least squares. 

.. _Chap:MultigridInterpolation:

Ghost cell interpolation
________________________

With AMR, multigrid requires ghost cells on the refinement boundary.
The interior stencils for the Helmholtz operator have a radius of one and thus only require a single layer of ghost cells (and no corner ghost cells).
These ghost cells are filled using a finite-difference stencil, see :numref:`Fig:MultigridInterpolation`.

.. _Fig:MultigridInterpolation:
.. figure:: /_static/figures/MultigridInterpolation.png
   :width: 40%
   :align: center

   Standard finite-difference stencil for ghost cell interpolation (open circle).
   We first interpolate the coarse-grid cells to the centerline (diamond).
   The coarse-grid interpolation is then used together with the fine-grid cells (filled circles) for interpolation to the ghost cell (open circle). 

Embedded boundaries introduce many pathologies for multigrid:

1. Cut-cell stencils may have a large radius (see :numref:`Fig:StencilRegion`) and thus require more ghost cell layers.
2. The EBs cut the grid in arbitrary ways, leading to multiple pathologies regarding cell availability. 

The pathologies mean that standard finite differencing fails near the EB, mandating a more general approach.
Our way of handling ghost cell interpolation near EBs is to reconstruct the solution (to specified order) in the ghost cells, using the available cells around the ghost cell (see :ref:`Chap:LeastSquares` for details). 
As per conventional wisdom regarding multigrid interpolation, this reconstruction does *not* use coarse-level grid cells that are covered by the fine level.

Figure :numref:`Fig:EBMultigridInterpolation` shows a typical interpolation stencil for the stencil in :numref:`Fig:StencilRegion`.
Here, the open circle indicates the ghost cell to be interpolated, and we interpolate the solution in this cell using neighboring grid cells (closed circles).
For this particular case there are 10 nearby grid cells available, which is sufficient for second order interpolation (which requirse at least 6 cells in 2D). 
   
.. _Fig:EBMultigridInterpolation:
.. figure:: /_static/figures/EBMultigridInterpolation.png
   :width: 40%
   :align: center

   Multigrid interpolation for refinement boundaries away from and close to an embedded boundary.

.. note::

   ``chombo-discharge`` implements a fairly general ghost cell interpolation scheme near the EB. The ghost cell values can be reconstructed to specified order (and with specified least squares weights). 
   

Relaxation methods
__________________

The Helmholtz equation is solved using multigrid, with various smoothers available on each grid level.
The currently supported smoothers are:

1. Standard point Jacobi relaxation. 
2. Red-black Gauss-Seidel relaxation in which the relaxation pattern follows that of a checkerboard. 
3. Multi-colored Gauss-Seidel relaxation in which the relaxation pattern follows quadrants in 2D and octants in 3D. 

Users can select between the various smoothers in solvers that use multigrid.

.. note::

   Multi-colored Gauss-Seidel usually provide the best convergence rates.
   However, the multi-colored kernels are twice as expensive as red-black Gauss-Seidel relaxation in 2D, and four times as expensive in 3D. 


Multiphase Helmholtz equation
-----------------------------

``chombo-discharge`` also supports a *multiphase version* where data exists on both sides of the embedded boundary.
The most common case is that involving discontinuous coefficients, e.g. for

.. math::

   \nabla\cdot\left[b\left(\mathbf{x}\right)\nabla\Phi\left(\mathbf{x}\right)\right] = 0. 

where :math:`b\left(\mathbf{x}\right)` is only piecewise constant.

Jump conditions
_______________

For the case of discontinous coefficients there is a jump condition on the interface between two materials:

.. math::

   b_1\partial_{n_1}\Phi + b_2\partial_{n_2}\Phi = \sigma,

where :math:`b_1` and :math:`b_2` are the Helmholtz equation coefficients on each side of the interface, and :math:`n_1 = -n_2` are the normal vectors pointing away from the interface in each phase.
:math:`\sigma` is a jump factor.

.. _Fig:JumpCondition:
.. figure:: /_static/figures/JumpCondition.png
   :width: 40%
   :align: center

   Example of cells and stencils that are involved in discretizing the jump condition. Open and filled circles indicate cells in separate phases.

Discretization
______________

To incorporate the jump condition in the Helmholtz discretization, we use a gradient reconstruction to obtain a solution to :math:`\Phi` on the boundary, and use this value to impose a Dirichlet boundary condition during multigrid relaxation.
Recalling the gradient reconstruction :math:`\frac{\partial\Phi}{\partial n} = w_{\textrm{B}}\Phi_{\textrm{B}} + \sum_{\mathbf{i}} w_{\mathbf{i}}\Phi_{\mathbf{i}}`, the matching condition (see :numref:`Fig:JumpCondition`) can be written as

.. math::

   b_1\left[w_{\textrm{B},1}\Phi_{\textrm{B}} + \sum_{\mathbf{i}} w_{\mathbf{i},1}\Phi_{\mathbf{i},1}\right] + b_2\left[w_{\textrm{B},2}\Phi_{\textrm{B}} + \sum_{\mathbf{i}} w_{\mathbf{i},2}\Phi_{\mathbf{i},2}\right] = \sigma.

This equation can be solved for the boundary value :math:`\Phi_{\textrm{B}}`, which can then be used to compute the finite-volume fluxes into the cut-cells.

.. note::

   For discontinuous coefficients the gradient reconstruction on one side of the EB does not reach into the other (since the solution is not differentiable across the EB).

AMRMultiGrid
------------

``AMRMultiGrid`` is the ``Chombo`` implementation of the Martin-Cartwright multigrid algorithm.
It takes an "operator factory" as an argument, and the factory can generate objects (i.e., operators) that encapsulate the discretization on each AMR level.

``chombo-discharge`` runs its own operator, and the user can use either of:

1. ``EBHelmholtzOpFactory`` for single-phase problems.
2. ``MFHelmholtzOpFactory`` for multi-phase problems.

The source code for these are located in :file:`$DISCHARGE_HOME/Source/Elliptic`.

Bottom solvers
______________

Chombo provides (at least) three bottom solvers which can be used with ``AMRMultiGrid``.

1. A regular smoother (e.g., point Jacobi).
2. A biconjugate gradient stabilized method (BiCGStab)
3. A generalized minimal residual method (GMRES).

The user can select between these for the various solvers that use multigrid.    
