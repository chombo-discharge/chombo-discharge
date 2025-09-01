.. _Chap:VV:

Verification and validation
===========================

We strive to include convergence testing (verification), and in some cases comparison with various types of experimental (validation).
Below, we discuss our approach to spatial and temporal convergence testing.

.. _Chap:SpatialConvergence:

Spatial convergence
-------------------

Assume that we have some evolution problem which provides a solution on the mesh as :math:`\phi_{\mathbf{i}}^k\left(\Delta x, \Delta t\right)` where :math:`\Delta x` is a uniform grid resolution and :math:`\Delta t` is the time step used for evolving the state from :math:`t=0` to :math:`t = k\Delta t`.

To estimate the spatial order of convergence for the discretization we can use Richardson extrapolation to estimate the error, using the results from a finer grid resolution as the "exact" solution.
We solve the problem on grids with resolutions :math:`\Delta x_c` and a finer resolution :math:`\Delta x_f` and estimate the error in the coarse-grid solution as

.. math::

   E_{\mathbf{i}}^k\left(\Delta x_c\right) = \phi_{\mathbf{i}}^k\left(\Delta x_c,\Delta t\right) - \left\{\mathcal{A}_{\Delta x_f\rightarrow \Delta x_c}\left[\phi^k\left(\Delta x_f,\Delta t\right)\right]\right\}_{\mathbf{i}}.

where :math:`\mathcal{A}_{\Delta x_f\rightarrow \Delta x_c}` is an averaging operator which coarsens the solution from the fine grid (:math:`\Delta x_f`) to the coarse grid (:math:`\Delta x_c`).
The error norms are computed from :math:`E_{\mathbf{i}}^k\left(\Delta x_c; \Delta x_f\right)`.
Specifically:

.. math::

   L_\infty\left[\phi^k\left(\Delta x_c, \Delta t\right)\right] &= \max\left|E_{\mathbf{i}}^k\left(\Delta x_c\right)\right|, \\
   L_1\left[\phi^k\left(\Delta x_c, \Delta t\right)\right] &= \frac{1}{\sum_{\mathbf{i}}}\sum_{\mathbf{i}}\left|E_{\mathbf{i}}^k\left(\Delta x_c\right)\right|, \\
   L_2\left[\phi^k\left(\Delta x_c, \Delta t\right)\right] &= \sqrt{\frac{1}{\sum_{\mathbf{i}}}\sum_{\mathbf{i}}\left|E_{\mathbf{i}}^k\left(\Delta x_c\right)\right|^2}.

where the sums run over the grid points.

.. tip::
   
   If one has an exact solution :math:`\phi_e(\mathbf{x},T)` available, one can replace :math:`\phi^k_{\mathbf{i}}\left(\Delta x, \Delta t\right)` by :math:`\phi_e\left(\mathbf{i}\Delta x, k\Delta t\right)`.

.. _Chap:TemporalConvergence:   

Temporal convergence
--------------------

For temporal convergence we compute the errors in the same way as for the spatial convergence, replacing :math:`\Delta t` and :math:`\Delta x` as fixed parameters.
The solution error is computed as

.. math::
   
   E_{\mathbf{i}}^k\left(\Delta t_c\right) = \phi_{\mathbf{i}}^k\left(\Delta x,\Delta t_c\right) - \phi_{\mathbf{i}}^{k^\prime}\left(\Delta x,\Delta t_f\right),

where :math:`k\Delta t_c = T` and :math:`k^\prime \Delta t_f = T`.
Temporal integration errors are computed as

.. math::
   
   L_\infty\left[\phi^k\left(\Delta x, \Delta t_c\right)\right] &= \max\left|E_{\mathbf{i}}^k\left(\Delta t_c\right)\right|, \\
   L_1\left[\phi^k\left(\Delta x, \Delta t_c\right)\right] &= \frac{1}{\sum_{\mathbf{i}}}\sum_{\mathbf{i}}\left|E_{\mathbf{i}}^k\left(\Delta t_c\right)\right|, \\
   L_2\left[\phi^k\left(\Delta x, \Delta t_c\right)\right] &= \sqrt{\frac{1}{\sum_{\mathbf{i}}}\sum_{\mathbf{i}}\left|E_{\mathbf{i}}^k\left(\Delta t_c\right)\right|^2}.

.. tip::

   If an exact solution as available, one can replace :math:`\phi^{k^\prime}_{\mathbf{i}}\left(\Delta x, \Delta t_f\right)` by :math:`\phi_e\left(\mathbf{i}\Delta x, k^\prime\Delta t_f\right)`.
