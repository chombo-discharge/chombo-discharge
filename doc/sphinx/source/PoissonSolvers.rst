.. _Chap:PoissonSolvers:

Poisson solvers
---------------

Currently, we only have one available Poisson solver represented by :doxy:`poisson_multifluid_gmg <classpoisson__multifluid__gmg>`, which derives from a more general interface.

The current Poisson solver uses a second order multifluid discretization, and the Poisson equation is solved by using geometric multigrid. The interface into :doxy:`poisson_multifluid_gmg <classpoisson__multifluid__gmg>` is by means if modifying solver settings:

.. literalinclude:: links/poisson_multifluid_gmg.options

Here, ``gmg_pre_smooth`` denotes the number of relaxations in the downsweep of the V-cycle; ``gmg_post_smooth`` denotes the number of relaxations in the upsweed, and ``gmg_bottom_smooth`` denotes the number of relaxations before entering the bottom solver on the bottom of the V-cycle. The tolerance limit is given by ``gmg_tolerance``; the solver will exit multigrid if the zero residual has been decreased by this amount. Likewise, ``hang`` denotes convergence rate limit; the solver will exit multigrid if the convergence rate drops below the ``hang`` threshold.

The ``bottom_drop`` flag specifies when :doxy:`poisson_multifluid_gmg <classpoisson__multifluid__gmg>` should drop to the bottom solver. The value in this setting specifies that we drop to the bottom solver when the number of cells in the :math:`x` coordinate equals this number. For example, users may specify a coarsest domain of :math:`128 \times 128`, but decide that multigrid should not coarsen all the way down. If the user wants to enter the bottom solver routine when the domain equals :math:`8\times8` cells, he must set ``bottom_drop`` equal to 8. For complex geometries, it does not always pay off to coarsen all the way down (this is double true when the domain also contains dielectrics). Typically, optimization of the Poisson solver for new simulation cases involves examination of optimal values of the number of smoothings and ``bottom_drop``.

At the bottom of the V-cycle, the user may specify the desired bottom solver type. For multifluid cases (i.e. with dielectrics), we only have relaxations solvers on the bottom. To use these solvers, the user must set ``gmg_bottom_solver`` to *simple*. For single-fluid cases (i.e. geometries with only electrodes), the user can use a more exact biconjugate stabilized gradient solver. To use this solver, set ``gmg_bottom_solver`` to *bicgstab*. It (almost) always pays off to use the gradient solver on the V-cycle bottom since it is more exact than the relaxation sweep. 
