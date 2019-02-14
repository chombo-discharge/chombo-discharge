.. _Chap:Equations:

The ``PlasmaC`` equation set
========================

``PlasmaC`` aims at being a moderately flexible framework for fluid plasma simulations. There are several abstractions in place that ensure that the code covers non-trivial geometries, multiple time stepping schemes, and sophisticated plasma-kinetic couplings. The equation set that ``PlasmaC`` (currently) solves is

.. math::
   :nowrap:

   \begin{align}
   &\nabla\cdot\left(\epsilon_r\nabla\cdot\Phi\right) = -\frac{\rho}{\epsilon_0}, \\[1ex]
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\[1ex]
   &\kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},\\[1ex]
   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n\right) = S.
   \end{align}

This must be supported by additional boundary conditions on electrodes and insulating surfaces (on which the surface charge :math:`\sigma` lives). The number of advected species and radiative transport equations is arbitrary; the user can select the coupling through our interfaces: He can even turn off radiative transport altogether.

The coupling that is (currently) available in ``PlasmaC`` is

.. math::
   :nowrap:

   \begin{align}
      \epsilon_r &= \epsilon_r(\mathbf{x}), (\textrm{can additionally be discontinuous}), \\[1ex]
      \mathbf{v} &= \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      D &= \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      S &= S(t, \mathbf{x}, \mathbf{E}, \nabla\mathbf{E}, n, \nabla n, \Psi), \\[1ex]
      \eta &= \eta\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      F &= F(t, \mathbf{x}, \mathbf{E}, n),
   \end{align}

where :math:`F` is the boundary flux on insulators or electrodes (which must be separately implemented). 


``PlasmaC`` works by embedding the equations above into an abstract C++ framework that the user must implement or reuse existing pieces of, and then compile into a *mini-application*. For most users, this will mostly include implementing a new geometry, a new plasma-kinetic scheme, or new functions for deciding when to coarsen or refine a certain spatial region. It is our goal that the user does not need to worry about temporal or spatial discretization of these equations. 

.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

``PlasmaC`` uses spatial discretization provided by Chombo. In ``PlasmaC``, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries. 

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit. 

.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for ``PlasmaC`` follows that of Chombo. In Chombo, the geometries are generated from a function :math:`f(\mathbf{x})` that describes the level-set surface. This is done by first constructing a set of boxes that covers the finest AMR level. If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box. The box is allocated in full, so using a smaller box will reduce the memory consumption (but increase run time). Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function. There are no graphs in boxes that are all-covered or all-regular. Furthermore, geometric data describes by the graph only exists in the cut cells themselves, so that this data is truly sparse. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible. The finest AMR level always dominates the memory consumption for the EB mesh, and arbitrarily fine grids are not possible. Essentially, one must expect that one dense box is allocated for each box that intersects the graph at the finest AMR level. In other words, the available system memory must be sufficiently so that one can refine the entire embedded boundary at the finest AMR level. For example, if your finest AMR grid contains :math:`10^4` cut-cell boxes each of size :math:`64^3` (both are realistic numbers), the graph will contain roughly :math:`2.6\times 10^9` graph nodes (one node per cell in each cut-cell box). The memory consumption per graph various, but if we estimate this at :math:`512` bytes per node on average, then the graph memory consumption is :math:`1.34\,\textrm{TB}`. 
