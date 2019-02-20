.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

PlasmaC uses spatial discretization provided by Chombo. In PlasmaC, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries. 

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit. 

.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for PlasmaC follows that of Chombo. In Chombo, the geometries are generated from a function :math:`f(\mathbf{x})` that describes the level-set surface. This is done by first constructing a set of boxes that covers the finest AMR level. If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box. The box is allocated in full, so using a smaller box will reduce the memory consumption (but increase run time). Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function. There are no graphs in boxes that are all-covered or all-regular. Furthermore, geometric data describes by the graph only exists in the cut cells themselves, so that this data is truly sparse. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible. The finest AMR level always dominates the memory consumption for the EB mesh, and arbitrarily fine grids are not possible. Essentially, one must expect that one dense box is allocated for each box that intersects the graph at the finest AMR level. In other words, the available system memory must be sufficiently so that one can refine the entire embedded boundary at the finest AMR level. For example, if your finest AMR grid contains :math:`10^4` cut-cell boxes each of size :math:`64^3` (both are realistic numbers), the graph will contain roughly :math:`2.6\times 10^9` graph nodes (one node per cell in each cut-cell box). The memory consumption per graph various, but if we estimate this at :math:`512` bytes per node on average, then the graph memory consumption is :math:`1.34\,\textrm{TB}`. 

Hyperbolic discretization
-------------------------

The hyperbolic discretization in :verb:`PlasmaC` 
