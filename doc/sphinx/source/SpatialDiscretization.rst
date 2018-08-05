.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

PlasmaC uses spatial discretization provided by Chombo. In PlasmaC, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries; there is no need for voxelization of the domain.

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit. 
