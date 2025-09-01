.. _Chap:AmrMesh:

AmrMesh
=======

``AmrMesh`` handles virtually the entire AMR and cut-cell infrastructure in ``chombo-discharge``, and has the following responsibility:

#. Generate grids form a set of cell tags (potentially more than one set of grids).
#. Load balance the grids if necessary.
#. Instantiate the cut-cell information on grids.
#. Provide a mask which tells the user if a cell is covered by a finer grid.
#. Generate operators that are required for handling AMR data, e.g.,
   
   * Ghost cell interpolators.
   * Coarsening operators.
   * Stencils for interpolation and extrapolation near the embedded boundaries.
   * Fine-to-coarse grid interpolators.
   * Infrastructure for particle-mesh operators.
   * ... and many others.

In addition to these, ``AmrMesh`` contains function for actually allocating data with user-defined centering (e.g., cell, face, EB) on a specified :ref:`Chap:Realm`.
It also has responsibility for allocating particle containers.

One thing that ``AmrMesh`` is *not*, is a numerical discretization holder for PDEs, which is a responsibility that is generally deferred to solvers.
In some cases ``AmrMesh`` certainly holds operators that have an underlying discretization, e.g., gradient operators which have a specific discretization based on least squares reconstruction around the embedded boundaries.
But in virtually all cases one should simply think of these operators as AMR operators that handle specific types of common data operations on mesh and particle data. 

``AmrMesh`` is an integral part of ``chombo-discharge``, and users will never have the need to modify it unless they are implementing something entirely new.
The behavior of :ref:`Chap:AmrMesh` is modified through its available input parameters.

.. tip::

   Here is the `AmrMesh C++ API <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classAmrMesh.html>`_

Class options
-------------

Various class options are available for adjusting the behavior of the ``Driver`` class.
Below, we include the current template options file for the ``Driver`` class.

.. literalinclude:: ../../../../Source/AmrMesh/CD_AmrMesh.options
   :emphasize-lines: 9-15,21-26
   :caption: Template input options for the ``AmrMesh`` class. Runtime adjustable options are highlighted.   

For users, modifying the behavior of ``AmrMesh`` is comparatively easy.
We have listed the most commonly adjusted grid parameters above, which include the specification of the physical domain, the maximum number of AMR levels that are permitted, as well as how grids are generated.

Coarse-grid decomposition
-------------------------

``AmrMesh.lo_corner`` and ``AmrMesh.hi_corner`` are the physical corners of the simulation domain.
``AmrMesh.coarsest_domain`` is the number of grid cells on the coarsest grid level (i.e., without AMR).
It is important that cell sizes are uniform, so one must always have :math:`\Delta x = \Delta y = \Delta z`.
Usually, this means that ``AmrMesh.lo_corner``, ``AmrMesh.hi_corner``, and ``AmrMesh.coarsest_domain`` must all be consistently defined. Moreover, it is normally desirable to make ``AmrMesh.coarsest_domain`` a factor of 2 (e.g., 64, 128, 256, etc.), since this permits arbitrarily deep multigrid coarsening.
This is not a requirement, however, although we do note that ``AmrMesh.coarsest_domain`` *must* be divisible by ``AmrMesh.blocking_factor``.

Domain decomposition
--------------------

With Cartesian AMR, each grid level is decomposed into grid blocks of constant size, or sizes that potentially vary between some min/max size along each dimension.
In ``chombo-discharge`` this is encapsulated by ``AmrMesh.blocking_factor``, which is the smallest grid box that can be generated when meshing the domain.
Likewise, ``AmrMesh.max_box_size`` is the maximum box size that can be produced, but usage of a constant box size is common an increased requirement in ``chombo-discharge``.

.. tip::

   Use fixed box sizes where ``AmrMesh.blocking_factor`` and ``AmrMesh.max_box_size`` are the same.

The flag ``AmrMesh.max_ebis_box`` indicates essentially the blocking factor when generating the cut-cell information at the start of a simulation.
It may happen for very large simulations that one has to increase the box size (e.g., to 32) in order to trim grid metadata.

``AmrMesh.ref_rat`` determines the refinement ratio between grid levels, and factors of 2 and 4 are supported.
We want to point out that mixed refinement factors are supported.

.. tip::

   If ``AmrMesh.max_amr_depth`` is greater than the number of refinement ratios specified in ``AmrMesh.ref_rat``, ``chombo-discharge`` will automatically fill in the remaining refinement ratios (padding with the last entry). I.e., one obtains factor 2 refinement every if ``AmrMesh.ref_rat = 2``.

Two gridding algorithms are supported, called Tiled mesh refinement the classical Berger-Rigoutsos refinement algorithm.
These are discussed in :ref:`Chap:MeshGeneration`, and the user can specify which one to use by setting ``AmrMesh.grid_algorithm`` accordingly.
In general, the tiled algorithm is exceedingly more performant at larger scales.

Ghost cells
-----------

It is normally not necessary to adjust the number of ghost cells in ``chombo-discharge`` simulations since most discretizations require at most 2 ghost cells.
Substantial efforts have been made to avoid increasing the required number of ghost cells.
Regardless, the number of ghost cells can be adjusted by setting ``AmrMesh.num_ghost`` to a specified value.
A companion parameter is ``AmrMesh.eb_ghost``, which specifies the maximum number of ghost cells that are used when computing the cut-cell discretization.

.. warning::
   
   One must always have ``AmrMesh.eb_ghost`` > ``AmrMesh.num_ghost``.

Finally, it is possible to evaluate implicit functions on the mesh (e.g., for figuring out how far a cell is from a boundary).
It can be useful to include a larger ghost region in these data holders, which is adjusted by ``AmrMesh.lsf_ghost``.

Multigrid interpolation
-----------------------

Multigrid interpolation is done using least squares reconstruction, as discussed in :ref:`Chap:MultigridInterpolation`.
The user can set the order, radius, and least squares weighting for this interpolation by setting ``AmrMesh.mg_interp_order``, ``AmrMesh.mg_interp_radius``, and ``AmrMesh.mg_interp_weight``.
Note that specifying ``AmrMesh.mg_interp_radius`` > ``AmrMesh.eb_ghost`` is sure to lead to a run-time error (possibly a segfault).

Interpolation near the EB
-------------------------

Most data in ``chombo-discharge`` is cell-centered, although data can be interpolated or extrapolated to cell centroids and EB centroids.
The interpolation type is done by setting ``AmrMesh.centroid_interp`` and ``AmrMesh.eb_interp`` accordingly.
Currently, we support the following options:

#. ``constant`` i.e., use the cell centered value.
#. ``linear``, using bi/tri-linear interpolation.
#. ``taylor``, using the Taylor series evaluated at the cell center.
#. ``lsq``, using a first order unweighted least squares polynomial.
#. ``pwl``, using piecewise linear reconstruction.
#. ``minmod``, using a minmod slope limiter.
#. ``superbee``, using a superbee slope limiter.
#. ``monotonized_central``, using a van Leer slope limiter.

Typically, one can just leave this option at ``minmod``.


