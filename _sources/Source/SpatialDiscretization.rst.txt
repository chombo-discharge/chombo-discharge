.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

Cartesian AMR
_____________

``chombo-discharge`` uses patch-based structured adaptive mesh refinement (AMR) provided by ``Chombo`` :cite:`chombo`.
In patch-based AMR the domain is subdivided into a collection of hierarchically nested grid levels, see :numref:`Fig:PatchBasedAMR`.
With Cartesian AMR each patch is a Cartesian block of grid cells. 
A *grid level* is composed of a union of grid patches sharing the same grid resolution, with the additional requirement that the patches on a grid level are *non-overlapping*.
With AMR, such levels can be hierarchically nested; finer grid levels exist on top of coarser ones.
In patch-based AMR there are only a few fundamental requirements on how such grids are constructed.
For example, a refined grid level must exist completely within the bounds of it's parent level. 
In other words, grid levels :math:`l-1` and :math:`l+1` are spatially separated by a non-zero number of grid cells on level :math:`l`.

.. _Fig:PatchBasedAMR:
.. figure:: /_static/figures/PatchBasedAMR.png
   :width: 60%
   :align: center

   Cartesian patch-based refinement showing two grid levels.
   The fine-grid level lives on top of the coarse level, and consists of two patches (red and blue colors) with two layers of ghost cells (dashed lines and orange shaded region). 

The resolution on level :math:`l+1` is typically finer than the resolution on level :math:`l` by an integer (usually power of two).
However,

.. important::
   
   ``chombo-discharge`` only supports refinement factors of 2 and 4.

Embedded boundaries
___________________

``chombo-discharges`` uses an embedded boundary (EB) formulation for describing complex geometries.
With EBs, the Cartesian grid is directly intersected by the geometry.
This is fundamentally different from unstructured grid where one generates a volume mesh that conforms to the surface mesh of the input geometry.
Since EBs are directly intersected by the geometry, there is no fundamental need for a surface mesh for describing the geometry.
Moreover, Cartesian EBs have a data layout which remains (almost) fully structured.
The connectivity of neighboring grid cells is still trivially found by fundamental strides along the data rows/columns, which allows extending the efficiency of patch-based AMR to complex geometries.
Figure :numref:`Fig:ComplexPatches` shows an example of patch-based grid refinement for a complex surface.

.. _Fig:ComplexPatches:
.. figure:: /_static/figures/ComplexPatches.png
   :width: 50%
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit.

Since EBs are directly intersected by the geometry, pathological cases can arise where a Cartesian grid cell consists of multiple volumes.
One can easily envision this case by intersecting a thin body with a Cartesian grid, as shown in :numref:`Fig:MultiCells`.
This figure shows a thin body which is intersected by a Cartesian grid, and this grid is then coarsened.
At the coarsened level, one of the grid cells has two cell fragments on opposite sides of the body.
Such multi-valued cells (a.k.a *multi-cells*) are fundamentally important for EB applications.
Note that there is no fundamental difference between single-cut and multi-cut grid cells.
This distinction exists primarily due to the fact that if all grid cells were single-cut cells the entire EB data structure would fit in a Cartesian grid block (say, of :math:`N_x \times N_y \times N_z` grid cells).
Because of multi-cells, EB data structures are not purely Cartesian.
Data structures need to live on more complex graphs that describe support multi-cells and, furthermore, describe the cell connectivity.
Without multi-cells it would be impossible to describe most complex geometries.
It would also be extremely difficult to obtain performant geometric multigrid methods (which rely on this type of coarsening). 

.. _Fig:MultiCells:
.. figure:: /_static/figures/MultiCells.png
   :width: 50%
   :align: center

   Example of how multi-valued cells occur during grid coarsening.
   Left: Original grid.
   Right: Coarsened grid.


.. _Chap:GeometryRepresentation:

Geometry representation
_______________________

``chombo-discharge`` uses (approximations to) signed distance functions (SDFs) for describing geometries.
Signed distance fields are functions :math:`f: \mathbb{R}^3\rightarrow \mathbb{R}` that describe the distance from the object.
These functions are also *implicit functions*, i.e. :math:`f\left(\mathbf{x}\right)=0` describes the surface of the object, :math:`f\left(\mathbf{x}\right) > 0` decribes a point inside the object and :math:`f\left(\mathbf{x}\right) < 0` describes a point outside the object.

Many EB applications only use the implicit function formulation, but ``chombo-discharge`` requires (an approximation to) the signed distance field.
There are two reasons for this:

#. The SDF can be used for robustly load balancing the geometry generation with orders of magnitude speedup over naive approaches. 
#. The SDF is useful for resolving particle collisions with boundaries, using e.g. simple ray tracing of particle paths.

To illustrate the difference between an SDF and an implicit function, consider the implicit functions for a sphere at the origin with radius :math:`R`:

.. math::
   :nowrap:
   
   \begin{align}
   d_1\left(\mathbf{x}\right) &= R - \left|\mathbf{x}\right|, \\
   d_2\left(\mathbf{x}\right) &= R^2 - \mathbf{x}\cdot\mathbf{x}.
   \end{align}

Here, only :math:`d_1\left(\mathbf{x}\right)` is a signed distance function.    

In ``chombo-discharge``, SDFs can be generated through analytic expressions, constructive solid geometry, or by supplying polygon tesselation.
NURBS geometries are, unfortunately, not supported.
Fundamentally, all geometric objects are described using ``BaseIF`` objects from ``Chombo``, see :ref:`Chap:BaseIF`.

Constructive solid geometry (CSG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Constructive solid geometry can be used to generate complex shapes from geometric primitives.
For example, to describe the union between two SDFs :math:`d_1\left(\mathbf{x}\right)` and :math:`d_2\left(\mathbf{x}\right)`:

.. math::

   d\left(\mathbf{x}\right) = \textrm{min}\left(d_1\left(\mathbf{x}\right), d_2\left(\mathbf{x}\right)\right)

Note that the resulting is an implicit function but is *not* an SDF.
However, the union typically approximates the signed distance field quite well near the surface.
``Chombo`` natively supports many ways of performing CSG.

EBGeometry
^^^^^^^^^^

While functions like :math:`R - \left|\mathbf{x}\right|` are quick to compute, a polygon surface may consist of hundreds of thousands of primitives (e.g., triangles).
Generating signed distance function from polygon tesselations is quite involved as it requires computing the signed distance to the closest feature, which can be a planar polygon (e.g., a triangle), edge, or a vertex. 
``chombo-discharge`` supports such functions through the `EBGeometry <https://github.com/rmrsk/EBGeometry>`_ package.

.. warning::

   The signed distance function for a polygon surface is only well-defined if it is manifold-2, i.e. it is watertight and does not self-intersect.
   ``chombo-discharge`` should nonetheless compute the distance field as best as it can, but the final result may not make sense in an EB context. 

Searching through all features (faces, edge, vertices) is unacceptably slow, and ``EBGeometry`` therefore uses a bounding volume hierarchy for accelerating these searches.
The bounding volume hierarchy is top-down constructed, using a root bounding volume (typically a cube) that encloses all triangles.
Using heuristics, the root bounding volume is then subdivided into two separate bounding volumes that contain roughly half of the primitives each.
The process is then recursed downwards until specified recursion criteria are met.
Additional details are provided in the `EBGeometry documentation <https://rmrsk.github.io/EBGeometry/>`_.

.. figure:: /_static/figures/Armadillo.png
   :width: 50%
   :align: center

   Example of an SDF reconstruction and cut-cell grid from a surface tesselation in ``chombo-discharge``.

.. _Chap:GeometryGeneration:

Geometry generation
___________________

``Chombo`` approach
^^^^^^^^^^^^^^^^^^^

The default geometry generation method in ``Chombo`` is to locate cut-cells on the finest AMR level first and then generate the coarser levels cells through grid coarsening.
This will look through all cells on the finest level, so for a domain which is effectively :math:`N\times N\times N` cells there are :math:`\mathcal{O}\left(N^3\right)` implicit function queries (in 2D, the complexity is :math:`\mathcal{O}\left(N^2\right)`). 
Note that as :math:`N` becomes large, say :math:`N=10^5`, geometric queries of this type become a bottleneck.

``chombo-discharge`` pruning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``chombo-discharge`` has made modifications to the geometry generation routines in ``Chombo``, resolving a few bugs and, most importantly, using the signed distance function for load balancing the geometry generation step.
This modification to ``Chombo`` yields a reduction of the original :math:`\mathcal{O}\left(N^3\right)` scaling in ``Chombo`` grid generation to an :math:`\mathcal{O}\left(N^2\right)` scaling in ``chombo-discharge``.
Typically, we find that this makes geometry generation computationally trivial (in the sense that it is very fast compared to the simulation). 

To understand this process, note that the SDF satisfies the Eikonal equation

.. math::
   :nowrap:
      
   \begin{equation}
   \left|\nabla f\right| = 1, 
   \end{equation}
   
and so it is well-behaved for all :math:`\mathbf{x}`.
The SDF can thus be used to prune large regions in space where cut-cells don't exist. 
For example, consider a Cartesian grid patch with cell size :math:`\Delta x` and cell-centered grid points :math:`\mathbf{x}_{\mathbf{i}} = \left(\mathbf{i} + \mathbf{\frac{1}{2}}\right)\Delta x` where :math:`\mathbf{i} \in \mathbb{Z}^3` are grid cells in the patch, as shown in :numref:`Fig:Pruning`.
We know that cut cells do not exist in the grid patch if :math:`\left|f\left(\mathbf{x}_{\mathbf{i}}\right)\right| > \frac{1}{2}\Delta x` for all :math:`\mathbf{i}` in the patch.
One can use this to perform a quick scan of the SDF on a *coarse* grid level first, for example on :math:`l=0`, and recurse deeper into the grid hierarchy to locate cut-cells on the other levels. 
Typically, a level is decomposed into Cartesian subregions, and each subregion can be scanned independently of the other subregions (i.e. the problem is embarassingly parallel).
Subregions that can't contain cut-cells are designated as *inside* or *outside*, depending on the sign of the SDF.
There is no point in recursively refining these to look for cut-cells at finer grid levels, owing to the nature of the SDF they can be safely pruned from subsequent scans at finer levels. 
The subregions that did contain cut-cells are refined and decomposed into sub-subregions. 
This procedure recurses until :math:`l=l_{\text{max}}`, at which point we have determined all sub-regions in space where cut-cells can exist (on each AMR level), and pruned the ones that don't.
This process is shown in :numref:`Fig:Pruning`. 
Once all the grid patches that contain cut-cells have been found, these patches are distributed (i.e., load balanced) to the various MPI ranks for computing the discrete grid information.

.. _Fig:Pruning:
.. figure:: /_static/figures/Pruning.png
   :width: 75%
   :align: center

   Pruning cut-cells with the signed distance field.
   Red-colored grid patches are grid patches entirely contained inside the EB.
   Green-colored grid patches are entirely outside the EB, while blue-colored grid patches contain cut-cells.

The above load balancing strategy is very simple, and it reduces the original :math:`O(N^3)` complexity in 3D to :math:`O(N^2)` complexity (in 2D the complexity is reduced from :math:`O(N^2)` to :math:`O(N)`).
The strategy works for all SDFs although, strictly speaking, an SDF is not fundamentally needed.
If a well-behaved Taylor series can be found for an implicit function, the bounds on the series can also be used to infer the location of the cut-cells, and the same algorithm can be used.
For example, generating compound objects with CSG are typically sufficiently well behaved (provided that the components are SDFs). 
However, implicit functions like :math:`d\left(\mathbf{x}\right) = R^2 - \mathbf{x}\cdot\mathbf{x}` must be used with caution. 

.. _Chap:MeshGeneration:

Mesh generation
_______________

``chombo-discharge`` supports two algorithm for AMR grid generation:

#. The classical Berger-Rigoutsos algorithm :cite:`Berger1991`.
#. A *tiled* algorithm :cite:`Gunney2016`.
   
Both algorithms work by taking a set of flagged cells on each grid level and generating new boxes that cover the flags.
Only *properly nested* grids are generated, in which case two grid levels :math:`l-1` and :math:`l+1` are separated by a non-zero number of grid cells on level :math:`l`.
This requirement is not fundamentally required for quad- and oct-tree grids, but is nevertheless usually imposed. 
For patch based AMR, the rationale for this requirement is that stencils on level :math:`l+1` should should only reach into grid cells on levels :math:`l` and :math:`l+1`. 
For example, ghost cells on level :math:`l+1` can then be interpolated from data only on levels :math:`l` and :math:`l+1`.

Berger-Rigoutsos algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^

The Berger-Rigoustous grid algorithm is implemented in ``Chombo`` and is called by ``chombo-discharge``.
The classical Berger-Rigoustous algorithm is inherently serial in the sense that is collects the flagged cells onto each MPI rank and then generates the boxes, see :cite:`Berger1991` for implementation details. 
Typically, it is not used at large scale in 3D due to its memory consumption. 

.. _BRMeshRefine:
.. figure:: /_static/figures/BRMeshRefine.png
   :width: 25%
   :align: center

   Classical cartoon of patch-based refinement. Bold lines indicate entire grid blocks.

Tiled mesh refinement
^^^^^^^^^^^^^^^^^^^^^
``chombo-discharge`` also supports a tiled algorithm where the grid boxes on each block are generated according to a predefined tiled pattern.
If a tile contains a single tag, the entire tile is flagged for refinement.
The tiled algorithm produces grids that are visually similar to octrees, but is slightly more general since it also supports refinement factors other than 2 and is not restricted to domain extensions that are an integer factor of 2 (e.g. :math:`2^{10}` cells in each direction).
Moreover, the algorithm is extremely fast and has low memory consumption even at large scales. 

.. _TiledMeshRefine:
.. figure:: /_static/figures/TiledMeshRefine.png
   :width: 25%
   :align: center

   Classical cartoon of tiled patch-based refinement. Bold lines indicate entire grid blocks. 

.. _Chap:RefinementPhilosophy:

Cell refinement philosophy
__________________________

``chombo-discharge`` can flag cells for refinement using various methods:

#. Refine all embedded boundaries down to a specified refinement level.
#. Refine embedded boundaries based on estimations of the surface curvature in the cut-cells.
#. Manually add refinement flags (by specifying boxes where cells will be refined).
#. Physics-based or data-based refinement where the user fetches data from solver classes (e.g., discretization errors, the electric field) and uses that for refinement.

The first two cases are covered by the ``Driver`` class in ``chombo-discharge`` (see :ref:`Chap:Driver`). 
In the first case the ``Driver`` class will simply fetch arguments from an input script which specifies the refinement depth for the embedded boundaries. 
In the second case, the ``Driver`` class will visit every cut-cell and check if the normal vectors in neighboring cut-cell deviate by more than a specified threshold angle. 
Given two normal vectors :math:`\mathbf{n}` and :math:`\mathbf{n}^\prime`, the cell is refined if

.. math::

   \mathbf{n}\cdot\mathbf{n}^\prime \geq \cos\theta_c,

where :math:`\theta_c` is a threshold angle for grid refinent. 

The other two cases are more complicated, and are covered by the :ref:`Chap:GeoCoarsener` and :ref:`Chap:CellTagger` classes.
