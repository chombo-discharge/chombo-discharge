.. _Chap:Model:

The ``PlasmaC`` model
=====================

This chapter discusses the overall ``PlasmaC`` model, including the spatial decomposition and a summary of supported solvers.

Main functionality
------------------

In this section we summarize how components in ``PlasmaC`` are connected, such that users may understand more readily how the code is designed, and how it can be extended.

There are four major components `PlasmaC`:

1. A computational geometry which describes a level-set geometry consisting of electrodes and possibly also dielectrics.
   This functionality is encapsulated by :ref:`Chap:computational_geometry`. 
2. An AMR mesh which contains the grids, grid generation routines, and functionality for handling data coarsening and refinement.
   This functionality is encapsulated by :ref:`Chap:amr_mesh`.
3. A time stepper which advances the equations of motion (whatever they are).
   The :ref:`Chap:time_stepper` class is abstract.
   In order to actually use ``PlasmaC`` for anything the user must either write his own derived :ref:`Chap:time_stepper` class,
   or use one of the provided physics modules. 
4. A cell tagger which flags cells for refinement and coarsening.
   This functionality is encapsulated by the :ref:`Chap:cell_tagger` class. 

Instantiations of the above four classes are fed into the :ref:`Chap:driver` class which contains the functionality for running a simulation.
The reason for the above division of labor is that we have wanted to segregate responsibilities in order to increase flexibility.
For that reason, the computational geometry does not have any view of the actual AMR grids; it only contains the level-set functions and some meta-information (such as the permittivity of a dielectric).
Likewise, the :ref:`Chap:amr_mesh` class only acts a centralized repository of useful functions for AMR simulations.
These functions include algorithms for generating AMR grids, allocating data across AMR, and synchronizing AMR levels (e.g. interpolating ghost cells).

All the physics is encapsulated by the :ref:`Chap:time_stepper` class. 
This class will have direct ownership of all the solvers and the functions required to advance them over a time step.
Instantiations of the class will also contain the routines for setting up a simulation, e.g. instantiating solvers, setting up boundary conditions. 

The :ref:`Chap:driver` class is only responsible for *running* a simulation.
This class will call for regrids at certain intervals, call the :ref:`Chap:time_stepper` for writing plot and checkpoint data, and call for the :ref:`Chap:time_stepper` to advance the equations of motion.
In order to understand how ``PlasmaC`` runs a simulation, it will be useful to first understand how :ref:`Chap:driver` works. 

Simulation inputs
_________________

``PlasmaC`` simulations take their input from a single simulation file, possibly appended with overriding options on the command line.
Simulations may consist of several hundred possible switches for altering the behavior of a simulation, and physics models in ``PlasmaC`` are therefore equipped with Python setup tools that collect all such options in a single file.
Generally, these input parameters are fetched from the dependencies of each class or module in a simulation.
For example, all numerical solvers have independent adjustment of output.
The input parameters for each solver class is appended in a separate file named :file:`<solver>.options`.
For example, the input parameters for the default Poisson solver defined in :file:`/src/poisson/poisson_multifluid_gmg.H` is contained in a file :file:`/src/poisson/poisson_multifluid_gmg.options`.
This means that any time a user wishes to use such a solver, he may fetch all the input options *for that class* directly from the supplementary options file.

Simulation outputs
__________________

Mesh data from ``PlasmaC`` simulations is by default written to HDF5 files.
Users that wish to write or output other types of data must supply code themselves. 
   
.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

`PlasmaC` uses structured adaptive mesh refinement (SAMR provided by Chombo :cite:`ebchombo`.
SAMR exists in two separate categories, patch-based and tree-based AMR.
Patch-based AMR is the more general type and contain tree-based grids as a subset; they can use refinement factors other than 2, as well as accomodate anisotropic resolutions and non-cubic patches.
In patch-based AMR the domain is subdivided into a collection of hierarchically nested overlapping patches (or boxes).
Each patch is a rectangular block of cells which, in space, exists on a subdomain of the union of patches with a coarser resolution.
Patch-based grids generally do not have unique parent-children relations: A fine-level patch may have multiple coarse-level parents.
An obvious advantage of a patch-based approach is that entire Cartesian blocks are sent into solvers, and that the patches are not restricted to squares or cubes that align with the coarse-grid boundary.
A notable disadvantage is that additional logic is required when updating a coarse grid level from the overlapping region of a finer level.
Tree-based AMR use quadtree or octree data structures that describe a hierarchy of unique parent-children relations throughout the AMR levels: Each child has exactly one parent, whereas each parent has multiple children (4 in 2D, 8 in 3D).
In ``PlasmaC`` and Chombo, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4.
On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells.
This allows us to combine the efficiency of patch-based AMR with complex geometries.
However, there is significant overhead with the embedded boundary approach and, furthermore, arbitrarily complex geometries are not possible.

.. _Chap:MeshGeneration:

Mesh generation
_______________

`PlasmaC` offers two algorithm for AMR grid generation.
Both algorithms work by taking a set of flagged cells on each grid level and generating new boxes that cover the flags.
The first algorithm that we support is the classical Berger-Rigoustous grid algorithm that ships with Chombo, see the figure below.
The classical Berger-Rigoustous algorithm is serial-like in the sense that is collects the flagged cells onto each MPI rank and then generates the boxes.
The algorithm is typically not used at large scale because of its memory consumption. 

.. figure:: figures/amr.png
   :width: 240px
   :align: center

   Classical cartoon of patch-based refinement. Bold lines indicate entire grid blocks.

As an alternative, we also support a tiled algorithm where the grid boxes on each block are generated according to a predefined tiled pattern.
If a tile contains a single tag, the entire tile is flagged for refinement.
The tiled algorithm produces grids that are similar to octrees, but it is more general since it also supports refinement factors other than 2, and is not restricted to domain extensions that are an integer factor of 2 (e.g. :math:`2^{10}` cells in each direction). 

.. figure:: figures/tiled.png
   :width: 360px
   :align: center

   Classical cartoon of tiled patch-based refinement. Bold lines indicate entire grid blocks. 
	   
.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for ``PlasmaC`` follows that of Chombo. In Chombo, the geometries are generated from an implicit function :math:`f(\mathbf{x}) = 0` that describes the level-set surface. 

In `Chombo`, geometry generation is done by first constructing a set of boxes that covers the finest AMR level.
If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box.
The box is allocated in full, so using a smaller box will reduce the memory consumption, but since ghost cells are used there is a limitation to how much one can reduce the memory.
Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function.
There should be no graphs in boxes that are all-covered or all-regular. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible.
Memory consumption generally depends on the complexity of the geometry, and arbitrarily fine grids with cut-cell geometries are not possible.
Consider for example a cubic domain of :math:`(16384)^3` cells which is decomposed into :math:`(32)^3` cell size patches.
This yields :math:`(512)^3` possible patches in total.
Now consider that this domain is cut in half by a plane with normal vector :math:`\mathbf{n} = \hat{\mathbf{x}}`.
This surface will require allocation of :math:`512\times512\times 1` patches for the geometry.
If each patch is padded with 4 ghost cells, this yields :math:`512^2\times(40)^3 \approx 1.6\times 10^{10}` cells.
Inside each cell we must store volume fractions, area fractions, cell centroids positions and so one.
Although the surface is simple, the required memory easily ranges in the terabyte range. 

The default load-balancing for geometry generation in `Chombo` is an even division of the uniform finest-level grid among all the available.
This is a reasonable approach for porous media where the cut-cells distribute evenly through the computational domain, but the approach is not scalable to large domain sizes. 

To achieve scalable geometry generation, we have changed how `Chombo` generates the geometry generation on the various levels.
Our new approach first generates a map on a *coarse* level which is specified by the user.
On the specified level the domain is broken up into equal-sized chunks and cut-cell boxes are located.
Uncut and cut boxes are load balanced among the various ranks.
We then proceed towards the next finer level where the cut-cell boxes are identified by a refinement of the box distribution on the previous level.
Boxes that resulted from a refinement of the coarse level cut boxes are again broken up into equal-sized chunks, whereas the uncut boxes are not.
This is again followed by load-balancing of the cut boxes, and this process is repeated recursively down to the finest AMR level.
In essence, the geometry generation is load balanced based on where the cut cells are going to be.
For the user, he will be able to switch between the `Chombo` and ``PlasmaC`` approaches to geometry generation load balancing by flipping a flag in an input script.
