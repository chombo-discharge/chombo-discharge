.. _Chap:Model:

The ``PlasmaC`` code
====================

The ``PlasmaC`` code is a loosely coupled code targeted at solving plasma problems. 
The code uses an embedded boundary (EB) adaptive mesh refinement (AMR) formalism where the grids frequently change and are adapted to the solution as simulations progress.
By design, the ``PlasmaC`` does not subcycle and all the grids are advanced using the same time step.
``PlasmaC`` also supports the concept of a :ref:`Chap:Realm`, which in short means that ``PlasmaC`` supports using one set of grids for Eulerian solvers and a different set of grids for Lagrangian solvers. 

The core functionality is centered around a set of solvers, for example a Poisson solver and a convection-diffusion-reaction solver, and then using the built-in solver functionality to advance the equations of motion.
This is done through a class ``time_stepper``, which is an abstract class that advances the equations of motion within the ``PlasmaC`` framework.
The ``time_stepper`` can instantiate an arbitrary number of solvers, and allows developers to use a fairly high-level description of their problem.
For example, all *solvers* have functions like ``write_plot_data(...)`` that the user may use within the ``time_stepper`` output routines.

Although many abstractions are in place so that user can describe a new set of physics, or write entirely new solvers into ``PlasmaC`` and still use the EBAMR formalism, ``PlasmaC`` also provide some physics modules for describing various types of problems.
These modules reside in :file:`/physics` and they are intended to both be problem-solving physics modules, and as well acting like benchmarks, regression tests, and examples for extension to new types of physics modules in the future. 


Main functionality
------------------

The main functionality in ``PlasmaC`` is centered around the concept of a :ref:`Chap:solver`. 

In this section we summarize how components in ``PlasmaC`` are connected, such that users may understand more readily how the code is designed. 

There are four major components `PlasmaC`:

1. A computational geometry which describes a level-set geometry consisting of electrodes and possibly also dielectrics.
   This functionality is encapsulated by :ref:`Chap:computational_geometry`. 
2. An AMR mesh which contains the grids, grid generation routines, and functionality for handling data coarsening and refinement.
   This functionality is encapsulated by :ref:`Chap:amr_mesh`.
   The ``amr_mesh`` class acts as a centralized repository for grid generation, performing AMR operations like filling ghost cells, allocating data over the AMR hierarchy and so on.
   ``amr_mesh`` is a standalone class - it has no view over the rest of ``PlasmaC``. 
3. A time stepper which advances the equations of motion (whatever they are).
   This class has been made abstract with a public interface that is used by the ``driver`` class (see below). 
   In order to actually use ``PlasmaC`` for anything the user must either write his own derived :ref:`Chap:time_stepper` class, or use one of the pre-defined physics modules.
   The ``time_stepper`` is purposefully quite general so that the whole ``PlasmaC`` framework can be set up to solve completely new sets of equations without affecting the rest of the framework. 
4. A cell tagger which flags cells for refinement and coarsening.
   This functionality is encapsulated by the :ref:`Chap:cell_tagger` class and it, too, is abstract. 

Instantiations of the above four classes are fed into the :ref:`Chap:driver` class which contains calling functions for generation the geometry, having the time integrator to perform and time step, performing, I/O, setting checkpoint/restart and so on. 
The reason for the above division of labor is that we have wanted to segregate responsibilities in order to increase flexibility.
For that reason, the computational geometry does not have any view of the actual AMR grids; it only contains the level-set functions and some meta-information (such as the permittivity of a dielectric).
Likewise, the :ref:`Chap:amr_mesh` class only acts a centralized repository of useful functions for AMR simulations.
These functions include algorithms for generating AMR grids, allocating data across AMR, and synchronizing AMR levels (e.g. interpolating ghost cells).

All the physics is encapsulated by the :ref:`Chap:time_stepper` class. 
This class will have direct ownership of all the solvers and the functions required to advance them over a time step.
Instantiations of the class will also contain the routines for setting up a simulation, e.g. instantiating solvers, setting up boundary conditions.
Typically, implementation new physics consists of writing a new ``time_stepper`` class that allocates the relevant solvers, and then implement the time integration algorithms that advances them.
The folder :file:`/physics` contains implementation of a few different physics modules.
Since problems within a physics module tend to be conceptually similar, all of these modules also have a Python setup script so that users can quickly set up new types of problems within the same module. 

The :ref:`Chap:driver` class is only responsible for *running* a simulation, and it uses ``time_stepper`` to do so.
The :ref:`Chap:driver` class will call for regrids at certain intervals, call the :ref:`Chap:time_stepper` for writing plot and checkpoint data, and also call for the :ref:`Chap:time_stepper` to advance the equations of motion through a function ``advance(...)``. 
In order to understand how ``PlasmaC`` runs a simulation, it will be useful to first understand how :ref:`Chap:driver` works.

Solvers
_______

Various solvers are implemented in ``PlasmaC``, see :ref:`Chap:SupportedSolvers`.
All solvers are designed to run through the ``time_stepper`` class.
Therefore, in order to run only a single solver (e.g. advection-diffusion or Poisson), one must have a ``time_stepper`` implementation that allocates the appropriate solver, sets it up, and runs it.
Currently, there are separate physics modules for each type of solver such that users may see how they are set up and run.
These are located in :file:`/physics/`.

The solvers may be abstract or non-abstract.
All solvers that are *not* abstract are supplemented by an options file that contain all the possible run-time configurations that can be made to the solver.
Such options can include multigrid parameters, how to handle particle deposition with refinement boundaries, slope limiters, etc.
For example, all numerical solvers have independent adjustment of output.
The input parameters for each solver class is included in a separate file named :file:`<solver>.options` that resides in the same folder as the solver.
For example, the input parameters for the default Poisson solver defined in :file:`/src/poisson/poisson_multifluid_gmg.H` is contained in a file :file:`/src/poisson/poisson_multifluid_gmg.options`.

Simulation inputs
_________________

``PlasmaC`` simulations take their input from a single simulation input file, possibly appended with overriding options on the command line.
Simulations may consist of several hundred possible switches for altering the behavior of a simulation, and physics models in ``PlasmaC`` are therefore equipped with Python setup tools that collect all such options in a single file.
Generally, these input parameters are fetched from the options file of each class that is used in a simulation.
Simulation options usually consist of a prefix, a suffix, and a configuration value.
For example, the configuration options that adjusts the number of time steps that will be run in a simulation is

.. code-block:: bash

   driver.max_steps = 100



Simulation outputs
__________________

Mesh data from ``PlasmaC`` simulations is by default written to HDF5 files.
Users that wish to write or output other types of data must supply code themselves.

In addition to plot files, MPI ranks can output information to separate files so that the simulation progress can be tracked.
See :ref:`Chap:Control` for details. 
This is also useful for debugging purposes. 
   
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
The geometric data in the box is allocated sparsely so that memory consumption due to EB information storage is typically not very high. 
In general, there should be no graphs in boxes that are all-covered or all-regular. 

When EB information is first generated across the AMR hierarchy, one begins by computing the information on the finest grid level.
From there, coarser levels are generated through *coarsening* of the fine-information data. 
The default load-balancing for geometry generation in `Chombo` is an even division of the grid level among the ranks. 
This is a reasonable approach for porous media where the cut-cells distribute evenly through the computational domain.
However, most geometries consists of a small 2D surface in 3D space and the default Chombo approach wastes a lot of time looking for cut-cells where they don't exist. 

To achieve scalable geometry generation, we have changed how `Chombo` generates the geometry generation on the various levels.
Our new approach first generates a map on a *coarse* level which is specified by the user.
On the specified level the domain is broken up into equal-sized chunks and cut-cell boxes are located.
Uncut and cut boxes are load balanced among the various ranks.
We then proceed towards the next finer level where the cut-cell boxes are identified by a refinement of the box distribution on the previous level.
Boxes that resulted from a refinement of the coarse level cut boxes are again broken up into equal-sized chunks, whereas the uncut boxes are not.
This is again followed by load-balancing of the cut boxes, and this process is repeated recursively down to the finest AMR level.
In essence, the geometry generation is load balanced based on where the cut cells are going to be.
For the user, he will be able to switch between the `Chombo` and ``PlasmaC`` approaches to geometry generation load balancing by flipping a flag in an input script.
