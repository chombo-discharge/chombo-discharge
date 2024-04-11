.. _Chap:Overview:

Overview
========

History
-------

``chombo-discharge`` is aimed at solving discharge problems.
It was originally developed at SINTEF Energy Research between 2015-2018, and aimed at simulating discharges in high-voltage engineering.
Further development was started in 2021, where much of the code was redesigned for improved modularity and performance. 

Key functionalities
-------------------

``chombo-discharge`` uses a Cartesian embedded boundary (EB) grid formulation and adaptive mesh refinement (AMR) formalism where the grids frequently change and are adapted to the solution as simulations progress.
Key functionalities are provided in :numref:`Tab:KeyCapabilities`. 

.. important::

   ``chombo-discharge`` is **not** a black-box model for discharge applications.
   It is an evolving research code with occasionally expanded capabilities, API changes, and performance improvements.
   Although problems can be set up through our Python toos, users should nonetheless be willing to take time to understand how the code works.
   In particular, developers should invest some effort in understanding the data structures and ``Chombo`` basics (see :ref:`Chap:Basics`). 

.. _Tab:KeyCapabilities:
.. list-table:: Key capabilities.
   :widths: 20 50
   :header-rows: 1

   * - Capabilities
     - Supported?
   * - Grids
     - **Fundamentally Cartesian**.
   * - Parallelized?
     - **Yes**, using OpenMP, MPI, or MPI+OpenMP.
   * - Load balancing?
     - **Yes**, with support for individual particle and fluid load balancing. 
   * - Complex geometries?
     - **Yes**, using embedded boundaries (i.e., cut-cells). 
   * - Adaptive mesh refinement?
     - **Yes**, using patch-based refinement.
   * - Subcycling in time?
     - **No**, only global time step methods.
   * - Computational particles?
     - **Yes.**
   * - Linear solvers?
     - **Yes**, using geometric multigrid in complex geometries. 
   * - Time discretizations?
     - **Mostly explicit.**
   * - Parallel IO?
     - **Yes**, using HDF5.
   * - Checkpoint-restart?
     - **Yes**, for both fluid and particles. 

An early version of ``chombo-discharge`` used sub-cycling in time, but this has now been replaced with global time stepping methods. 
That is, all the AMR levels are advanced using the same time step.       
``chombo-discharge`` has also incorporated many changes to the EB functionality supplied by ``Chombo``.
This includes much faster grid generation, support for polygon surfaces, and many performance optimizations (in particular to the EB formulation).

``chombo-discharge`` supports both fluid and particle methods, and can use multiply parallel distributed grids (see :ref:`Chap:Realm`) for individually load balancing particle and fluid kernels. 
Although many abstractions are in place so that user can describe a new set of physics, or write entirely new solvers into ``chombo-discharge`` and still use the embedded boundary formalism, ``chombo-discharge`` also provides several physics modules for describing various types of problems.

Organization
------------

The ``chombo-discharge`` source files are organized as follows:

.. list-table:: Code organization.
   :widths: 10 50
   :header-rows: 1

   * - Folder
     - Explanation
   * - :file:`Source`
     -  Source files for the AMR core, solvers, and various utilities.
   * - :file:`Physics`
     - Various implementations that can run the ``chombo-discharge`` source code.
   * - :file:`Geometries`
     - Various geometries.
   * - :file:`Submodules`
     - Git submodule dependencies.
   * - :file:`Exec`
     - Various executable applications. 
