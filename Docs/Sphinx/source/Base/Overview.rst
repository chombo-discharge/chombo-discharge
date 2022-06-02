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
     - **Yes**, using flat MPI.
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
See e.g. :ref:`Chap:ImplementedModels`.

       
Design approach
---------------

A fundamental design principle in ``chombo-discharge`` is the division between the AMR core, geometry, solvers, physics coupling, and user applications. 
As an example, the fundamental time integrator class ``TimeStepper`` in ``chombo-discharge`` is a just an abstraction, i.e. it only presents an API which application codes must use. 
Because of that, ``TimeStepper`` can be used for solving completely unrelated problems. 
We have, for example, implementations of ``TimeStepper`` for solving radiative transfer equations, advection-diffusion problems, electrostatic problems, or for plasma problems.

The division between computational concepts (e.g., AMR functionality and solvers) exists so that users will be able to solve problems across a range of geometries, add new solvers functionality, or write entirely new applications, without requiring deep changes to ``chombo-discharge``.
:numref:`Fig:Design` shows the basic design sketch of the ``chombo-discharge`` code.
To the right in this figure we have the AMR core functionality, which supplies the infrastructure for running the solvers. 
In general, solvers may share common features (such as elliptic discretizations) or be completely disjoint.
For this reason numerical solvers are asked to *register* AMR requirements.
For example, elliptic solvers need functionality for interpolating ghost cells over the refinement boundary, but pure particle solvers have no need for such functionality.
A consequence of this is that the numerical solvers are literally asked (during their instantiation) to register what type of AMR infrastructure they require. 
In return, the AMR core will allocate this infrastructure and make it available to solver, as illustrated in :numref:`Fig:Design`. 

.. _Fig:Design:
.. figure:: /_static/figures/Design.png
   :width: 70%
   :align: center

   Concept design sketch for ``chombo-discharge``. 

``chombo-discharge`` also uses *loosely coupled* solvers as a foundation for the code design, where a *solver* indicates a piece of code for solving an equation.
For example, solving the Laplace equation :math:`\nabla^2\Phi = 0` is encapsulated by one of the ``chombo-discharge`` solvers.
Some solvers in ``chombo-discharge`` have a null-implemented API, i.e. we have enforced a strict separation of the solver interface and the solver implementation.
This constraint exists because while new features may be added to a discretization, we do not want such changes to affect upstream application code.
An example of this is the ``FieldSolver``, which conceptualizes a numerical solver for solving for electrostatic field problems.
The ``FieldSolver`` is an API with no fundamental discretization -- it only contains high-level routines for understanding the type of solver being dealt with. 
Yet, it is the ``FieldSolver`` API which is used by application code (and not the implementation class).

All numerical solvers interact with a common AMR core that encapsulates functionality for running the solvers.
All solvers are also compatible with mesh refinement and complex geometries, but they can only run through *application codes*, i.e. *physics modules*. 
These modules encapsulate the time advancement of either individual or coupled solvers.
Solvers only interact with one another through these modules, and these modules usually advance the equations of motion using the method-of-lines. 
