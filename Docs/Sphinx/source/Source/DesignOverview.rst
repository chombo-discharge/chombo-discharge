.. _Chap:DesignOverview:

Overview
========

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

The major design components in ``chombo-discharge`` are:

#. :ref:`Chap:Driver` for running simulations.
#. :ref:`Chap:AmrMesh` for encapsulating (almost) all AMR and EB functionality in a common core class.
#. :ref:`Chap:TimeStepper` for integrating the equations of motion.
#. :ref:`Chap:ComputationalGeometry` for representing computational geometries (such as electrodes and dielectrics).
#. :ref:`Chap:CellTagger` for flagging cells for refinement and coarsening.
#. :ref:`Chap:GeoCoarsener` for manually removing refinement flags along the EB surface.

.. caution::

   In the current version of ``chombo-discharge``, :ref:`Chap:GeoCoarsener` has become redundant (due to improvements in the algorithms).
   It will be removed in future versions. 

