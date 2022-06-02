.. _Chap:ElectrostaticsModel:

Electrostatics model
====================

The electrostatics model solves

.. math::

   \nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0}

subject to the constraints and boundary conditions given in :ref:`Chap:FieldSolver`.

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/Electrostatics` is the simplest way.
To see available setup options, run

.. code-block:: bash

   ./setup.py --help

For example, to set up a new problem in :file:`$DISCHARGE_HOME/MyApplications/MyElectrostaticsProblem` for a cylinder geometry, run

.. code-block:: bash

   ./setup.py -base_dir=MyApplications -app_name=MyElectrostaticsProblem -geometry=Cylinder

This will set up a new problem in a cylinder geometry (defined in :file:`Geometries/Cylinder`).

Example programs
----------------

Some example programs for this module are given in

* :file:`$DISCHARGE_HOME/Exec/Examples/Electrostatics/Armadillo`
* :file:`$DISCHARGE_HOME/Exec/Examples/Electrostatics/Mechshaft`
* :file:`$DISCHARGE_HOME/Exec/Examples/Electrostatics/ProfiledSurface`  


Verification
------------

Spatial convergence tests for this module (and consequently the underlying numerical discretization) are given in

* :file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C1`
* :file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C2`
* :file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C3`

All tests operate by computing solutions on grids with resolutions :math:`\Delta x` and :math:`\Delta x/2` and then computing the solution error using the approach outlined in :ref:`Chap:SpatialConvergence`. 


C1: Coaxial cable (2D)
______________________

C2: Profiled surface (2D)
_________________________

C3: Mechanical shaft (3D)
_________________________

