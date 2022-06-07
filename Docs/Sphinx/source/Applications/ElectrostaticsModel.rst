.. _Chap:ElectrostaticsModel:

Electrostatics model
====================

The electrostatics model solves

.. math::

   \nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0}

subject to the constraints and boundary conditions given in :ref:`Chap:FieldSolver`.
The implementation is by a class

.. code-block:: c++

   class FieldStepper : public TimeStepper

and is found in :file:`$DISCHARGE_HOME/Physics/Electrostatics`.
The C++ API is found `here <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhysics_1_1Electrostatics_1_1FieldStepper.html>`_.

Solvers
-------

This module uses the following solvers:

#. :ref:`Chap:FieldSolver`.

Time stepping
-------------

``FieldStepper`` is a class with only stationary solves.
This is enforced by making the ``TimeStepper::advance`` method throw an error, while the field solve itself is done in the initialization procedure in ``TimeStepper::postInitialize``.

.. important::
   To run the solver, one must set ``Driver.max_steps = 0``.

Setting the space charge
------------------------

Default behavior
________________

By default, the initial space charge for this problem is given by Gaussian

.. math::

   \rho\left(\mathbf{x}\right) = \rho_0\exp\left(-\frac{\left|\mathbf{x}-\mathbf{x}_0\right|^2}{2R^2}\right),

where :math:`\rho_0` is an amplitude, :math:`\mathbf{x}_0` is the center and :math:`R` is the Gaussian radius.
These are set by the input options

.. code-block:: none

   FieldStepper.init_rho     = 1.0  
   FieldStepper.rho_center   = 0 0 0
   FieldStepper.rho_radius   = 1.0


Custom value
____________

For a more general way of specifying the space charge, ``FieldStepper`` has a public member function

.. code-block:: c++

   void setRho(const std::function<Real(const RealVect& a_pos)>& a_rho, const phase::which_phase a_phase) noexcept;

Setting the surface charge
--------------------------

By default, the initial surface charge is set to a constant.
This constant is given by

.. code-block:: none

   FieldStepper.init_sigma = 1.0

Custom value
____________

For a more general way of specifying the surface charge, ``FieldStepper`` has a public member function

.. code-block:: c++

   void setSigma(const std::function<Real(const RealVect& a_pos)>& a_sigma) noexcept;   

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


C1: Coaxial cable
_________________

:file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C1` is a spatial convergence test in a coaxial cable geometry. 
Figure :numref:`Fig:ElectrostaticsC1_1` shows the field distribution.
Note that there is a dielectric embedded between the two cylindrical shells. 

.. _Fig:ElectrostaticsC1_1:
.. figure:: /_static/figures/ElectrostaticsC1_1.png
   :width: 40%
   :align: center

   Field distribution for a coaxial cable geometry on a :math:`512^2` grid. 

The computed convergence rates are given in :numref:`Fig:ElectrostaticsC1_2`.
We find second order convergence in all three norms. 

.. _Fig:ElectrostaticsC1_2:
.. figure:: /_static/figures/ElectrostaticsC1_2.png
   :width: 45%
   :align: center

   Spatial convergence rates for two-dimensional coaxial cable geometry. 


C2: Profiled surface
____________________

:file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C2` is a 2D spatial convergence test for an electrode and a dielectric slab with surface profiles.
Figure :numref:`Fig:ElectrostaticsC2_1` shows the field distribution. 

.. _Fig:ElectrostaticsC2_1:
.. figure:: /_static/figures/ElectrostaticsC2_1.png
   :width: 45%
   :align: center

   Field distribution for a profiled surface geometry on a :math:`2048^2` grid. 

The computed convergence rates are given in :numref:`Fig:ElectrostaticsC2_2`.
We find second order convergence in all three norms. 

.. _Fig:ElectrostaticsC2_2:
.. figure:: /_static/figures/ElectrostaticsC2_2.png
   :width: 45%
   :align: center

   Spatial convergence rates for two-dimensional dielectric surface profile. 

C3: Dielectric shaft
____________________

:file:`$DISCHARGE_HOME/Exec/Convergence/Electrostatics/C3` is a 3D spatial convergence test for a dielectric shaft perpendicular to the background field. 
Figure :numref:`Fig:ElectrostaticsC3_1` shows the field distribution for a :math:`256^3` grid. 

.. _Fig:ElectrostaticsC3_1:
.. figure:: /_static/figures/ElectrostaticsC3_1.png
   :width: 25%
   :align: center

   Field distribution for a profiled surface geometry on a :math:`256^3` grid. 

The computed convergence rates are given in :numref:`Fig:ElectrostaticsC3_2`.
We find second order convergence in :math:`L_1` and :math:`L_2` on all grids, and find second order convergence in the max-norm on sufficiently fine grids.
On coarser grids, the reduced convergence rate in the max-norm is probably due to under-resolution of the geometry. 

.. _Fig:ElectrostaticsC3_2:
.. figure:: /_static/figures/ElectrostaticsC3_2.png
   :width: 45%
   :align: center

   Spatial convergence rates.

   
     

   
