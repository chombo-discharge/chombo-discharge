.. _Chap:StreamerInceptionModel:

Streamer inception
==================

The streamer inception model solves the integral

.. math::

   K\left(\mathbf{x}\right) = \int_{\mathbf{x}}^{\alpha \leq 0} \alpha(E)\text{d}l

where :math:`E = |\mathbf{E}|` and where :math:`\alpha(E)` is the effective ionization coefficient.
The coefficient runs along electric field lines.

``chombo-discharge`` solves for :math:`K = K\left(\mathbf{x}\right)` for all :math:`\mathbf{x}`.
I.e., the inception integral is evaluated over the entire volume.
This differs from the conventional approach where the user will first extract electric field lines for post-processing.

In addition to the above, the user can specified a critical threshold value for :math:`K_c` which is used for computing

* The critical volume :math:`V_c = \int_{K>K_c} \textrm{d}V`.
* The inception voltage :math:`U_c`.

Solvers
-------

The streamer inception model uses

* The :ref:`Chap:FieldSolver` for computing the electric field.
* The tracer particle module (:ref:`Chap:TracerParticleSolver`) for reconstructing the inception integral.

Implementation
--------------

Ionization coefficient
----------------------

To set the ionization coefficient, use the member function

.. code-block:: c++

   setAlpha(const std::function<Real(const Real E)>&) noexcept;

For example:

.. code-block::

   StreamerInceptionStepper<> inceptionStepper;

   const auto alpha = [](const Real E) -> Real {
      return 1/E;
   }

   inceptionStepper.setAlpha(alpha);

Voltage levels
--------------

By default, the streamer inception time stepper will read voltage levels from the input script.
These are in the format

.. code-block:: bash

   StreamerInceptionStepper.voltage_lo    = 1.0   # Low voltage multiplier
   StreamerInceptionStepper.voltage_hi    = 10.0  # Highest voltage multiplier
   StreamerInceptionStepper.voltage_steps = 3     # Number of voltage steps

Here, ``voltage_lo`` is the lowest voltage that we solve for, while ``voltage_hi`` is the highest voltage we solve for.

Inception threshold
-------------------

Use ``StreamerInceptionStepper.K_inception`` for setting the inception threshold.
For example:

.. code-block:: bash

   StreamerInceptionStepper.K_inception   = 18		

Background ionization rate
--------------------------

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/StreamerInception` is the simplest way.
To see available setup options, run

.. code-block:: bash

   ./setup.py --help

For example, to set up a new problem in :file:`$DISCHARGE_HOME/MyApplications/MyStreamerInception` for a cylinder geometry, run

.. code-block:: bash

   ./setup.py -base_dir=MyApplications -app_name=MyStreamerInception -geometry=Cylinder

This will set up a new problem in a cylinder geometry (defined in :file:`Geometries/Cylinder`).

Example programs
----------------

Example programs that use the streamer inception model are given in

* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/Cylinder`.
* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/Armadillo`.
