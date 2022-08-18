.. _Chap:StreamerInceptionModel:

Streamer inception
==================

The streamer inception model solves the avalanche integral

.. math::

   K\left(\mathbf{x}\right) = \int_{\mathbf{x}}^{\alpha \leq 0} \alpha(E)\text{d}l

where :math:`E = |\mathbf{E}|` and where :math:`\alpha(E)` is the effective ionization coefficient.
The coefficient runs along electric field lines.

``chombo-discharge`` solves for :math:`K = K\left(\mathbf{x}\right)` for all :math:`\mathbf{x}`.
I.e., the inception integral is evaluated over the entire volume.
This differs from the conventional approach where the user will first extract electric field lines for post-processing.

In addition to the above, the user can specify a critical threshold value for :math:`K_c` which is used for computing

* The critical volume :math:`V_c = \int_{K>K_c} \textrm{d}V`.
* The inception voltage :math:`U_c`.
* The drift times :math:`t_d`.

Solvers
-------

The streamer inception model uses

* The :ref:`Chap:FieldSolver` for computing the electric field.
* The tracer particle module (:ref:`Chap:TracerParticleSolver`) for reconstructing the inception integral.

Implementation
--------------

``chombo-discharge`` uses a Particle-In-Cell method to solve the inception integral. A particle is placed within each cell and integrated along the electric field lines until the particle exits the domain or the effective ionization coefficient :math:`\alpha(E)` becomes negative. The integration is executed with time step and integration algorithm specified from user input, the latter either Euler or trapezoidal integration. 


Effective ionization coefficient
----------------------

To set the effective ionization coefficient, use the member function

.. code-block:: c++

   setAlpha(const std::function<Real(const Real E)>&) noexcept;

For example:

.. code-block::

   StreamerInceptionStepper<> inceptionStepper;

   const auto alpha = [](const Real E) -> Real {
      return 1/E;
   }

   inceptionStepper.setAlpha(alpha);

Inception algorithm
----------------------

``StreamerInceptionStepper.inception_alg`` sets the inception algorithm parameters. The first
input is the integration algorithm, either Euler ("euler") or trapezoidal (trapz) integration.
The second input is the step algorithm, which decides whether the integration steps are relative
("dx") or fixed ("fixed") compared to the grid resolution. 
The third input is the integration step size.

For example:

.. code-block:: bash

		StreamerInceptionStepper.inception_alg = trapz fixed 0.2

Print report
-------------

Use ``StreamerInceptionStepper.print_report`` to print out values at the end of the simulation.
For example:

.. code-block:: bash

   StreamerInceptionStepper.print_report = true

Plot variables
---------------

``StreamerInceptionStepper.plt_vars`` sets which variables are plotted in the simulation.
The options are:

.. code-block:: bash

   poisson   # Electric field
   tracer    # particles
   neg_ions  #
   K         # Inception integral
   Uinc      # Inception voltage
   bg_rate   # Background ionization rate
   emission  # Field emission
   alpha     # Effective ionization coefficient
   eta       # Eta coefficient


* "poisson" - electric field
* "tracer" - particles
* "neg_ions"
* "K" - inception integral
* "Uinc" - inception voltage
* "bg_rate" - background ionization
* "emission" - field emission
* "alpha" - effective ionization coefficient
* "eta" - eta coefficient

For example:

.. code-block:: bash

		StreamerInceptionStepper.plt_vars = poisson neg_ions K emission Uinc
		
Static mode
------------
   
Voltage levels
^^^^^^^^^^^^^^^

By default, the streamer inception time stepper will read voltage levels from the input script.
These are in the format

.. code-block:: bash

   StreamerInceptionStepper.voltage_lo    = 1.0   # Low voltage multiplier
   StreamerInceptionStepper.voltage_hi    = 10.0  # Highest voltage multiplier
   StreamerInceptionStepper.voltage_steps = 3     # Number of voltage steps

Here, ``voltage_lo`` is the lowest voltage that we solve for, while ``voltage_hi`` is the highest voltage we solve for.

Inception threshold
^^^^^^^^^^^^^^^^^^^^

Use ``StreamerInceptionStepper.K_inception`` for setting the inception threshold.
For example:

.. code-block:: bash

   StreamerInceptionStepper.K_inception   = 18

Dynamic mode
-------------

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

* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/ElectrodeRoughness`.
* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/Armadillo`.

Caveats
-------

The model is intended to be used with a nearest-grid-point deposition scheme (which is also volume-weighted).
When running the model, ensure that the the :ref:`Chap:TracerParticleSolver` flags are set as follows:

.. code-block:: bash

   TracerParticleSolver.deposition   = ngp 
   TracerParticleSolver.volume_scale = true
