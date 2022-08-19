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
* The inception probability :math:`dP(t,t+\Delta t) = [1-P(t)] \int_{V_{c}} \frac{dn_e}{dt}(1-\frac{\eta}{\alpha}}) dV \Delta t`.

Solvers
-------

The streamer inception model uses

* The :ref:`Chap:FieldSolver` for computing the electric field.
* The tracer particle module (:ref:`Chap:TracerParticleSolver`) for reconstructing the inception integral.

Implementation
--------------

``chombo-discharge`` uses a Particle-In-Cell method to solve the inception integral. A particle is placed within each cell in the grid and integrated along the electric field lines until the particle exits the domain or the effective ionization coefficient :math:`\alpha(E)` becomes negative. The running integration tracker is stored locally for each particle and is added to for each integration step. When the particle exits the domain or has a negative :math:`\alpha` it is flagged and its integration is finished. 
The integration is executed for both polarities (+/-) with time step and integration algorithm specified from user input, the latter either Euler or trapezoidal integration.

The Euler algorithm gives:

.. math::

   K += \alpha(E(x)) * dx,

while the trapezoidal algorithm gives:

.. math::

   K += 0.5 * dx * [\alpha(E(x)) + \alpha(E(x+dx))].

The critical volume is calculated by adding the volumes of the cells where :math:`K>=K_c`.

The inception voltage is solved by linear interpolation between the :math:`K` values for the different voltage levels. It is computed for every :math:`\mathbf{x}`, i.e. every cell in the grid, where :math:`min(K)<=K_c` and :math:`max(K)>K_c`. Inception voltages are only computed when there are at least two voltage levels, for obvious reasons.

The probability of inception, i.e. that an electron appears within the critical volume in an interval :math:`[t, t + \Delta t]`, is computed assuming the background ionization rate is only affected by field emission and detachment of electrons from negatively charged ions.

Effective ionization coefficient
---------------------------------

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
input is the integration algorithm, either Euler (``euler``) or trapezoidal (``trapz``) integration.
The second input is the step algorithm, which decides whether the integration steps are relative
(``dx``) or fixed (``fixed``) compared to the grid resolution. 
The third input is the integration step size.

For example:

.. code-block:: bash

		StreamerInceptionStepper.inception_alg = trapz fixed 0.2

Print report & output file
---------------------------

Use ``StreamerInceptionStepper.print_report`` to save the values for voltages, maximum :math:`K`, critical volume, and Rdot (time to appearance of electrons within the critical volume)  at the end of the simulation.
For example:

.. code-block:: bash

   StreamerInceptionStepper.print_report = true

The report is stored to the file specified by ``StreamerInceptionStepper.output_file``, for example:

.. code-block:: bash

   StreamerInceptionStepper.output_file = report.txt


Plot variables
---------------

``StreamerInceptionStepper.plt_vars`` sets which variables are plotted in the simulation.
The options are:

* ``poisson``  - Electric field
* ``tracer``   - Particles
* ``neg_ions`` - Negative ions
* ``K``        - Inception integral
* ``Uinc``     - Inception voltage
* ``bg_rate``  - Background ionization rate
* ``emission`` - Field emission
* ``alpha``    - Effective ionization coefficient
* ``eta``      - Eta coefficient

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

The background ionization rate is calculated assuming contributions from detachment of electrons from negative ions and field emission. 

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

The figure below shows an example of the avalanche integral solved for an |SF6| gas with an irregular electrode surface:

.. |SF6| replace:: SF\ :sub:`6`

.. image:: chombo-discharge/Exec/Examples/StreamerInception/ElectrodeRoughness/field.png
   :width: 400

  
Caveats
-------

The model is intended to be used with a nearest-grid-point deposition scheme (which is also volume-weighted).
When running the model, ensure that the the :ref:`Chap:TracerParticleSolver` flags are set as follows:

.. code-block:: bash

   TracerParticleSolver.deposition   = ngp 
   TracerParticleSolver.volume_scale = true
