.. _Chap:DischargeInceptionModel:

Discharge inception model
*************************

Overview
========

The discharge inception model computes the inception voltage and probability of discharge inception for arbitrary geometries and voltage forms.
For estimating the streamer inception, the module resolves the electron avalanche integral

.. math::
   :label: K_integral
	   
   K\left(\mathbf{x}\right) = \int_{\mathbf{x}}^{\text{until }\alpha_{\text{eff}} \leq 0} \alpha_{\text{eff}}(E,\mathbf{x}^\prime)\text{d}l

where :math:`E = |\mathbf{E}|` and where :math:`\alpha_{\text{eff}}(E,\mathbf{x}) = \alpha(E,\mathbf{x}) - \eta(E,\mathbf{x})` is the effective ionization coefficient.
The integration runs along electric field lines.

The discharge inception model solves for :math:`K = K\left(\mathbf{x}\right)` for all :math:`\mathbf{x}`, i.e., the inception integral is evaluated for all starting positions of the first electron.
This differs from the conventional approach where the user will first extract electric field lines for post-processing.
Note that :math:`K\left(\mathbf{x}\right) = K(\mathbf{x}; U)` where :math:`U` is the applied voltage.

In addition to the above, the user can specify a critical threshold value for :math:`K_c` which is used for computing

* The *critical volume* :math:`V_c = \int_{K\geq K_c} \textrm{d}V`.
* The inception voltage :math:`U_c`.
* The probability of having the first electron in the critical volume, :math:`dP(t,t+\Delta t)`.

The discharge inception model is implemented through the following files:

* ``DischargeInceptionStepper``, which implements :ref:`Chap:TimeStepper`.
* ``DischargeInceptionTagger``, which implements :ref:`Chap:CellTagger` and flags cells for refinement and coarsening.
* ``DischargeInceptionSpecies``, which implements :ref:`Chap:CdrSpecies` for the negative ion transport description.

.. tip::

   The source code for the discharge inception model is given in :file:`$DISCHARGE_HOME/Physics/DischargeInception`.
   See `DischargeInceptionStepper <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhysics_1_1DischargeInception_1_1DischargeInceptionStepper.html>`_ for the C++ API for this time stepper.

Townsend criterion
------------------

One may also solve for the Townsend inception criterion, which is formulated as follows:

.. math::

   T\left(\mathbf{x}\right) = \gamma\left(\exp\left[K\left(\mathbf{x}\right)\right]-1\right) \geq 1.

The interpretation of this criterion is that each starting electron produces :math:`\exp\left[K\left(\mathbf{x}\right)\right]-1` secondary electron-ion pairs.
The residual ions will drift towards cathode surfaces and generate secondary ionization with a user-supplied efficiency :math:`\gamma=\gamma\left(E,\mathbf{x}\right)`.

The discharge inception model can be run in two modes:

* A stationary mode, where one only calculations :math:`K\left(\mathbf{x}\right)` for a range of voltages, see :ref:`Chap:StationaryMode`.
* In transient mode, where :math:`K\left(\mathbf{x}\right)` is computed dynamically according to a user-supplied voltage shape.
  This mode can also be used to evaluate the inception probability for a given voltage curve, see :ref:`Chap:TransientMode`.

Implementation
--------------

The discharge inception model is implemented in :file:`$DISCHARGE_HOME/Physics/DischargeInception` as

.. literalinclude:: ../../../../Physics/DischargeInception/CD_DischargeInceptionStepper.H
   :language: c++
   :lines: 72-79
   :dedent: 4
	   
The template template parameters indicate which types of solvers to used within the compound algorithm.
The template parameters indicate the following: 

* ``typename P`` is the tracer particle solver (see :ref:`Chap:TracerParticleSolver`) for reconstructing the inception integral.
* ``typename F`` is the field solver (see :ref:`Chap:FieldSolver`) that is used when computing the electric field.
* ``typename C`` is the convection diffusion reaction solver (see :ref:`Chap:CdrSolver`) to use when resolving negative ion transport.

.. _Chap:StationaryMode:

Stationary mode
===============

The stationary mode resolves :eq:`K_integral` for a range of voltages, including both positive and negative polarities.
This procedure occurs through the following steps:

#. Solve for the background field with an applied voltage :math:`U = 1` on all live electrodes.
#. Obtain the critical field by obtaining the root :math:`\alpha_{\text{eff}}(E) = 0`.
   This step is performed using Brent's method.
#. Scale the voltage to the lowest voltage that exceeds the critical field.
#. Solve for :math:`K(\mathbf{x})` from the lowest voltage that exceeds the critical field until a user-specified threshold.
   These calculations are performed for both polarities.    

On output, the user is provided with the following:

* :math:`K` values for positive and negative voltages, as well as the Townsend criterion.
* Critical volumes and surfaces.
* Inception voltage, both for streamer inception and Townsend inception, for both polarities. 
* The *critical position*, i.e., the position corresponding to :math:`K(\mathbf{x})`.

.. _Chap:TransientMode:

Transient mode
==============

In transient mode the user can apply a voltage curve :math:`U = U(t)` reconstruct the :math:`K` value at each time step, and recompute the critical volume so that we obtain

.. math::

   K &= K(\mathbf{x}, t) \\
   T &= T(\mathbf{x}, t) \\
   V_c &= V_c(t).

The rationale for computing the above quantities is that one may describe the availability of negative ions within the critical volumes by augmenting the model with an ion transport description.
Furthermore, if one also has a description of the electron detachment rate, one may then evaluate the probability that an electron detaches from a negative ion inside of the critical region. 

Negative ion transport
----------------------

We also assume that ions move as drifting Brownian walkers in the electric field (see :ref:`Chap:ItoDiffusion`).
This can be written in the fluctuating hydrodynamics limit as an evolution equation for the ion density

.. math::
   
   \frac{\partial n_-}{\partial t} = -\nabla\cdot\left(\mathbf{v} n_-\right) + \nabla\cdot\left(D\nabla  n_-\right) + \sqrt{2D n_-}\mathbf{Z},

where :math:`\mathbf{Z}` represents uncorrelated Gaussian white noise.
Note that the above equation is a mere rewrite of the Ito process for a collection of particles; it is not really useful per se since it is a tautology for the original Ito process. 

However, we are interested in the average ion distribution over many experiments, so by taking the ensemble average we obtain a regular advection-diffusion equation for the evolution of the negative ion distribution :math:`\langle n_-\rangle`:

.. math::
   
   \frac{\partial \langle n_-\rangle}{\partial t} = -\nabla\cdot\left(\mathbf{v} \langle n_-\rangle\right) + \nabla\cdot\left(D\nabla \langle n_-\rangle\right).

This equation is sensible only when :math:`\langle n_-\rangle` is interpreted as an ion density distribution (over many identical experiments).

Inception probability
---------------------

The probability of discharge inception in a time interval :math:`[t,t+\text{d}t]` is given by

.. math::
   
   \text{d}P(t) = \left[1-P\left(t\right)\right]\lambda(t) \text{d}t,

where :math:`\lambda(t)` is a placeholder for the electron generation rate within the critical volume.
The exact solution for :math:`P(t)` is

.. math::

   P(t) = 1 - \exp\left(-\int_{-\infty}^t \lambda\left(t^\prime\right)\text{d}t^\prime\right).

An expression for the electron generation rate can be given by

.. math::

    \lambda(t) = \int_{V_c(t)}\left\langle\frac{\partial n_e}{\partial t}\right\rangle\left(1 - \frac{\eta}{\alpha}\right)\text{d}V + \int_{A_c(t)}\frac{\left\langle j_e\right\rangle}{e}\left(1-\frac{\eta}{\alpha}\right)\text{d} A,

Inserting the expression for :math:`\lambda` and integrating for :math:`P(t)` yields

.. math::
   :label: DischargeInceptionProbability
	   
   P(t) = 1  - \exp\left[-\int_0^t\left(\int_{V_c(t^\prime)}\left\langle\frac{dn_{\text{e}}}{dt^\prime}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V + \int_{A_c(t^\prime)}\frac{j_e}{q_{\text{e}}}\left(1-\frac{\eta}{\alpha}\right) \text{d}A\right)\text{d}t^\prime\right].

Here, :math:`\left\langle\frac{d n_{\text{e}}}{dt}\right\rangle` is the electron production rate from both background ionization and electron detachment, i.e.

.. math::

   \left\langle\frac{d n_{\text{e}}}{dt}\right\rangle = S_{\text{bg}} + k_d \left\langle n_-\right\rangle,

where :math:`S_{\text{bg}}` is the background ionization rate set by the user, :math:`k_d` is the negative ion detachment rate, and :math:`\left\langle n_-\right\rangle` is the negative ion distribution.
The second integral is due to electron emission from the cathode and into the critical volume.
Note that, internally, we always ensure that :math:`j_{\text{e}} dA` evaluates to zero on anode surfaces.

Statistical time lag
--------------------

We also compute the probability of a first electron appearing in the time interval :math:`[t, t+\Delta t]`, given by

.. math::
   :label: DischargeInceptionProbability2
   
   \Delta P(t, t+\Delta t) = \left[1-P(t)\right] \left(\int_{V_c(t^\prime)}\left\langle\frac{dn_{\text{e}}}{dt^\prime}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V + \int_{A_c(t^\prime)}\frac{j_e}{q_{\text{e}}}\left(1-\frac{\eta}{\alpha}\right) \text{d}A\right)\Delta t

When running in transient mode the user must set the voltage curve, and pay particular caution to setting the initial ion density, mobility, and detachment rates.

The statistical time lag, or average waiting time for the first electron, is available from the computed data, and is given by integrating :math:`t \text{d}P`, which yields

.. math::

   \tau = \frac{1}{P(t)}\int_0^\infty t\left[1-P(t)\right]\lambda(t)\text{d}t.

Other derived values (such as the standard deviation of the waiting time) is also available, and can be calculated from the :math:`P(t)` and :math:`\lambda(t)` similar to the procedure above.
Numerically, this is calculated using the trapezoidal rule. 

.. _Chap:DischargeInceptionInputData:

Input data
==========

The input to the discharge inception model are:

#. Streamer inception threshold.
#. Townsend ionization and attachment coefficient.
#. Background ionization rate (e.g., from cosmic radiation).
#. Electron detachment rate from negative ions.
#. Negative ion mobility and diffusion coefficients. 
#. Initial negative ion density.
#. Secondary emission coefficients.
#. Voltage curve (for transient simulations).
#. Space and surface charge if resolving a space-charge influenced field.

The input data to the discharge inception model is done by passing in C++-functions to the class.
These functions are mainly in the forms of an ``std::function``, so the user can pass in function pointers or lambdas for configuring the behavior of the model.
The relevant user API for setting the above variables are listed below.

.. literalinclude:: ../../../../Physics/DischargeInception/CD_DischargeInceptionStepper.H
   :language: c++
   :lines: 291-386
   :dedent: 4

.. tip::

   It is relatively simple to use tabulated data input for setting the input (e.g., the transport coefficients).
   See :ref:`Chap:LookupTable`.

Algorithms
==========

The discharge inception model uses a combination of electrostatic field solves, Particle-In-Cell, and fluid advection for resolving the necessary dynamics.
The various algorithms involved are discussed below.

Field solve
-----------

Since the background field scales linearly with applied voltage, we require only a single field solve at the beginning of the simulation.
This field solve is done with an applied voltage of :math:`U = 1\,\text{V}` and the electric field is then later scaled by the appropriate voltage.

Inception integral
------------------

We use a Particle-In-Cell method for computing the inception integral :math:`K\left(\mathbf{x}\right)` for an arbitrary electron starting position.
All grid cells where :math:`\alpha_{\textrm{eff}} > 0` are seeded with one particle on the cell centroid and the particles are then tracked through the grid.
The particles move a user-specified distance along field lines :math:`\mathbf{E}`,  using first or second order integration. 
If a particle leaves through a boundary (EB or domain boundary), or enters a region :math:`\alpha_{\text{eff}} \leq 0`, the integration for that particle is stopped.
Once all the particle integrations have halted, we rewind the particles back to their starting position and deposit their weight on the mesh, which provides us with :math:`K = K\left(\mathbf{x}\right)`.

.. note::

   When tracking positive ions for evaluation of the Townsend criterion, the same algorithms are used.

Euler
_____

For the Euler rule the particle weight for a particle :math:`p` the update rule is

.. math::

   \mathbf{X}_p^{k+1} = \mathbf{X}_p^k - \mathbf{\hat{E}}\left(\mathbf{X}_p^k\right)\Delta x
   
   w_p^{k+1} = w_p^k + \alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{X}_p^k\right)\right|,\mathbf{X}_p^k\right)\Delta x,

where :math:`\Delta x` is an integration length.

Trapezoidal
___________

With the trapezoidal rule we first update

.. math::

   \mathbf{X}_p^\prime = \mathbf{X}_p^k - \mathbf{\hat{E}}\left(\mathbf{X}_p^k\right)\Delta x,

which is followed by

.. math::

      \mathbf{X}_p^{k+1} = \mathbf{X}_p^k + \frac{\Delta x}{2}\left[\mathbf{\hat{E}}\left(\mathbf{X}_p^k\right) + \mathbf{\hat{E}}\left(\mathbf{X}_p^\prime\right)\right].

      w_p^{k+1} = w_p^k + \frac{\Delta x}{2}\left[\alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{X}_p^k\right)\right|,\mathbf{X}_p^k\right) + \alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{X}_p^\prime\right)\right|,\mathbf{X}_p^\prime\right)\right]

.. _Chap:DischargeInceptionStepSizeSelection:

Step size selection
___________________

The permitted tracer particle step size is controlled by user-specified maximum and minimum space steps:

.. math::

   &\Delta_{\textrm{phys,min}}: \textrm{Minimum physical step size}  \\
   &\Delta_{\textrm{phys,max}}: \textrm{Maximum physical step size}  \\
   &\Delta_{\textrm{grid,min}}: \textrm{Minimum grid-cell step size}  \\
   &\Delta_{\textrm{grid,max}}: \textrm{Maximum grid-cell step size}  \\
   &\Delta_{\alpha}:            \textrm{Avalanche length}  \\
   &\Delta_{\nabla\alpha}:      \textrm{Rate-of-change of ionization coefficient}.

The particle integration step size is then selected according to the following heuristic:

.. math::

   \Delta X &= \min\left(\Delta_\alpha\frac{1}{\overline{\alpha}}, \Delta_{\nabla\alpha}\frac{\overline{\alpha}}{\left|\nabla\overline{\alpha}\right|}\right) \\
   \Delta X &= \min\left(\Delta X, \Delta_{\textrm{phys,max}}\right) \\
   \Delta X &= \max\left(\Delta X, \Delta_{\textrm{phys,min}}\right) \\
   \Delta X &= \min\left(\Delta X, \Delta_{\textrm{grid,max}}\Delta x\right) \\
   \Delta X &= \max\left(\Delta X, \Delta_{\textrm{grid,min}}\Delta x\right).

These parameters are implemented through the following input options:

.. code-block:: text

   # Particle integration controls
   DischargeInceptionStepper.min_phys_dx      = 1.E-10		## Minimum permitted physical step size
   DischargeInceptionStepper.max_phys_dx      = 1.E99		## Maximum permitted physical step size
   DischargeInceptionStepper.min_grid_dx      = 0.5		## Minimum permitted grid step size
   DischargeInceptionStepper.max_grid_dx      = 5.0		## Maximum permitted grid step size
   DischargeInceptionStepper.alpha_dx         = 5.0		## Step size relative to avalanche length
   DischargeInceptionStepper.grad_alpha_dx    = 0.1		## Maximum step size relative to alpha/grad(alpha)
   DischargeInceptionStepper.townsend_grid_dx = 2.0		## Space step to use for Townsend tracking		

Note that the input variable ``townsend_grid_dx`` determines the spatial step (relative to the grid resolution :math:`\Delta x`) when tracking ions for the Townsend region reconstruction.


Critical volume
---------------

The critical volume is computed as

.. math::

   V_c = \int_{K\left(\mathbf{x}\right) \geq K_c \cup \gamma\exp\left[K\left(\mathbf{x}\right)\right] \ge 1} \text{d}V.

Note that the critical volume is both voltage and polarity dependent.

Critical surface
----------------

The critical surface is computed as

.. math::

   A_c = \int_{K\left(\mathbf{x}\right) \geq K_c \cup \gamma\left(\exp\left[K\left(\mathbf{x}\right)\right]-1\right) \ge 1} \text{d}A.

Note that the critical surface is both voltage and polarity dependent, and is non-zero only on cathode surfaces.

Inception voltage
-----------------

The inception voltage for starting a critical avalanche can be computed in the stationary solver mode (see :ref:`Chap:StationaryMode`).
This is done separately for the streamer and Townsend inception voltages.

Streamer inception
__________________

For streamer inception we use :math:`K\left(\mathbf{x}; U\right)` for a range of voltages :math:`U \in U_1, U_2, \ldots` and (linearly) interpolate between these values.
If two values of the :math:`K` integral bracket :math:`K_c`, i.e.

.. math::

   K_a &= K\left(\mathbf{x}; U_a\right) \leq K_c \\
   K_b &= K\left(\mathbf{x}; U_b\right) \geq K_c

then we can estimate the inception voltage for a starting electron at position :math:`\mathbf{x}` through linear interpolation as

.. math::

   U_{\text{inc, streamer}}\left(\mathbf{x}\right) = U_a + \frac{K_c - K_a}{K_b - K_a}\left(U_b - U_a\right)

Townsend inception
__________________

A similar method to the one used above is used for the Townsend criterion, using e.g. :math:`T\left(\mathbf{x}; U\right) = \gamma\left(\exp\left[K\left(\mathbf{x}; U\right)\right]-1\right)`, then if

.. math::
   
   T_a &= T\left(\mathbf{x}; U_a\right) \leq 1, \\
   T_b &= T\left(\mathbf{x}; U_b\right) \ge 1,

then we can estimate the inception voltage for a starting electron at position :math:`\mathbf{x}` through linear interpolation as

.. math::

   U_{\text{inc, Townsend}}\left(\mathbf{x}\right) = U_a + \frac{1 - T_a}{T_b - T_a}\left(U_b - U_a\right)


Minimum inception voltage
_________________________

The minium inception voltage is the minimum voltage required for starting a critical avalanche (or Townsend process) for an arbitrary starting electron.
For any position :math:`\mathbf{x}`, then

.. math::

   U_{\text{inc}}\left(\mathbf{x}\right) = \min\left[U_{\text{inc, streamer}}\left(\mathbf{x}\right), U_{\text{inc, Townsend}}\left(\mathbf{x}\right)\right]
   
From the above, this is simply

.. math::

   U_{\text{inc}}^{\text{min}} = \min_{\forall \mathbf{x}} \left[U_{\text{inc}}\left(\mathbf{x}\right)\right].

From the above we also determine

.. math::

   \mathbf{x}_{\text{inc}}^{\text{min}} \leftarrow \mathbf{x}\text{ that minimizes } U_{\text{inc}}\left(\mathbf{x}\right) \forall \mathbf{x},

which is the position of the first electron that enables a critical avalanche at the minimum inception voltage.

Inception probability
---------------------

The inception probability is given by :eq:`DischargeInceptionProbability` and is computed using straightforward numerical quadrature:

.. math::

   \int_{V_c}\left\langle\frac{dn_{\text{e}}}{dt}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V \approx \sum_{\mathbf{i}\in K_\mathbf{i} > K_c} \left(\left\langle\frac{dn_{\text{e}}}{dt}\right\rangle\right)_{\mathbf{i}}\left(1 - \frac{\eta_{\mathbf{i}}}{\alpha_{\mathbf{i}}}\right)\kappa_{\mathbf{i}}\Delta V_{\mathbf{i}},

and similarly for the surface integral.

.. important::

   The integration runs over *valid cells*, i.e. grid cells that are not covered by a finer grid.

Advection algorithm
-------------------

The advection algorithm for the negative ion distribution follows the time stepping algorithms described in the advection-diffusion model, see :ref:`Chap:AdvectionDiffusionModel`.

Adaptive mesh refinement
========================

The discharge inception model runs its own mesh refinement routine, which refines the mesh if

.. math::

   \alpha_{\text{eff}}\left(\left|\mathbf{E}\right|, \mathbf{x}\right)\Delta x > \lambda,

where :math:`\lambda` is a user-specified refinement criterion.
The user can control refinement buffers and criterion through the following input options (see :ref:`Chap:DischargeInceptionConfiguration`).

* ``DischargeInceptionTagger.buffer`` Adds a buffer region around tagged cells.
* ``DischargeInceptionTagger.max_voltage`` Maximum voltage that will be simulated.
* ``DischargeInceptionTagger.ref_alpha`` Sets the refinement criterion :math:`\lambda` as above.

.. tip::

   When using transient mode, it may be useful to simply set the maximum voltage to the peak voltage of the voltage curve.

   For stationary solves it might be difficult because the range of voltages is determined automatically during run-time.
   It may be beneficial to run the program twice, first using ``max_voltage = 1``, and then running the program again using the peak voltage from the output file.

   
.. _Chap:DischargeInceptionConfiguration:

Solver configuration
====================

The ``DischargeInceptionStepper`` class come with user-configurable input options which are given below.

.. literalinclude:: ../../../../Physics/DischargeInception/CD_DischargeInceptionStepper.options
   :language: text

In the above options, the user can select the integration algorithms, the mode, and where to place the output file (which contains, e.g., the values of the ionization integral). 
The user can also include the following data in the HDF5 output files, by setting the ``plt_vars`` configuration option:

* ``field``    - Potential, field, and charge distributions.
* ``K``        - Inception integral.
* ``T``        - Townsend criterion.  
* ``Uinc``     - Inception voltage.
* ``alpha``    - Effective ionization coefficient.
* ``eta``      - Eta coefficient.  
* ``bg_rate``  - Background ionization rate.
* ``emission`` - Field emission.
* ``poisson``  - Poisson solver.
* ``tracer``   - Tracer particle solver.
* ``cdr``      - CDR solver.
* ``ions``     - Ion solver.

.. important::

   The interface for setting the transport data (e.g., ionization coefficients) occurs via the C++ interface.              

.. _Chap:DischargeInceptionSetup:

Setting up a new problem
========================

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/DischargeInception` is the simplest way.
A full description is available in the ``README.md`` file contained in the folder:

.. literalinclude:: ../../../../Physics/DischargeInception/README.md
   :language: markdown
	      
To see available setup options, use

.. code-block:: bash

   python setup.py --help


Example programs
================

Example programs that use the discharge inception model are given in

High-voltage vessel
-------------------

* :file:`$DISCHARGE_HOME/Exec/Examples/DischargeInception/Vessel`.
  This program is set up in 2D (stationary) and 3D (transient) for discharge inception in atmospheric air.
  The input data is computed using BOLSIG+.

Electrode with surface roughness
--------------------------------

* :file:`$DISCHARGE_HOME/Exec/Examples/DischargeInception/ElectrodeRoughness`.
  This program is set up in 2D (stationary) and 3D (transient) for discharge inception on an irregular electrode surface. 
  We use SF6 transport data as input data, computed using BOLSIG+.

Electrode with surface roughness
--------------------------------

* :file:`$DISCHARGE_HOME/Exec/Examples/DischargeInception/ElectrodeRoughness`.
  This program is set up in 2D and 3D (stationary) mode, and includes the influence of the Townsend criterion.
