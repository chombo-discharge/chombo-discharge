.. _Chap:StreamerInceptionModel:

Streamer inception model
========================

Overview
--------

The streamer inception model computes the inception voltage and probability of streamer inception for arbitrary geometries and voltage forms.

For estimating the streamer inception, the module solves the electron avalanche integral

.. math::

   K\left(\mathbf{x}\right) = \int_{\mathbf{x}}^{\text{until }\alpha_{\text{eff}} \leq 0} \alpha_{\text{eff}}(E,\mathbf{x}^\prime)\text{d}l

where :math:`E = |\mathbf{E}|` and where :math:`\alpha_{\text{eff}}(E,\mathbf{x}) = \alpha(E,\mathbf{x}) - \eta(E,\mathbf{x})` is the effective ionization coefficient.
The integration runs along electric field lines.

The streamer inception model solves for :math:`K = K\left(\mathbf{x}\right)` for all :math:`\mathbf{x}`, i.e., the inception integral is evaluated for all starting positions of the first electron.
This differs from the conventional approach where the user will first extract electric field lines for post-processing.

In addition to the above, the user can specify a critical threshold value for :math:`K_c` which is used for computing

* The critical volume :math:`V_c = \int_{K>K_c} \textrm{d}V`.
* The inception voltage :math:`U_c`.
* The probability of having the first electron in the critical volume, :math:`dP(t,t+\Delta t)`.

Note that since the :math:`K\left(\mathbf{x}\right) = K(\mathbf{x}; U)` where :math:`U` is the applied voltage, these values are computed for a user-specified range of voltages. 
This "range of voltages" can be a series of discrete values, or a voltage curve (e.g., lightning impulse).

In addition to this, one may also track positive ions and solve for the Townsend inception criterion, which is formulated as

.. math::

   T\left(\mathbf{x}\right) = \gamma \exp\left[K\left(\mathbf{x}\right)\right] \geq 1.

The interpretation of this criterion is that each starting electron grows into :math:`\exp\left[K\left(\mathbf{x}\right)\right]` electron-ion pairs.
The residual ions will drift towards cathode surfaces and generate secondary ionization with a user-supplied efficiency :math:`\gamma=\gamma\left(E,\mathbf{x}\right)`.

The streamer inception model can be run in two modes:

* :ref:`Chap:StationaryMode`.
* :ref:`Chap:TransientMode`.

These are discussed below.

Implementation
______________

The streamer inception model is implemented in :file:`$DISCHARGE_HOME/Physics/StreamerInception` as

.. code-block:: c++

   template <typename P, typename F, typename C>
   class StreamerInceptionStepper : public TimeStepper

The template complexity is due to flexibility in what solver implementations we use.
The following solvers are used within the model:

* The tracer particle solver (:ref:`Chap:TracerParticleSolver`) for reconstructing the inception integral.
* The :ref:`Chap:FieldSolver` for computing the electric field.
* The convection diffusion reaction module (:ref:`Chap:CdrSolver`) for modeling negative ion transport. 

.. _Chap:StationaryMode:

Stationary mode
_______________

The stationary mode solves the :math:`K` integral for a range of input voltages.
The computation is done for both polarities so that the user obtains:

* :math:`K` values for positive and negative voltages, as well as the Townsend criterion.
* Critical volumes :math:`V_c` for positive and negative voltages.
* Inception voltage, using a combined Townsend-streamer criterion.

.. _Chap:TransientMode:

Transient mode
______________

In the transient mode we apply a voltage curve :math:`U = U(t)` reconstruct the :math:`K` value at each time step and recompute the critical volume so that we obtain

.. math::

   K = K(t) \\
   T = T(t) \\
   V_c = V_c(t)

We also assume that ions move as drifting Brownian walkers in the electric field (see :ref:`Chap:ItoDiffusion`).
This can be written in the fluctuating hydrodynamics limit as an evolution equation for the ion distribution (not density!) as

.. math::
   
   \frac{\partial \langle n_-\rangle}{\partial t} = -\nabla\cdot\left(\mathbf{v} \langle n_-\rangle\right) + \nabla\cdot\left(D\nabla \langle n_-\rangle\right) + \sqrt{2D\langle n_-\rangle}\mathbf{Z},

where :math:`\mathbf{Z}` represents uncorrelated Gaussian white noise.
Note that the above equation is a mere rewrite of the Ito process for a collection of particles; it is not really useful per se since it is a tautology for the original Ito process. 

However, we are interested in the average ion distribution over many experiments, so by taking the ensemble average we obtain a regular advection-diffusion equation for the evolution of the negative ion distribution (note that we redefine :math:`\langle n_-\rangle` to be the ensemble average).

.. math::
   
   \frac{\partial \langle n_-\rangle}{\partial t} = -\nabla\cdot\left(\mathbf{v} \langle n_-\rangle\right) + \nabla\cdot\left(D\nabla \langle n_-\rangle\right).

This equation is sensible only when :math:`\langle n_-\rangle` is interpreted as an ion density distribution (over many identical experiments). 

The above quantities are then used for computing the probability of streamer inception in a time interval :math:`[t,t+\text{d}t]`, which is

.. math::
   \text{d}P(t) = \left[1-P\left(t\right)\right]\lambda(t) \text{d}t,

where :math:`\lambda(t)` is a placeholder for the electron generation rate, given by

.. math::

    \lambda(t) = \int_{V_c(t)}\left\langle\frac{\partial n_e}{\partial t}\right\rangle\left(1 - \frac{\eta}{\alpha}\right)\text{d}V + \int_{A_c(t)}\frac{\left\langle j_e\right\rangle}{e}\left(1-\frac{\eta}{\alpha}\right)\text{d} A,

Inserting the expression for :math:`\lambda` and integrating for :math:`P(t)` yields

.. math::
   :label: StreamerInceptionProbability
	   
   P(t) = 1  - \exp\left[-\int_0^t\left(\int_{V_c(t^\prime)}\left\langle\frac{dn_{\text{e}}}{dt^\prime}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V + \int_{A_c(t^\prime)}\frac{j_e}{q_{\text{e}}}\left(1-\frac{\eta}{\alpha}\right) \text{d}A\right)\text{d}t^\prime\right].

Here, :math:`\left\langle\frac{d n_{\text{e}}}{dt}\right\rangle` is the electron production rate from both background ionization and electron detachment, i.e.

.. math::

   \left\langle\frac{d n_{\text{e}}}{dt}\right\rangle = S_{\text{bg}} + k_d \left\langle n_-\right\rangle,

where :math:`S_{\text{bg}}` is the background ionization rate set by the user, :math:`k_d` is the negative ion detachment rate, and :math:`\left\langle n_-\right\rangle` is the negative ion distribution.
The second integral is due to electron emission from the cathode and into the critical volume.
Note that, internally, we always ensure that :math:`j_{\text{e}} dA` evaluates to zero on anode surfaces.

We also compute the probability of a first electron appearing in the time interval :math:`[t, t+\Delta t]`, given by

.. math::
   :label: StreamerInceptionProbability2
   
   \Delta P(t, t+\Delta t) = \left[1-P(t)\right] \left(\int_{V_c(t^\prime)}\left\langle\frac{dn_{\text{e}}}{dt^\prime}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V + \int_{A_c(t^\prime)}\frac{j_e}{q_{\text{e}}}\left(1-\frac{\eta}{\alpha}\right) \text{d}A\right)\Delta t

When running in transient mode the user must set the voltage curve (see :ref:`Chap:StreamerInceptionVoltageCurve`) and pay particular caution to setting the initial ion density, mobility, and detachment rates.

The statistical time lag, or average waiting time for the first electron, is available from the computed data, and is given by integrating :math:`t \text{d}P`, which yields

.. math::

   \tau = \int_0^\infty t\left[1-P(t)\right]\lambda(t)\text{d}t.

Other derived values (such as the standard deviation of the waiting time) is also available, and can be calculated from the :math:`P(t)` and :math`\lambda(t)` similar to the procedure above.
Numerically, this is calculated using the trapezoidal rule. 

.. _Chap:StreamerInceptionInputData:

Input data
----------

The input to the streamer inception model are:

#. Streamer inception threshold.
#. Townsend ionization coefficient.
#. Townsend attachment coefficients.
#. Background ionization rate (e.g., from cosmic radiation).
#. Electron detachment rate (from negative ions).
#. Negative ion mobility.
#. Negative ion diffusion coefficient.   
#. Initial negative ion density.
#. Secondary emission coefficients.
#. Voltage curve (for transient simulations).

The input data to the streamer inception model is mostly done by passing in C++-functions to the class.
These functions are mainly in the forms

.. code-block:: c++

   std::function<Real(const Real& E)>
   std::function<Real(const Real& E, const RealVect& x)>

The user can specify analytic fields or use tabulated data, and pass these in through a C++ lambda function.
An example of defining an analytic input function is

.. code-block:: c++

   auto alphaCoeff = [](const Real& E, const RealVect& x) -> void {
      return 1/E.
   };

Tabulated data (see :ref:`Chap:LookupTable1D`) can also be used as follows,

.. code-block:: c++
		
   LookupTable1D<2> tableData;
   
   auto alphaCoeff = [tableData](const Real& E, const RealVect& x) -> void {
      return tableData.getEntry<1>(E);
   };

.. note::

   The :math:`K` integral is determined only by the Townsend ionization and attachment coefficients.
   The Townsend criterion is then a derived value of :math:`K` and the secondary electron emission coefficient :math:`\gamma` .
   The remaining transport data is used for calculating the inception probability (appearance of a first electron in the critical volume). 
   

Inception threshold
___________________

Use in class input value ``StreamerInceptionStepper.K_inception`` for setting the inception threshold.

For example:

.. code-block:: text

   StreamerInceptionStepper.K_inception   = 12.0

Townsend ionization coefficient
_______________________________

To set the Townsend ionization coefficient, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setAlpha(const std::function<Real(const RealVect& E, const RealVect& x)>& a_alpha) noexcept;


Townsend attachment coefficient
_______________________________

To set the Townsend attachment coefficient, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setEta(const std::function<Real(const Real& E, const RealVect& x)>& a_eta) noexcept;
   

Negative ion mobility
_____________________

To set the negative ion mobility, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setIonMobility(const std::function<Real(const Real& E)>& a_mobility) noexcept;
   

Negative ion diffusion coefficient
__________________________________

To set the negative ion diffusion coefficient, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setIonDiffusion(const std::function<Real(const Real& E)>& a_diffCo) noexcept;   


Negative ion density
____________________

To set the negative ion density, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setIonDensity(const std::function<Real(const RealVect x)>& a_density) noexcept;

Secondary emission
__________________

To set the secondary emission efficiency at cathodes, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setSecondaryEmission(const std::function<Real(const Real& E, const RealVect& x)>& a_coeff) noexcept;

This efficiency is position-dependent so that the user can set different efficiencies for different materials (or different positions in a single material).

   
Background ionization rate
__________________________

The background ionization rate describes the appearance of a first electron from a background contribution, e.g. through cosmic radiation, decay of radioactive isotopes, etc.

To set the background ionization rate, use the member function

.. code-block:: c++

   StreamerInceptionStepper::setBackgroundRate(const std::function<Real(const Real& E, const RealVect& x)>& a_backgroundRate) noexcept;

Detachment rate
_______________

The detachment rate from negative describes the apperance of electrons through the equation

.. math::

   \left\langle\frac{dn_{\text{e}}}{dt}\right\rangle = k_d \left\langle n_-\right\rangle

where :math:`\left\langle n_-\right\rangle` is the negative ion density in units of :math:`m^{-3}` (or strictly speaking the negative ion probability density). 
This is used when calculating the inception probability, and the user sets the detachment rate :math:`k_d` through

.. code-block:: c++
		
   StreamerInceptionStepper::setDetachmentRate(const std::function<Real(const Real& E, const RealVect& x)>& a_backgroundRate) noexcept;

Field emission
______________

To set the field emission current, use the function

.. code-block:: c++

   StreamerInceptionStepper::setFieldEmission(const std::function<Real(const Real& E, const RealVect& x)>& a_currentDensity) noexcept;

This will set a field-dependent emission rate from cathodes given by the input function.
Note that, under the hood, the function indicates a general cathode emission current which can be the sum of several contributions (field emission, photoelectric effect etc.).

.. important::

   The input function should provide the surface current density :math:`j_e` (in units of :math:`\text{C}\cdot\text{m}^{-2}\cdot \text{s}^{-1}`).

Input voltages
______________

By default, the model will always read voltage levels from the input script.
These are in the format

.. code-block:: text

   StreamerInceptionStepper.voltage_lo    = 1.0   # Low voltage multiplier
   StreamerInceptionStepper.voltage_hi    = 10.0  # Highest voltage multiplier
   StreamerInceptionStepper.voltage_steps = 3     # Number of voltage steps



.. _Chap:StreamerInceptionVoltageCurve:

Voltage curve
_____________

To set the voltage curve, use the member function

.. code-block:: c++

   StreamerInceptionSteppersetVoltageCurve(const std::function<Real(const Real& time)>& a_voltageCurve) noexcept;

This is relevant only when running a transient simulation. 

Algorithms
----------

The streamer inception model uses a combination of electrostatic field solves, Particle-In-Cell, and fluid advection for resolving the necessary dynamics.
The various algorithms involved are discussed below.

Field solve
___________

Since the background field scales linearly with applied voltage, we require only a single field solve at the beginning of the simulation.
This field solve is done with an applied voltage of :math:`U = 1\,\text{V}` and the electric field is then simply later scaled by the actual voltage.

Inception integral
__________________

We use a Particle-In-Cell method for computing the inception integral :math:`K\left(\mathbf{x}\right)` for an arbitrary electron starting position.
All grid cells where :math:`\alpha_{\textrm{eff}} > 0` are seeded with one particle on the cell centroid and the particles are then tracked through the grid.
The particles move a user-specified distance along field lines :math:`\mathbf{E}` and the particle weights are updated using first or second order integration.
If a particle leaves through a boundary (EB or domain boundary), or enters a region :math:`\alpha_{\text{eff}} \leq 0`, the integration is stopped.
Once the particle integration halts, we rewind the particles back to their starting position and deposit their weight on the mesh, which provides us with :math:`K = K\left(\mathbf{x}\right)`.

Euler
^^^^^

For the Euler rule the particle weight for a particle :math:`p` the update rule is

.. math::

   \mathbf{x}_p^{k+1} = \mathbf{x}_p^k - \mathbf{\hat{E}}\left(\mathbf{x}_p^k\right)\Delta x
   
   w_p^{k+1} = w_p^k + \alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{x}_p^k\right)\right|,\mathbf{x}_p^k\right)\Delta x,

where :math:`\Delta x` is a user-specified integration length.

Trapezoidal
^^^^^^^^^^^

With the trapezoidal rule the update is first

.. math::

   \mathbf{x}_p^\prime = \mathbf{x}_p^k - \mathbf{\hat{E}}\left(\mathbf{x}_p^k\right)\Delta x

followed by


.. math::

      \mathbf{x}_p^{k+1} = \mathbf{x}_p^k + \frac{\Delta x}{2}\left[\mathbf{\hat{E}}\left(\mathbf{x}_p^k\right) + \mathbf{\hat{E}}\left(\mathbf{x}_p^\prime\right)\right].

      w_p^{k+1} = w_p^k + \frac{\Delta x}{2}\left[\alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{x}_p^k\right)\right|,\mathbf{x}_p^k\right) + \alpha_{\text{eff}}\left(\left|\mathbf{E}\left(\mathbf{x}_p^\prime\right)\right|,\mathbf{x}_p^\prime\right)\right]

.. note::

   When tracking positive ions for evaluation of the Townsend criterion, the same algorithms are used.
   The exception is that the positive ions are simply tracked along field lines until they strike a cathode, so that there is no integration with respect to :math:`\alpha_{\text{eff}}`.


Critical volume
_______________

The critical volume is computed as

.. math::

   V_c = \int_{K\left(\mathbf{x}\right) > K \cup \gamma\exp\left[K\left(\mathbf{x}\right)\right] \ge 1} \text{d}V.

Note that the critical volume is both voltage and polarity dependent.

Critical surface
________________

The critical surface is computed as

.. math::

   A_c = \int_{K\left(\mathbf{x}\right) > K \cup \gamma\exp\left[K\left(\mathbf{x}\right)\right] \ge 1} \text{d}A.

Note that the critical surface is both voltage and polarity dependent, and is non-zero only on cathode surfaces.

Inception voltage
_________________

Arbitrary starting electron
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The inception voltage for starting a critical avalanche can be computed in the stationary solver mode.
In this case we compute :math:`K\left(\mathbf{x}; U\right)` for a range of voltages :math:`U \in U_1, U_2, \ldots`.

If two values of the :math:`K` integral bracket :math:`K_c`, i.e.

.. math::

   K_a = K\left(\mathbf{x}; U_a\right) \leq K_c

   K_b = K\left(\mathbf{x}; U_b\right) \geq K_c

then we can estimate the inception voltage for a starting electron at position :math:`\mathbf{x}` through linear interpolation as

.. math::

   U_{\text{inc, streamer}}\left(\mathbf{x}\right) = U_a + \frac{K_c - K_a}{K_b - K_a}\left(U_b - U_a\right)

A similar method is used for the Townsend criterion, using e.g. :math:`T\left(\mathbf{x}; U\right) = \gamma\exp\left[K\left(\mathbf{x}; U\right)\right]`, then if


.. math::
   
   T_a = T\left(\mathbf{x}; U_a\right) \leq 1

   T_b = T\left(\mathbf{x}; U_b\right) \ge 1

then we can estimate the inception voltage for a starting electron at position :math:`\mathbf{x}` through linear interpolation as

.. math::

   U_{\text{inc, Townsend}}\left(\mathbf{x}\right) = U_a + \frac{1 - T_a}{T_b - T_a}\left(U_b - U_a\right)

The inception voltage for position :math:`\mathbf{x}` is then

.. math::

   U_{\text{inc}} = \min\left[U_{\text{inc, streamer}}\left(\mathbf{x}\right), U_{\text{inc, Townsend}}\left(\mathbf{x}\right)\right]
   

Minimum inception voltage
^^^^^^^^^^^^^^^^^^^^^^^^^

The minium inception voltage is the minimum voltage required for starting a critical avalanche (or Townsend process) for an arbitrary starting electron.
From the above, this is simply

.. math::

   U_{\text{inc}}^{\text{min}} = \min_{\forall \mathbf{x}} \left[U_{\text{inc}}\left(\mathbf{x}\right)\right].

From the above we also determine

.. math::

   \mathbf{x}_{\text{inc}}^{\text{min}} \leftarrow \mathbf{x}\text{ that minimizes } U_{\text{inc}}\left(\mathbf{x}\right) \forall \mathbf{x},

which is the position of the first electron that enables a critical avalanche at the minimum inception voltage.

.. note::

   The minimum inception voltage is the minimum voltage required for starting a critical avalanche.
   However, as :math:`U \rightarrow U_{\text{inc}}^{\text{min}}` we also have :math:`V_c \rightarrow 0`, requires the a starting electron *precisely* in :math:`\mathbf{x}_{\text{inc}}^{\text{min}}`.

Inception probability
_____________________

The inception probability is given by :eq:`StreamerInceptionProbability` and is computed using straightforward numerical quadrature:

.. math::

   \int_{V_c}\left\langle\frac{dn_{\text{e}}}{dt}\right\rangle\left(1-\frac{\eta}{\alpha}\right) \text{d}V \approx \sum_{\mathbf{i}\in K_\mathbf{i} > K_c} \left(\left\langle\frac{dn_{\text{e}}}{dt}\right\rangle\right)_{\mathbf{i}}\left(1 - \frac{\eta_{\mathbf{i}}}{\alpha_{\mathbf{i}}}\right)\kappa_{\mathbf{i}}\Delta V_{\mathbf{i}},

and similarly for the surface integral.

.. important::

   The integration runs over *valid cells*, i.e. grid cells that are not covered by a finer grid.

Advection algorithm
___________________

The advection algorithm for the negative ion distribution follows the time stepping algorithms described in the advection-diffusion model, see :ref:`Chap:AdvectionDiffusionModel`.


Simulation control
------------------

Here, we discuss simulation controls that are available for the streamer inception model.
These all appear in the form ``StreamerInceptionStepper.<option>``.

verbosity
_________

The ``verbosity`` input option controls the model chattiness (to the ``pout.*`` files).
Usually we have

.. code-block:: text

   StreamerInceptionStepper.verbosity = -1

mode
____

The mode flag switches between stationary and transient solves.
Accepted values are ``stationary`` and ``transient``, e.g.,

.. code-block:: text

   StreamerInceptionStepper.mode = stationary

.. important::

   When running in stationary mode, set ``Driver.max_steps=0``. 


inception_alg
_____________

Controls the streamer inception algorithm (for computing the :math:`K` integral).
This should be specified in the form

.. code-block:: text

   StreamerInceptionStepper.inception_alg = <algorithm> <mode> <value>

These indicate the following:

* ``<algorithm>`` indicates the integration algorithm.
  Currently supported is ``trapz`` (trapezoidal rule) and ``euler``.

* ``mode`` indicates the integration step size selection.
  This can be the following:
  
  * ``fixed`` for a spatially fixed step size.
  * ``dx`` for step sizes relative to the grid resolution :math:`\Delta x`.
  * ``alpha`` For setting the step size relative to the avalanche length :math:`1/\alpha`.and mode is either ``fixed`` or ``dx``.

  Normally, ``alpha`` will yield the best integration results since the step size is adaptively selected, taking large steps where :math:`\alpha` is small and smaller steps where :math:`\alpha` is large.

* ``value`` indicates a step size, and has a different interpretation for the various modes.
  * If using ``fixed`` integration, ``value`` indicates the physical length of the step size.
  * If using ``dx`` integration,  ``value`` indicates the step size relative to the grid cell resolution.
  * If using ``alpha`` integration, ``value`` indicates the step size relative to the avalanche length :math:`1/\alpha`.

For example, the following will set an Euler integration with an adaptive step size:

.. code-block:: text

   StreamerInceptionStepper.inception_alg = euler alpha 0.5

   
full_integration
________________

Normally, it will not necessary to integrate the particles beyond :math:`w > K_c` since this already implies inception.
The flag ``full_integration`` can be used to turn on/off integration beyond :math:`K_c`.
If the flag is set to false, the particle integration routine will terminate once a particle weight reaches :math:`K_c`.
If the flag is set to true, the particle integration routine will proceed until the particles leave the domain or ionization volumes. 

.. tip::

   Setting ``full_integration`` to false can lead to large computational savings when the ionization volumes are large.


output_file
___________

Controls the overall report file for stationary and transient solves.
The user specifies a filename for a file which will be created (in the same directory as the application is running), containing a summary of the most important simulation output variables.

.. warning::

   Running a new simulation will overwrite the specified ``output_file``. 

For example:

.. code-block:: text

   StreamerInceptionStepper.output_file = report.txt

K_inception
___________

Controls the critical value of the :math:`K` integral.
E.g.,

.. code-block:: text

   StreamerInceptionStepper.K_inception = 12

eval_townsend
_____________

Controls whether or not the Townsend criterion is also evaluated.

E.g.,

.. code-block:: text

   StreamerInceptionStepper.eval_townsend = true

Will turn on the Townsend criterion. 

plt_vars
________

Controls plot variables that will be written to HDF5 outputs in the :file:`plt` folder. 
Valid options are

* ``K``        - Inception integral
* ``Uinc``     - Inception voltage
* ``bg_rate``  - Background ionization rate
* ``emission`` - Field emission
* ``alpha``    - Effective ionization coefficient
* ``eta``      - Eta coefficient  
* ``poisson``  - Poisson solver
* ``tracer``   - Tracer particle solver
* ``ions``     - Ion solver
* ``T``        - Townsend criterion

For example:

.. code-block:: text

   StreamerInceptionStepper.plt_vars = K Uinc bg_rate emission ions

For stationary mode
____________________

For the stationary mode the following input flags are required:

* ``voltage_lo`` Lowest simulated voltage. 
* ``voltage_hi`` High simulated voltage. 
* ``voltage_steps`` Extra voltage steps between ``voltage_lo`` and ``voltage_hi``.

These voltages levels are used when running a stationary solve.   
For example:

.. code-block:: text

   StreamerInceptionStepper.voltage_lo    = 10E3
   StreamerInceptionStepper.voltage_hi    = 30E3
   StreamerInceptionStepper.voltage_steps = 5

For transient mode
__________________

For the transient mode the following input options must be set:

* ``ion_transport`` For turning on/off ion transport.
* ``transport_alg`` For controlling the transport algorithm.
  Valid options are ``euler``, ``heun``, or ``imex`` (for semi-implicit with corner transport upwind).
* ``cfl`` Which controls the ion advection time step.
* ``min_dt`` For setting the minimum time step used.
* ``max_dt`` For setting the maximum time step used.

For example,

.. code-block:: text
		
   StreamerInceptionStepper.ion_transport = true 
   StreamerInceptionStepper.transport_alg = imex  
   StreamerInceptionStepper.cfl           = 0.8  
   StreamerInceptionStepper.min_dt        = 0.0  
   StreamerInceptionStepper.max_dt        = 1E99 

.. warning::

   The ``ctu`` option exists because the default advection solver for the streamer inception model is the corner transport upwind solver (see :ref:`Chap:CdrCTU`).
   Ensure that ``CdrCTU.use_ctu = true`` if using ``StreamerInceptionStepper.transport_alg = ctu`` algorithm and set ``CdrCTU.use_ctu = false`` otherwise.

  
Caveats
_______

The model is intended to be used with a nearest-grid-point deposition scheme (which is also volume-weighted).
When running the model, ensure that the :ref:`Chap:TracerParticleSolver` flag is set as follows:

.. code-block:: text

   TracerParticleSolver.deposition   = ngp    

Adaptive mesh refinement
------------------------

The streamer inception model runs its own mesh refinement routine, which refines the mesh if

.. math::

   \alpha_{\text{eff}}\left(\left|\mathbf{E}\right|, \mathbf{x}\right)\Delta x > \lambda,

where :math:`\lambda` is a user-specified refinement criterion.

This is implemented in a class

.. code-block:: c++

   class StreamerInceptionTagger : public CellTagger

and is automatically included in simulations when setting up the application through the Python setup tools (see :ref:`Chap:StreamerInceptionSetup`).
The user can control refinement buffers and criterion through the following input options:

* ``StreamerInceptionTagger.buffer`` Adds a buffer region around tagged cells.
* ``StreamerInceptionTagger.max_voltage`` Maximum voltage that will be simulated.
* ``StreamerInceptionTagger.ref_alpha`` Sets the refinement criterion :math:`\lambda` as above.

For example:

.. code-block:: text
		
   StreamerInceptionTagger.buffer      = 4  
   StreamerInceptionTagger.max_voltage = 30E3
   StreamerInceptionTagger.ref_alpha   = 2.0

.. _Chap:StreamerInceptionSetup:

Setting up a new problem
------------------------

To set up a new problem, using the Python setup tools in :file:`$DISCHARGE_HOME/Physics/StreamerInception` is the simplest way.
To see available setup options, run

.. code-block:: text

   python3 setup.py --help

For example, to set up a new problem in :file:`$DISCHARGE_HOME/MyApplications/MyStreamerInception` for a cylinder geometry, run

.. code-block:: text

   python3 setup.py -base_dir=MyApplications -app_name=MyStreamerInception -geometry=Cylinder

This will set up a new problem in a cylinder geometry (defined in :file:`Geometries/Cylinder`).
The main file is named :file:`program.cpp`` and contains default implementations for the required input data (see :ref:`Chap:StreamerInceptionInputData`).


Example programs
----------------

Example programs that use the streamer inception model are given in

High-voltage vessel
___________________

* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/Vessel`.
  This program is set up in 2D (stationary) and 3D (transient) for streamer inception in atmospheric air.
  The input data is computed using BOLSIG+.

Electrode with surface roughness
________________________________

* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/ElectrodeRoughness`.
  This program is set up in 2D (stationary) and 3D (transient) for streamer inception on an irregular electrode surface. 
  We use SF6 transport data as input data, computed using BOLSIG+.

Electrode with surface roughness
________________________________

* :file:`$DISCHARGE_HOME/Exec/Examples/StreamerInception/ElectrodeRoughness`.
  This program is set up in 2D and 3D (stationary) mode, and includes the influence of the Townsend criterion.
