.. _Chap:CdrPlasmaModel:

CDR plasma model
****************

In the CDR plasma model we are solving


.. math::
   :label: CdrPlasmaEquations
      
   &\nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0},\\
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\
   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n\right) = S,
   
The above equations must be supported by additional boundary conditions on electrodes and insulating surfaces. 

Radiative transport can be done either in the diffusive approximation or by means of Monte Carlo methods (see :ref:`Chap:RtSolver`).
Diffusive RTE methods (see :ref:`Chap:EddingtonSP1`) involve solving

.. math::
   
   \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},
   
where :math:`\Psi` is the isotropic photon density, :math:`\kappa` is an absorption length and :math:`\eta` is an isotropic source term.
I.e., :math:`\eta` is the number of photons produced per unit time and volume.
The time dependent term can be turned off and the equation can be solved stationary.

The module also supports discrete photons where photon transport and absorption is done by sampling discrete photons (see :ref:`Chap:MonteCarloRTE`).
In general, discrete photon methods incorporate better physics (like shadows).
They can also be more easily adapted to scattering media.
They are, on the other hand, inherently stochastic which implies that some extra caution must be exercised when integrating the equations of motion.

The coupling that is (currently) available in ``chombo-discharge`` is

.. math::
   :label: CdrPlasmaCoupling

   \epsilon_r =& \epsilon_r(\mathbf{x}),\\
   \mathbf{v} =& \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right),\\
   D =& \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right),\\
   S =& S\left(t, \mathbf{x}, \mathbf{E}, \nabla\mathbf{E}, n, \nabla n, \Psi\right),\\
   \eta =& \eta\left(t, \mathbf{x}, \mathbf{E}, n\right),\\
   F =& F(t, \mathbf{x}, \mathbf{E}, n),


where :math:`F` is the boundary flux on insulators or electrodes (which must be separately implemented).

``chombo-discharge`` works by embedding the equations above into an abstract C++ framework (see :ref:`Chap:CdrPlasmaPhysics`) that the user must implement or reuse existing pieces of, and then compile into an executable.

.. tip::
   
   The CDR plasma model resides in :file:`/Physics/CdrPlasma` and describes plasmas in the drift-diffusion approximation.
   This physics model also includes the following subfolders:

   * :file:`/Physics/CdrPlasma/PlasmaModel` which contains various implementation of some plasma models that we have used.
   * :file:`/Physics/CdrPlasma/TimeSteppers` contains various algorithms for advancing the equations of motion. 
   * :file:`/Physics/CdrPlasma/CellTaggers` contains various algorithms for flagging cells for refinement and coarsening.

   See `CdrPlasmaStepper <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhysics_1_1CdrPlasma_1_1CdrPlasmaStepper.html>`_ for the overall C++ API.

This module uses the following solvers:

#. Advection-diffusion-reaction solver, :ref:`Chap:CdrSolver`.
#. Electrostatics solvers, :ref:`Chap:FieldSolver`.
#. Radiative transfer solver (either Monte-Carlo or continuum approximation), :ref:`Chap:RtSolver`.
#. Surface charge solver, see :ref:`Chap:SurfaceODESolver`.

Time discretizations
====================

Here, we discuss two discretizations of :eq:`CdrPlasmaEquations`.
Firstly, note that there are two layers to the time integrators:

#. A pure class ``CdrPlasmaStepper`` which inherits from ``TimeSteppers`` but does not implement an ``advance`` method.
   This class simply provides the base functionality for more easily developing time integrators.
   ``CdrPlasmaStepper`` contains methods that are necessary for coupling the solvers, e.g. calling the :ref:`Chap:CdrPlasmaPhysics` methods at the correct time.
#. Implementations of ``CdrPlasmaPhysics``, which implement the ``advance`` method and can thus be used for advancing models.
   
The supported time integrators are located in :file:`$DISCHARGE_HOME/CdrPlasma/TimeSteppers`.
There are two integrators that are commonly used.

* A Godunov operator splitting with either explicit or implicit diffusion.
  This integrator also supports semi-implicit formulations. 
* A spectral deferred correction (SDC) integrator with implicit diffusion.
  This integrator is an implicit-explicit.

Briefly put, the Godunov operator is our most stable integrator, while the SDC integrator is our most accurate integrator. 

..
   Time step limitations
   ---------------------

   For explicit advection the time step limit is

   .. math::
      :label: dtA

      \Delta t = \frac{\Delta x}{\sum_{i=1}^{\textrm{d}} |v_i|},

   where :math:`d = 2,3` is the spatial dimension. 
   For explicit diffusion with constant diffusion coefficient :math:`D` the time step limit is

   .. math::
      :label: dtD

      \Delta t = \frac{\Delta x^2}{2D \textrm{d}}.

   For explicit advection-diffusion the time step limit is

   .. math::
      :label: dtAD

      \Delta t \leq \frac{1}{\sum_{i=1}^{\textrm{d}} \frac{|v_i|}{\Delta x} + \frac{2D\textrm{d}}{\Delta x^2}}.

   Note that the resolution :math:`\Delta x` dictates whether or not one should run with implicit diffusion or not.
   Implicit diffusion requires solving at least one extra Helmholtz equation for the diffusive species, but the time step can also be larger. 

.. _Chap:godunov:

Godunov operator splitting
--------------------------

The ``CdrPlasmaGodunovStepper`` implements ``CdrPlasmaStepper`` and defines an operator splitting method between charge transport and plasma chemistry.
It has a formal order of convergence of one.
The source code is located in :file:`$DISCHARGE_HOME/Physics/CdrPlasma/TimeSteppers/CdrPlasmaGodunovStepper`.

.. warning::

   Splitting the terms yields *splitting errors* which can dominate for large time steps.
   Typically, the operator splitting discretization is not suitable for large time steps. 

The basic advancement routine for ``CdrPlasmaGodunovStepper`` is as follows:

#. Advance the charge transport :math:`\phi^k \rightarrow \phi^{k+1}` with the source terms set to zero.
#. Compute the electric field.
#. Advance the plasma chemistry over the same time step using the field computed above
   I.e., advance :math:`\partial_t\phi = S` over a time step :math:`\Delta t`.  
#. Advance the radiative transport part.
   This can also involve discrete photons.

The transport/field steps can be done in various ways:
The following transport algorithms are available:

* **Euler**, where everything is advanced with the Euler rule.
* **Semi-implicit**, where the Euler field/transport step is performed with an implicit coupling to the electric field.

In addition, diffusion can be treated

* **Explicitly**, where all diffusion advances are performed with an *explicit* rule.
* **Implicitly**, where all diffusion advances are performed with an *implicit* rule.
* **Automatically**, where diffusion advances are performed with an implicit rule only if time steps dictate it, and explicitly otherwise.

.. note::

   When setting up a new problem with the Godunov time integrator, the default setting is to use automatic diffusion and a semi-implicit coupling.
   These settings tend to work for most problems.

Specifying transport algorithm
______________________________

To specify the transport algorithm, modify the flag ``CdrPlasmaGodunovStepper.transport``, and set it to ``semi_implicit`` or ``euler``.
Everything else is an error.

Note that for the Godunov integrator, it is possible to center the advective discretization at the half time step.
That is, the advancement algorithm is

.. math::

   n^{k+1} = n^{k} - \nabla\cdot\left(n^{k+1/2}\mathbf{v}\right) + \nabla\cdot\left(D\nabla\phi^k\right),

where :math:`n^{k+1/2}` is obtained by also including transverse slopes (i.e., extrapolation in time).
See :cite:t:`trebotich2015` for details.
Note that the formal order of accuracy is still one, but the accuracy of the advective discretization is increased substantially.

Specifying diffusion
____________________

To specify how diffusion is treated, modify the flag ``CdrPlasmaGodunovStepper.diffusion``, and set it to ``auto``, ``explicit``, or ``implicit``.
In addition, the flag ``CdrPlasmaGodunovStepper.diffusion_thresh`` must be set to a number.

When diffusion is set to ``auto``, the integrator switches to implicit diffusion when

.. math::

   \frac{\Delta t_{\textrm{A}}}{\Delta t_{\textrm{AD}}} > \epsilon,

where :math:`\Delta t_{\textrm{A}}` is the advection-only limited time step and :math:`\Delta t_{\textrm{AD}}` is the advection-diffusion limited time step.

.. note::

   When there are multiple species being advected and diffused, the integrator will perform extra checks in order to maximize the time steps for the other species.

Time step limitations
_____________________

The basic time step limitations for the Godunov integrator are:

* Manually set maximum and minimum time steps
* Courant-Friedrichs-Lewy conditions, either on advection, diffusion, or both.
* The dielectric relaxation time.

The user is responsible for setting these when running the simulation.
Note when the semi-implicit scheme is used, it is not necessary to restrict the time step by the dielectric relaxation time.

Solver configuration
____________________

``CdrPlasmaGodunovStepper`` contains several options that can be configured when running the solver, which are listed below:

.. literalinclude:: ../../../../Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.options

.. _Chap:SISDC:

Spectral deferred corrections
-----------------------------

The ``CdrPlasmaImExSdcStepper`` uses implicit-explicit (ImEx) spectral deferred corrections (SDCs) to advance the equations.
This integrator implements the ``advance`` method for ``CdrPlasmStepper``, and is a high-order method with implicit diffusion.

SDC basics
__________

First, we provide a quick introduction to the SDC procedure. 
Given an ordinary differential equation (ODE) as

.. math::
   \frac{\partial u}{\partial t} = F(u,t), \quad u(t_0) = u_0,

the exact solution is

.. math::
   u(t) = u_0 + \int_{t_0}^tF\left(u,\tau\right)d\tau.

Denote an approximation to this solution by :math:`\widetilde{u}(t)` and the correction by :math:`\delta(t) = u(t) - \widetilde{u}(t)`. The measure of error in :math:`\widetilde{u}(t)` is then

.. math::
   R(\widetilde{u}, t) = u_0 + \int_{t_0}^tF(\widetilde{u}, \tau)d\tau - \widetilde{u}(t).

Equivalently, since :math:`u = \widetilde{u} + \delta`, we can write

.. math::
   \widetilde{u} + \delta = u_0 + \int_{t_0}^t F\left(\widetilde{u}+\delta, \tau\right)d\tau. 

This yields

.. math::
   \delta = \int_{t_0}^t\left[F\left(\widetilde{u}+\delta, \tau\right) - F\left(\widetilde{u}, \tau\right)\right]d\tau + R\left(\widetilde{u},t\right). 

This is called the correction equation. The goal of SDC is to iteratively solve this equation in order to provide a high-order discretization.

The ImEx SDC method in ``chombo-discharge`` uses implicit diffusion in the SDC scheme.
Coupling to the electric field is always explicit.
The user is responsible for specifying the quadrature nodes, as well as setting the number of sub-intervals in the SDC integration and the number of corrections.
In general, each correction raises the discretization order by one.

Time step limitations
_____________________

The ImEx SDC integrator is limited by

* The dielectric relaxation time.
* An advective CFL conditions.

In addition to this, the user can specify maximum/minimum allowed time steps.

..
   We now discuss the explicit-implicit SDC method.
   First, we apply the method of lines (MOL) such that

   .. math::
      :nowrap:

      \begin{eqnarray}
      \frac{d\phi_{\mathbf{i}}}{dt} &= \mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}\right) + \mathcal{F}_{\textrm{D}}\left(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}\right), \\
      \frac{d\sigma_{\mathbf{i}}}{dt} &= \mathcal{F}_{\sigma}\left(t, \phi_{\mathbf{i}}\right),
      \end{eqnarray}

   where :math:`\phi_{\mathbf{i}}` denotes a cell-averaged variable, :math:`\mathcal{F}_{\sigma}` is the surface charge flux, :math:`\mathcal{F}_{\textrm{AR}}` is the advection-reaction operator, and :math:`\mathcal{F}_{\textrm{D}}` is the diffusion operator.

   **SISDC predictor**

   In what follows, we suppress the index :math:`{\mathbf{i}}` as it is not explicitly needed.
   Given an interval :math:`[t_n, t_{n+1}]` on which a solution is sought, SDC methods divide this interval into :math:`p` subintervals :math:`t_n = t_{n,0} < t_{n,1} < \ldots < t_{n,p} = t_{n+1}`.
   Our discussion pertains only to the interval :math:`[t_n, t_{n+1}]` so we compress the notation to :math:`t_m\equiv t_{n,m}`.
   First, we obtain predictor solution :math:`\phi_{m}^0, m=0,1,\ldots,p` as the semi-implicit advance

   .. math::
      :nowrap:

      \begin{eqnarray}
      \phi_{m+1}^0 &= \phi_m^0 + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m,\phi_m^0\right) + \mathcal{F}_{\textrm{D}}\left(t_{m+1},\phi_{m+1}^0; \mathbf{E}_{m+1}^0\right)\right],\\
      \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m,\phi_m^0\right).
      \end{eqnarray}

   This defines a Helmholtz problem for :math:`\phi_{m+1}^0` through :math:`\mathcal{F}_{\textrm{D}}`. Generally, the upper subscript denotes an SDC iteration where subscript 0 is the SISDC predictor, and we also have :math:`\phi_0^0 = \phi(t_n)` and :math:`\sigma_0^0 = \sigma(t_n)`. This predictor treats advection and chemistry terms explicitly, and diffusion implicitly. Other types of semi-implicit or multi-implicit couplings are possible :cite:`Bourlioux2003,Layton2004,Nonaka2012`. SDC improves this solution by using deferred corrections: Given a numerical solution :math:`\phi_{m+1}^k`, we compute an error :math:`\delta_{m+1}^k` and obtain the next iterate :math:`\phi_{m+1}^{k+1} = \phi_{m+1}^k + \delta_{m+1}^k`. Each iteration raises the discretization order by one :cite:`Dutt2000,Minion2003`, to maximum order :math:`p+1`. Critical to the success of this approach is the precise evaluation of the numerical quadrature. 

   The parametric coupling of the electric field complicates things since the predictor contains :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0\right)`, implying that the Poisson equation and the diffusion advance require concurrent solves for the diffusion update. We simplify this system by using a weak coupling by first computing

   .. math::
      :nowrap:

      \begin{eqnarray}
      \phi_{m+1}^{0,\ast} &= \phi_m^0 + \Delta t_m\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^0\right), \\
      \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m, \phi_m^0\right),
      \end{eqnarray}

   Next, we will approximate :math:`\mathbf{E}_{m+1}^{0}` for use in the predictor. There are two choices for this coupling; one may either use :math:`\mathbf{E}_m^0` for computation of the diffusion coefficients, which we will refer to as the semi-implicit coupling, or one may use fixed-point iteration and compute :math:`\mathbf{E}_{m+1}^{0,\ast} = \mathbf{E}\left(\phi_{m+1}^{0, \ast}, \sigma_{m+1}^0\right)`, followed by the diffusion advance

   .. math::
      \phi_{m+1}^{0,\dagger} = \phi_{m+1}^{0,\ast} + \Delta t_m\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^0; \mathbf{E}_{m+1}^\ast\right),

   which we will refer to as the implicit coupling. This is e.g. the electric field coupling used in :cite:`Marskar2019`. This approximation can be improved by using more fixed-point iterations that computes :math:`\mathbf{E}_{m+1}^{0,\dagger} = \mathbf{E}\left(\phi_{m+1}^{0,\dagger}, \sigma_{m+1}^0\right)` and then re-solves the predictor equation with :math:`\mathbf{E}_{m+1}^{0,\dagger}` in place of :math:`\mathbf{E}_{m+1}^{0,\ast}`. The process can then be repeated for increased accuracy. Regardless of which coupling is used, we have now calculated :math:`\phi_{m+1}^0`, :math:`\sigma_{m+1}^0`, through which we obtain :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0, \sigma_{m+1}^0\right)`, and :math:`\Psi_{m+1}^0 = \Psi\left(\mathbf{E}_{m+1}^0, \phi_{m+1}^0\right)`. Finally, we remark that the SISDC predictor is a sequentially advanced semi-implicit Euler method, which is locally second order accurate and globally first order accurate. Each step of the predictor can be thought of as a Godunov splitting between the advective-reactive and diffusive terms. 

   SISDC corrector
   ^^^^^^^^^^^^^^^
   Next, the semi-implicit discretization of the correction equation is

   .. math::
      \begin{split}
      \delta_{m+1}^k &= \delta_m^k  + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k + \delta_m^k\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)\right.\\
      &+ \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k + \delta_{m+1}^k; \mathbf{E}_{m+1}^k\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k; \mathbf{E}_{m+1}^k\right)\right] - \left(R_{m+1}^k - R_{m}^k\right).
      \end{split}

   We furthermore define

   .. math::
      \begin{split}
      R_{m+1}^k - R_m^k &= \int_{t_m}^{t_{m+1}}\left[\mathcal{F}_{\textrm{AR}}\left(\phi^k\right) + \mathcal{F}_{\textrm{D}}\left(\phi^k; \mathbf{E}^k\right)\right]d\tau - \phi_{m+1}^k + \phi_m^k \\
      &\equiv I_m^{m+1}\left(\phi^k\right) - \phi_{m+1}^k + \phi_m^k. 
      \end{split} 

   Evaluation of :math:`I_m^{m+1}` yields :math:`p` quadrature rules and we may write

   .. math::
      I_m^{m+1}\left(\phi^k\right) = \sum_{l=0}^p q_m^l\left[\mathcal{F}_{\textrm{AR}}\left(t_l, \phi^k_l\right) + \mathcal{F}_{\textrm{D}}\left(t_l, \phi^k_l; \mathbf{E}_l^k\right)\right],

   where the weights :math:`q_m^l` are quadrature weights. The final update for :math:`\phi^{k+1}_{m+1}` is then

   .. math::
      \begin{split}
      \phi_{m+1}^{k+1} &= \phi_{m}^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) -\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right.\\
      & + \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \phi_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right] + I_{m}^{m+1}\left(\phi^k\right).
      \end{split}

   With the exception of :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right)`, all quantities on the right-hand are known and the correction equation is reduced to a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` with error :math:`\delta_{m+1}^k = \phi_{m+1}^{k+1} - \phi_{m+1}^k`. An analogous equation is found for :math:`\sigma_{m+1}^{k+1}`.

   The correction step has the same coupling to the electric field as the prediction step in that :math:`\mathbf{E}_{m+1}^{k+1}` appears in the update equation for :math:`\phi_{m+1}^{k+1}`. As for the prediction, we use a weak coupling through which we first compute

   .. math::
      :nowrap:

      \begin{eqnarray}
      \phi_{m+1}^{k+1,\ast} &= \phi_m^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right] + I_m^{m+1}\left(\phi^k\right),\\
      \sigma_{m+1}^{k+1} &= \sigma_m^{k+1} + \Delta t_m\left[F_\sigma\left(t_m, \phi_m^{k+1}\right) - F_\sigma\left(t_m, \phi_m^{k}\right)\right] + \Sigma_m^{m+1}\left(\phi^k\right). 
      \end{eqnarray}

   The solution for :math:`\sigma_{m+1}^{k+1}` is final since all charge is injected through the advection operator for :math:`\phi`. The term :math:`\Sigma_m^{m+1}` contains the injected charge through :math:`I_m^{m+1}\left(\phi^k\right)`, as was discussed in :ref:`Chap:SpatialDiscretization`. We then solve

   .. math::
      \phi_{m+1}^{k+1} = \phi_{m+1}^{k+1, \ast} + \Delta t_m\left[\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right],

   with some approximation for :math:`\mathbf{E}_{m+1}^{k+1}`. As before, this coupling can be made either semi-implicitly or implicitly. The corrector equation defines a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` using :math:`\phi_{m+1}^{k+1,\ast}` as the previous solution and :math:`-\mathcal{F}_{\textrm{D}}\left(\phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)` as a source term.

   Order, stability, and computational cost
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   For consistency with the literature, denote the SISDC method which uses :math:`P` nodes (i.e. :math:`P-1` subintervals) and :math:`K` total iterations (i.e. :math:`K-1` iterations of the correction equation) by :math:`\verb|SISDC|_P^K`. This method will have a global order of accuracy :math:`\min\left(K,P\right)` if the quadrature can be evaluated with appropriate accuracy. Order reductions may occur if the interpolating polynomial in the quadrature suffers from Runge's phenomenon. As we discuss below, uniformly spaced nodes have some computational advantage but is therefore also associated with some risk. Safer choices include Lobatto nodes or Chebyshev nodes (with inclusion of endpoints) to minimize the risk of order reductions. Implications on the choice of quadrature nodes can be found in :cite:`Layton2005`. 

   For explicit advection, the deferred correction procedure integrates the correction equation sequentially and therefore does not allow each substep :math:`\Delta t_m` to exceed the CFL-limited time step :math:`\Delta t_{\textrm{cfl}}`, i.e. :math:`\Delta t_m < \Delta t_{\textrm{cfl}} \forall m`. Since we have :math:`\Delta t = \sum_m\Delta t_m`, uniform nodes maximize :math:`\Delta t` subject to the CFL constraint. For example, an :math:`\verb|SISDC|_P^K` method with uniformly spaced nodes has a maximum possible time step :math:`\Delta t < (P-1)\Delta t_{\textrm{cfl}}`. For the same number of function evaluations, the allowed time step with for Lobatto or Chebyshev nodes is smaller. For :math:`P\leq 3`, the uniform nodes, Lobatto nodes, and Chebyshev nodes coincide. Larger time steps are possible with uniform nodes for :math:`P>3`, which has some computational consequence. The table below summarizes the largest possible time steps for the :math:`\verb|SISDC|_P^K` method with the various quadratures. Finally, note that :math:`\Delta t_m < \Delta t_{\textrm{cfl}}` does not guarantee stability since further restrictions may required for stability of the reaction terms.

   ==========  =================================== ====================================   ================================
    :math:`P`   Lobatto                             Chebyshev                             Uniform
   ==========  =================================== ====================================   ================================
   2           :math:`\Delta t_{\textrm{cfl}}`      :math:`\Delta t_{\textrm{cfl}}`       :math:`\Delta t_{\textrm{cfl}}`
   3           :math:`2\Delta t_{\textrm{cfl}}`     :math:`2\Delta t_{\textrm{cfl}}`      :math:`2\Delta t_{\textrm{cfl}}`
   4           :math:`2.26\Delta t_{\textrm{cfl}}`  :math:`1.73\Delta t_{\textrm{cfl}}`   :math:`3\Delta t_{\textrm{cfl}}`
   5           :math:`3.05\Delta t_{\textrm{cfl}}`  :math:`2.82\Delta t_{\textrm{cfl}}`   :math:`4\Delta t_{\textrm{cfl}}`
   6           :math:`3.50\Delta t_{\textrm{cfl}}`  :math:`3.29\Delta t_{\textrm{cfl}}`   :math:`5\Delta t_{\textrm{cfl}}`
   7           :math:`4.26\Delta t_{\textrm{cfl}}`  :math:`4.36\Delta t_{\textrm{cfl}}`   :math:`6\Delta t_{\textrm{cfl}}`
   ==========  =================================== ====================================   ================================

   For the predictor step, it is necessary to evaluate :math:`\mathcal{F}_{\textrm{AR}}\left(\phi_m^{k+1}\right)` and thus update the Poisson and radiative transfer equations at each node. In addition, it is necessary to solve the diffusion equation at every node except :math:`m=0` for every diffusive species, which may also require auxiliary updates of the electric field. The corrector step contains extra floating point operator due to the extra terms :math:`\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)` and :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k\right)` and the quadrature :math:`I_m^{m+1}`. The computational cost of adding in these terms is small compared to the cost of an Euler update of the advection-reaction equation since one must also computate source terms, drift velocities, and boundary conditions in addition to construction of the hybrid divergence. In short, the computational cost of the predictor and corrector steps are about the same.

   Next, we provide some remarks on the extra computational work involved for higher order methods. Broadly speaking, the total amount of floating point operations increases quadratically with the order. Each node requires evaluation of one advection-reaction operator, at least one electric field update, and one radiative transfer update. Likewise, each substep requires one diffusion solve. Thus, :math:`\verb|SISDC|_K^K` requires :math:`K^2` advection-reaction evaluations, :math:`(K-1)^2` diffusion solves, :math:`(K-1)^2` radiative transfer updates, and at least :math:`K^2` electric field updates. In these estimates we have assumed that the diffusion solve couples semi-implicitly to the electric field, thus each corrector iteration requires one electric field update per node, giving a total cost :math:`K^2`. Strictly speaking, the number of advection-reaction evaluations is slightly less since :math:`\mathcal{F}_{\textrm{AR}}\left(t_0, \phi_0^k\right)` does not require re-evaluation in the corrector, and :math:`\mathcal{F}_{\textrm{AR}}\left(t_p,\phi_p^{K-1}\right)` does not need to be computed for the final iteration since the lagged quadrature is not further needed. Nonetheless, the computational work is quadratically increasing, but this is partially compensated by allowance of larger time steps since the :math:`\verb|SISDC|_K^K` has a stability limit of :math:`(K-1)\Delta t_{\textrm{cfl}}` rather than :math:`\Delta t_{\textrm{cfl}}` for uniformly spaced nodes. For comparison with the predictor :math:`\verb|SISDC|_K^1` which is a first order method, the work done for integration over :math:`(K-1)\Delta t_{\textrm{cfl}}` amounts to :math:`K-1` advection-reaction updates, :math:`K-1` diffusion updates, :math:`K-1` radiative transfer updates, and :math:`K` electric field updates. If we take the electric field updates as a reasonable metric for the computational work, the efficiency of the :math:`K` th order method over the first order method is about :math:`K` for integration over the same time interval, i.e. it increases linearly rather than quadratically. However, this estimate is only valid if we do not take accuracy into account. In practice, the predictor does not provide the same accuracy as the corrector over the same integration interval. A fair comparison of the extra computational work involved would require that the accuracy of the two methods be the same after integration over a time :math:`(K-1)\Delta t_{\textrm{cfl}}`, which will generally require more substeps for the first order method. While we do not further pursue this quantification in this paper, the pertinent point is that the extra computational work involved for tolerance-bound higher order discretizations increases sub-linearly rather than quadratically when compared to lower-order equivalents.

   We have implemented the SISDC algorithm in the ``imex_sdc`` class in :file:`physics/CdrPlasma/time_steppers/imex_sdc`.      

.. _Chap:CdrPlasmaPhysics:

CdrPlasmaPhysics
================

Overview
________

:ref:`Chap:CdrPlasmaPhysics` is an abstract class which represents the plasma physics for the CDR plasma module, i.e. it provides the coupling functions in :eq:`CdrPlasmaCoupling`.
The source code for the class resides in :file:`/Physics/CdrPlasma/CD_CdrPlasmaPhysics.H`.
Note that the entire class is an interface, whose implementations are used by the time integrators that advance the equations.

There are no default input parameters for :ref:`Chap:CdrPlasmaPhysics`, as users must generally implement their own kinetics.
The class exists solely for providing the integrators with the necessary fundamentals for filling solvers with the correct quantities at the same time, for example filling source terms and drift velocities.

A successful implementation of :ref:`Chap:CdrPlasmaPhysics` has the following:

#. Instantiated a list of :ref:`Chap:CdrSpecies`.
   These become :ref:`Chap:CDR` solvers and contain initial conditions and basic transport settings for the convection-diffusion-reaction solvers.
  
#. Instantiated a list :ref:`Chap:RtSpecies`.
   These become :ref:`Chap:RadiativeTransfer` solvers and contain metadata for the radiative transport solvers.
  
#. Implemented the core functionality that couple the solvers together. 

``chombo-discharge`` automatically allocates the specified number of convection-diffusion-reaction and radiative transport solvers from the list of species the is intantiated.
For information on how to interface into the CDR solvers, see :ref:`Chap:CdrSpecies`.
Likewise, see :ref:`Chap:RtSpecies` for how to interface into the RTE solvers.

Implementation of the core functionality is comparatively straightforward, but can lead to boilerplate code.
For this reason we also provide an implementation layer :ref:`Chap:CdrPlasmaJSON` that provides a plug-and-play interface for specifying the plasma physics by using a JSON schema for description the physics.

Complete API
____________

The full API for the ``CdrPlasmaPhysics`` class is given below:

.. literalinclude:: ../../../../Physics/CdrPlasma/CD_CdrPlasmaPhysics.H
   :language: c++
   :dedent: 4
   :lines: 28-263

.. _Chap:CdrPlasmaJSON:

JSON interface
==============

Since implementations of :ref:`Chap:CdrPlasmaPhysics` are usually boilerplate, we provide a class ``CdrPlasmaJSON`` which can initialize and parse various types of initial conditions and reactions from a JSON input file.
This class is defined in ``$DISCHARGE_HOME/Physics/PlasmaModels/CdrPlasmaJSON``.

``CdrPlasmaJSON`` is a full implementation of ``CdrPlasmaPhysics`` which supports the definition of various species (neutral, plasma species, and photons) and methods of coupling them.
We expect that ``CdrPlasmaJSON`` provides the simplest method of setting up a new plasma model.
It is also comparatively straightforward to extend the class with further required functionality.

In the JSON interface, the radiative transfer solvers always solve for the number of photons that lead to photoionization events.
This means that the interpretation of :math:`\Psi` is the number of photoionization events during the previous time step.
This is true for both continuum and discrete radiative transfer models. 

Usage
-----

To use this plasma model, use ``-physics CdrPlasmaJSON`` when setting up a new plasma problem (see :ref:`Chap:CdrPlasmaNewProblem`).
When ``CdrPlasmaJSON`` is instantiated, the constructor will parse species, reactions, initial conditions, and boundary conditions from a JSON file that the user provides.
In addition, users can parse transport data or reaction rates from tabulated ASCII files that they provide.

To specify the input plasma kinetics file, include

Specifying input file
---------------------

``CdrPlasmaJSON`` will read a JSON file specified by the input variable ``CdrPlasmaJSON.chemistry_file``.

Discrete photons
----------------

There are two approaches when using discrete photons, and both rely on the user setting up the application with the Monte Carlo photon solver (rather than continuum solvers).
For an introduction to the particle radiative transfer solver, see :ref:`Chap:MonteCarloRTE`.

The user must use one of the following:

* Set the following class options:

   .. code-block:: text
		   
      CdrPlasmaJSON.discrete_photons = true
      
      McPhoto.photon_generation = deterministic
      McPhoto.source_type       = number

   When specifying ``CdrPlasmaJSON.discrete_photons = true``, ``CdrPlasmaJSON`` will do a Poisson sampling of the number of photons that are generated in each cell and put this in the radiative transfer solvers' source terms.
   This means that the radiative transfer solver source terms *contain the physical number of photons generated in one time step*. 
   To turn off sampling inside the radiative transfer solver, we specify ``McPhoto.photon_generation = stochastic`` and set ``McPhoto.source_type = number`` to let the solver know that the source contains the number of physical photons. 

* Alternatively, set the following class options:
   
   .. code-block:: text
		   
      CdrPlasmaJSON.discrete_photons = false
      
      McPhoto.photon_generation = stochastic
      McPhoto.source_type       = volume_rate
      
   In this case the ``CdrPlasmaJSON`` class will fill the solver source terms with the volumetric rate, i.e. the number of photons produced per unit volume and time.
   When ``McPhoto`` generates the photons it will compute the number of photons generated in a cell through Poisson sampling :math:`n = P\left(S_\gamma\Delta V\Delta t\right)` where :math:`P` indicates a Poisson sampling operator.

Fundamentally, the two approaches differ only in where the Poisson sampling is performed.
With the first approach, plotting the radiative transfer solver source terms will show the number of physical photons generated.
In the second approach, the source terms will show the volume photo-generation rate. 

Gas law and neutral background
------------------------------

General functionality
_____________________

To include the gas law and neutral species, include a JSON object ``gas`` with the field ``law`` specified.
Currently, ``law`` can be either ``ideal``, ``troposphere``, or ``table``.

The purpose of the gas law is to set the temperature, pressure, and neutral density of the background gas.
In addition, we specify the neutral species that are used through the simulation.
These species are *not* stored on the mesh; we only store function pointers to their temperature, density, and pressure. 

It is also possible to include a field ``plot`` which will then include the temperature, pressure, and density in plot files. 

Ideal gas
_________

To specify an ideal gas law, specify ideal gas law as follows:

.. code-block:: json

   {"gas":
     {
       "law": "ideal",
       "temperature": 300,
       "pressure": 1
     }
    }

In this case the gas pressure and temperatures will be as indicated, and the gas number density will be computed as

.. math::

   \rho = \frac{p_0^\prime N_{\textrm{A}}}{RT_0},

where :math:`p^\prime` is the pressure converted to Pascals.

Note that the input temperature should be specified in Kelvin, and the input pressure in atmospheres. 

Troposphere
___________

It is also possible to specify the pressure, temperature, and density to be functions of tropospheric altitude.
In this case one must specify the extra fields

* ``molar mass`` For specifying the molar mass (in :math:`\textrm{g}\cdot\textrm{mol}^{-1}`) of the gas.  
* ``gravity`` Gravitational acceleration :math:`g`.
* ``lapse rate`` Temperature lapse rate :math:`L` in units of :math:`\textrm{K}/\textrm{m}`.

In this case the gas temperature pressure, and number density are computed as

.. math::

   T(h) = T_0 - Lh

.. math::

   p(h) = p_0\left((1 - \frac{Lh}{T_0}\right)^{\frac{g M}{RL}}

.. math::

   \rho(h) = \frac{p^\prime(h) N_{\textrm{A}}}{RT(h)}

For example, specification of tropospheric conditions can be included by

.. code-block:: json
		
   {"gas":
     {
       "law": "troposphere",
       "temperature": 300,
       "pressure": 1,
       "molar_mass": 28.97,
       "gravity": 9.81,
       "lapse_rate": 0.0065,
       "plot": true       
     }
   }

Tabulated
_________

To specify temperature, density, and pressure as function of altitude, set ``law`` to ``table`` and incldue the following fields:

* ``file`` For specifying which file we read the data from.
* ``height`` For specifying the column where the height is stored (in meters).
* ``temperature`` For specifying the column where the temperature (in Kelvin) is stored.
* ``pressure`` For specifying the column where the pressure (in Pascals) is stored.
* ``density`` For specifying the column where the density (in :math:`\textrm{kg}\cdot\textrm{m}^{-3}`) is stored.
* ``molar mass`` For specifying the molar mass (in :math:`\textrm{g}\cdot\textrm{mol}^{-1}`) of the gas.
* ``min height`` For setting the minimum altitude in the ``chombo-discharge`` internal table.
* ``max height`` For setting the minimum altitude in the ``chombo-discharge`` internal table.
* ``res height`` For setting the height resolution in the ``chombo-discharge`` internal table.

For example, assume that our file ``MyAtmosphere.dat`` contains the following data:

.. code-block:: text

   # z [m]              rho [kg/m^3]    T [K]           p [Pa]
   0.0000000E+00	1.2900000E+00	2.7210000E+02	1.0074046E+05
   1.0000000E+03	1.1500000E+00	2.6890000E+02	8.8751220E+04
   2.0000000E+03	1.0320000E+00	2.6360000E+02	7.8074784E+04
   3.0000000E+03	9.2860000E-01	2.5690000E+02	6.8466555E+04
   4.0000000E+03	8.3540000E-01	2.4960000E+02	5.9844569E+04

If we want to truncate this data to altitude :math:`z \ in[1000\,\textrm{m}, 3000\,\textrm{m}]` we specify:

.. code-block:: json

   {"gas":
     {
       "law": "table",
       "file": "ENMSIS_Atmosphere.dat",
       "molar mass": 28.97,
       "height": 0,
       "temperature": 2,
       "pressure": 3,
       "density": 1,
       "min height": 1000,
       "max height": 3000,
       "res height": 10
     }
   }

Neutral species background
__________________________

Neutral species are included by an array ``neutral species`` in the ``gas`` object. 
Each neutral species must have the fields

* ``name`` Species name
* ``molar fraction`` Molar fraction of the species.

If the molar fractions do not add up to one, they will be normalized.

.. warning::
   
   Neutral species are *not* tracked on the mesh.
   They are simply stored as functions that allow us to obtain the (spatially varying) density, temperature, and pressure for each neutral species.
   If a neutral species needs to be tracked on the mesh (through e.g. a convection-diffusion-reaction solver) it must be defined as a plasma species.
   See :ref:`Chap:PlasmaSpeciesJSON`. 

For example, a standard nitrogen-oxygen atmosphere will look like:

.. code-block:: json
   
   {"gas":
     {
       "law": "ideal",
       "temperature": 300,
       "pressure": 1,
       "plot": true,
       "neutral species":
       [
	 {
	   "name": "O2",
	   "molar_fraction": 0.2
	 },
	 {
	   "name": "N2",
	   "molar_fraction": 0.8
	 }
     ]
 }

.. _Chap:PlasmaSpeciesJSON:

Plasma species
--------------

The list of plasma species is included by an array ``plasma species``.
Each entry *must* have the entries

* ``name`` (string) For identifying the species name.
* ``Z`` (integer) Species charge number.
* ``mobile`` (true/false) Mobile species or not.
* ``diffusive`` (true/false) Diffusive species or not.

Optionally, the field ``initial data``, can be included for providing initial data to the species
Details are discussed further below.

For example, a minimum version would look like

.. code-block:: json

   {"plasma species":
     [
       {"name": "N2+", "Z":  1, "mobile": false, "diffusive": false},
       {"name": "O2+", "Z":  1, "mobile": false, "diffusive": false},
       {"name": "O2-", "Z": -1, "mobile": false, "diffusive": false}
     ]
   }

Initial data
____________

Initial data can be provided with

* Function based densities.
* Computational particles (deposited using a nearest-grid-point scheme).

Density functions
^^^^^^^^^^^^^^^^^

To provide initial data one include ``initial data`` for each species.
Currently, the following fields are supported:

* ``uniform`` For specifying a uniform background density.
  Simply the field ``uniform`` and a density (in units of :math:`m^{-3}`)
* ``gauss2`` for specifying Gaussian seeds :math:`n = n_0\exp\left(-\frac{\left(\mathbf{x}-\mathbf{x_0}\right)^2}{2R^2}\right)`.
  ``gauss2`` is an array where each array entry must contain

  * ``radius``, for specifying the radius :math:`R`: 
  * ``amplitude``, for specifying the amplitude :math:`n_0`. 
  * ``position``, for specifying the seed position :math:`\mathbf{x}`.
    
  The position must be a 2D/3D array.

* ``gauss2`` for specifying Gaussian seeds :math:`n = n_0\exp\left(-\frac{\left(\mathbf{x}-\mathbf{x_0}\right)^4}{2R^4}\right)`.
  ``gauss4`` is an array where each array entry must contain

  * ``radius``, for specifying the radius :math:`R`: 
  * ``amplitude``, for specifying the amplitude :math:`n_0`. 
  * ``position``, for specifying the seed position :math:`\mathbf{x}`.
    
  The position must be a 2D/3D array.  

* ``height profile`` For specifying a height profile along :math:`y` in 2D, and :math:`z` in 3D.
  To include it, prepare an ASCII files with at least two columns.
  The height (in meters) must be specified in one column and the density (in units of :math:`m^{-3}`) in another.
  Internally, this data is stored in a lookup table (see :ref:`Chap:LookupTable1D`). 
  Required fields are
  
  * ``file`` , for specifying the file.
  * ``height``, for specifying the column that stores the height.
  * ``density``, for specifying the column that stores the density.
  * ``min height``, for trimming data to a minimum height.
  * ``max height``, for trimming data to a maximum height.
  * ``res height``, for specifying the resolution height in the ``chombo-discharge`` lookup tables.

  In addition, height and density columns can be scaled in the internal tables by including

  * ``scale height`` for scaling the height data.
  * ``scale density`` for scaling the density data.

.. note::

   When multiple initial data fields are specified, ``chombo-discharge`` takes the superposition of all of them.

Initial particles
^^^^^^^^^^^^^^^^^

Initial particles can be included with the ``initial particles`` field.
The current implementation supports

* ``uniform`` For drawing initial particles randomly distributed inside a box.
  The user must specify the two corners ``lo corner`` and ``hi corner`` that indicate the spatial extents of the box, and the ``number`` of computational particles to draw.
  The weight is specified by a field ``weight``.
  For example:
  
  .. code-block:: json

   {"plasma species":
     [
       {
         "name": "e",
         "Z":  -1,
	 "mobile": true,
	 "diffusive": true,
	 "initial particles": {
	   "uniform": {
	     "lo corner": [0,0,0],
	     "hi corner": [1,1,1],
	     "number": 100,
	     "weight": 1.0
	   }
	 }
       }
     ]
   }


* ``sphere`` For drawing initial particles randomly distributed inside a sphere.
  Mandatory fields are

  * ``center`` for specifying the sphere center.
  * ``radius`` for specifying the sphere radius.
  * ``number`` for the number of computational particles.
  * ``weight`` for the initial particle weight.

  .. code-block:: json

   {"plasma species":
     [
       {
         "name": "e",
         "Z":  -1,
	 "mobile": true,
	 "diffusive": true,
	 "initial particles": {
	   "sphere": {
	     "center": [0,0,0],
	     "radius": 1.0,
	     "number": 100,
	     "weight": 1.0
	   }
	 }
       }
     ]
   }

* ``copy`` For using an already initialized particle distribution.
  The only mandatory fields is ``copy``, e.g.

  .. code-block:: json
		  
     {"plasma species":
     [
       {
         "name": "e",
         "Z":  -1,
	 "mobile": true,
	 "diffusive": true,
	 "initial particles": {
	   "sphere": {
	     "center": [0,0,0],
	     "radius": 1.0,
	     "number": 100,
	     "weight": 1.0
	   }
	 }
       },
       {
         "name": "O2+",
         "Z": 1,
	 "mobile": true,
	 "diffusive": true,
	 "initial particles": {
	    "copy": "e"
	 }
       }
     ]}

  This will copy the particles from the species ``e`` to the species ``O2+``.

  .. warning::

     The species one copies from must be defined *before* the species one copies *to*.


Complex example
^^^^^^^^^^^^^^^

For example, a species with complex initial data that combines density functions with initial particles can look like:

.. code-block:: json

   {"plasma species":
     [
       {
         "name": "N2+",
         "Z":  1,
	 "mobile": false,
	 "diffusive": false,
	 "initial data": {
	   "uniform": 1E10,
	   "gauss2" :
	     [
	       {
	          "radius": 100E-6,
		  "amplitude": 1E18,
		  "position": [0,0,0]
	       },
	       {
	          "radius": 200E-6,
		  "amplitude": 2E18,
		  "position": [1,0,0]
	       }
	     ],
	    "gauss4":
	      [
	        {
	          "radius": 300E-6,
		  "amplitude": 3E18,
		  "position": [0,1,0]
		},
		{
	          "radius": 400E-6,
		  "amplitude": 4E18,
		  "position": [0,0,1]
		}
	      ],
	    "height profile": {
 	      "file": "MyHeightProfile.dat",
 	      "height": 0,
	      "density": 1,
	      "min height": 0,
	      "max height": 100000,
	      "res height": 10,
	      "scale height": 100,
	      "scale density": 1E6   
	    }
	 },
	 "initial particles": {
	   "sphere": {
	     "center": [0,0,0],
	     "radius": 1.0,
	     "number": 100,
	     "weight": 1.0
	   }
	 }	 
       }
     ]
   }

.. _Chap:CdrPlasmaJSONMobility:

Mobilities
__________

If a species is specified as mobile, the mobility is set from a field ``mobility``, and the field ``lookup`` is used to specify the method for computing it. 
Currently supported are:

* Constant mobility.
* Function-based mobility, i.e. :math:`\mu = \mu(E,N)`.
* Tabulated mobility, i.e. :math:`\mu = \mu(E,N)`.

The cases are discussed below. 

**Constant mobility**

Setting ``lookup`` to ``constant`` lets the user set a constant mobility.
If setting a constant mobility, the field ``value`` is also required.
For example:

.. code-block:: json
		  
   {"plasma species":
     [
       {"name": "e", "Z":  -1, "mobile": true, "diffusive": false,
	"mobility": {
	  "lookup" : "constant",
	  "value": 0.05,
	 }
       }
     ]
   }

**Function-based mobility**

Setting ``lookup`` to ``function E/N``  lets the user set the mobility as a function of the reduced electric field.
When setting a function-based mobility, the field ``function`` is also required.

Supported functions are: 

* ``ABC``, in which case the mobility is computed as

   .. math::

      \mu(E) = A \frac{E^B}{N^C}.

   The fields ``A``, ``B``, and ``C`` must also be specified.
   For example:

   .. code-block:: json
		  
      {"plasma species":
        [
	  {"name": "e", "Z":  -1, "mobile": true, "diffusive": false,
	   "mobility": {
	     "lookup" : "function E/N",
	     "function": "ABC",
	     "A": 1,
	     "B": 1,
	     "C": 1	   
	    }
	  }
	]
      }

**Tabulated mobility**

Specifying ``lookup`` to ``table E/N`` lets the user set the mobility from a tabulated value of the reduced electric field.
BOLSIG-like files can be parsed by specifying the header which contains the tabulated data, and the columns that identify the reduced electric field and mobilities.
This data is then stored in a lookup table, see :ref:`Chap:LookupTable1D`.

For example:

.. code-block:: json

   {"plasma species":
     [
       {"name": "e", "Z":  -1, "mobile": true, "diffusive": false,
        "mobility": {
	  "lookup" : "table E/N",
	  "file": "transport_file.txt",
	  "header": "# Electron mobility (E/N, mu*N)",
	  "E/N ": 0,
	  "mu*N": 1,
	  "min E/N": 10,
	  "max E/N": 1000,
	  "points": 100,
	  "spacing": "exponential",
	  "dump": "MyMobilityTable.dat"
	 }
       }
     ]
   }

In the above, the fields have the following meaning:

* ``file`` The file where the data is found.
  The data must be stored in rows and columns.
* ``header``, the contents of the line preceding the table data.
* ``E/N``, the column that contains :math:`E/N`.
* ``mu*N``, the column that contains :math:`\mu\cdot E`.
* ``min E/N``, for trimming the data range.
* ``max E/N``, for trimming the data range.
* ``points``, for specifying the number of points in the lookup table.
* ``spacing``, for specifying how to regularize the table.
* ``dump``, an optional argument (useful for debugging) which will write the table to file. 

Note that the input file does *not* need regularly spaced or sorted data.
For performance reasons, the tables are always resampled, see :ref:`Chap:LookupTable1D`.

Diffusion coefficients
______________________

Setting the diffusion coefficient is done *exactly* in the same was as the mobility.
If a species is diffusive, one must include the field ``diffusion`` as well as ``lookup``.
For example, the JSON input for specifying a tabulated diffusion coefficient is done by

.. code-block:: json

   {"plasma species":
     [
       {"name": "e", "Z":  -1, "mobile": false, "true": false,
        "diffusion": {
	  "lookup" : "table E/N",
	  "file": "transport_file.txt",
	  "header": "# Electron diffusion coefficient (E/N, D*N)",
	  "E/N ": 0,
	  "D*N": 1,
	  "min E/N": 10,
	  "max E/N": 1000,
	  "points": 1000,
	  "spacing": "exponential"
	 }
       }
     ]
   }

Temperatures
____________

Plasma species temperatures can set by including a field ``temperature`` for the plasma species.

.. warning::
   
   If the ``temperature`` field is omitted, the species temperature will be set to the gas temperature.

**Constant temperature**


To set a constant temperature, include the field ``temperature`` and set ``lookup`` to constant and specify the temperature through the field ``value`` as follows:

.. code-block:: json

   {"plasma species":
     [
       {
         "name": "O2",
         "Z":  0,
	 "mobile": false,
         "true": false,
	 "temperature": {
	   "lookup": "constant",
  	   "value": 300
	 }
	}
     ]
   }

**Tabulated temperature**

To include a tabulated temperature :math:`T = T(E,N)`, set ``lookup`` to ``table E/N``.
The temperature is then computed as

.. math::

   T = \frac{2 \epsilon}{3k_{\textrm{B}}},

where :math:`\epsilon` is the energy and :math:`k_{\textrm{B}}` is the Boltzmann constant. 

The following fields are required:

* ``file`` for specifying which file the temperature is stored.
* ``header`` for specifying where in the file the temperature is stored.
* ``E/N`` for specifying in which column we find :math:`E/N`.
* ``eV`` for specifying in which column we find the species energy (in units of electron volts).
* ``min E/N`` for trimming the data range.
* ``max E/N`` for trimming the data range.
* ``points`` for setting the number of points in the lookup table.
* ``spacing`` for setting the grid point spacing type.
* ``dump`` for writing the final table to file.

For a further explanation to these fields, see :ref:`Chap:CdrPlasmaJSONMobility`.

A complete example is:

.. code-block:: json

   {"plasma species":
     [
       {
         "name": "e",
         "Z":  -1,
	 "mobile": true,
         "true": true,
	 "temperature": {
	   "lookup": "table E/N",
	   "file": "transport_data.txt",
	   "header": "# Electron mean energy (E/N, eV)",
	   "E/N": 0,
	   "eV": 1,
	   "min E/N": 10,
	   "max E/N": 1000,
	   "points": 1000,
	   "spacing": "exponential",
	   "dump": "MyTemperatureTable.dat"
	 }
	}
     ]
   }

Photon species
--------------

As for the plasma species, photon species (for including radiative transfer) are included by an array ``photon species``.
For each species, the required fields are

* ``name`` For setting the species name.
* ``kappa`` For specifying the absorption coefficient. 

Currently, ``kappa`` can be either

* ``constant`` Which lets the user set a constant absorption coefficient. 
* ``helmholtz`` Computes the absorption coefficient as
  
  .. math::

     \kappa = \frac{p_X\lambda}{\sqrt{3}}

  where :math:`\lambda` is a specified input parameter and :math:`p_X` is the partial pressure of some species :math:`X`.

* ``stochastic A`` which samples a random absorption coefficient as

  .. math::

     \kappa = K_1 \left(\frac{K_2}{K_1}\right)^{\frac{f-f1}{f2-f1}}.

  Here, :math:`f_1` and :math:`f_2` are frequency ranges, :math:`K_1` and :math:`K_2` are absorption coefficients, and :math:`f` is a stochastically sampled frequency.
  Note that this method is only sensible when using discrete photons.


**Constant absorption coefficients**

When specifying a constant absorption coefficient, one must include a field ``value`` as well.
For example:

.. code-block:: json

 {"photon species":
    [
      {
        "name": "UVPhoton",
        "kappa": "constant",
        "value": 1E4
      }
    ]
 }

**Helmholtz absorption coefficients**

The interface for the Helmholtz-based absorption coefficients are inspired by :cite:t:`Bourdon2007` approach for computing photoionization.
This method only makes sense if doing a Helmholtz-based reconstruction of the photoionization profile as a relation:

.. math::
   
   \left[\nabla^2 - \left(p_{\textrm{O}_2} \lambda\right)^2\right]S_\gamma = -\left(A p_{\textrm{O}_2}^2\frac{p_q}{p + p_q}\xi\nu\right)S_i,

where

* :math:`S_\gamma` is the number of photoionization events per unit volume and time. 
* :math:`A` is a model coefficient.
* :math:`\frac{p_q}{p + p_q}` is a quenching factor.
* :math:`\xi` is a photoionization efficiency.
* :math:`\nu` is a relative excitation efficiency.
* :math:`S_i` is the electron impact ionization source term.

Since the radiative transfer solver is based on the Eddington approximation, the Helmholtz reconstruction can be written as

.. math::
   
   \kappa \Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla \Psi\right) = \frac{\eta}{c}

where the absorption coefficient is set as

.. math::
     
   \kappa(\mathbf{x}) = \frac{p_{\textrm{O}_2}\lambda}{\sqrt{3}}.

The photogeneration source term is still

.. math::

   \eta = \frac{p_q}{p + p_q}\xi\nu S_i,

but the photoionization term is

.. math::

   S_\gamma = \frac{c A p_{\textrm{O}_2}}{\sqrt{3}\lambda}\Psi.

Note that the photoionization term is, in principle, *not* an Eddington approximation.
Rather, the Eddington-like equations occur here through an approximation of the exact integral solution to the radiative transfer problem. 
In the pure Eddington approximation, on the other hand, :math:`\Psi` represents the total number of ionizing photons per unit volume, and we would have :math:`S_\gamma = \frac{\Psi}{\Delta t}` where :math:`\Delta t` is the time step.

When specifying the ``kappa`` field as ``helmholtz``, the absorption coefficient is computed as

.. math::

   \kappa(\mathbf{x}) = \frac{p_X\left(\mathbf{x}\right)\lambda}{\sqrt{3}}

where :math:`p_X` is the partial pressure of a species :math:`X` and :math:`\lambda` is the same input parameter as in the Helmholtz reconstruction. 
These are specified through fields ``neutral`` and ``lambda`` as follows:

.. code-block:: json

 {"photon species":
    [
      {
        "name": "UVPhoton",
        "kappa": "helmholtz",
        "lambda": 0.0415,
	"neutral": "O2"
      }
    ]
 }

This input will set :math:`\kappa\left(\mathbf{x}\right) = \frac{p_{\textrm{O}_2}\left(\mathbf{x}\right)\lambda}{\sqrt{3}}`.

.. note::
   
   The source term :math:`\eta` is specified when specifying the plasma reactions, see :ref:`Chap:CdrPlasmaReactionsJSON`.

**Stochastic sampling**

Setting the ``kappa`` field to ``stochastic A`` will stochastically sample the absorption length from

.. math::

   \kappa = K_1 \left(\frac{K_2}{K_1}\right)^{\frac{f-f1}{f2-f1}}.

where :math:`K_1 = p_X\chi_{\textrm{min}}`, :math:`K_1 = p_X\chi_{\textrm{max}}`, and :math:`f_1` and :math:`f_2` are frequency ranges.
Like above, :math:`p_X` is the partial pressure of some species :math:`X`.
Note that all input parameters are given in SI units. 

Stochastic sampling of the absorption length only makes sense when using discrete photons -- this particular method is inspired by the method in :cite:t:`Chanrion2008`.
For example:

.. code-block:: json

 {"photon species":
    [
      {
        "name": "UVPhoton",
        "kappa": "stochastic A",
        "neutral": "O2",	 
        "f1":   2.925E15,
        "f2":   3.059E15,
        "chi min": 2.625E-2,
        "chi max": 1.5
      }
    ]
 }  


.. _Chap:CdrPlasmaReactionsJSON:

Plasma reactions
----------------

Plasma reactions are reactions between charged and neutral species and are written in the form

.. math::

   A + B + \ldots \rightarrow C + D + \ldots.

Importantly, the left hand side of the reaction can only consist of charged or neutral species.
It is not permitted to put a photon species on the left hand side of these reactions; photo-ionization is handled separately by another set of reaction types (see :ref:`Chap:PhotoReactionsJSON`).
However, photon species *can* appear on the left hand side of the equation. 

When specifying reactions in this form, the reaction rate is computed as

.. math::

   R = k n_A n_B\ldots

When computing the source term for some species :math:`X`, we subtract :math:`R` for each time :math:`X` appears on the left hand side of the reaction and add :math:`R` for each time :math:`X` appears on the right-hand side of the reaction.

Specifying reactions
____________________

Reactions of the above type are handled by a JSON array ``plasma reactions``, with required fields:

* ``reaction`` (string) containing the reaction process.
* ``lookup`` (string) for determining how to compute the reaction rate. 

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + O2 -> e + e + O2+",
	 "lookup": "constant",
	 "rate": 1E-30
       }
     ]
   }

This adds a reaction :math:`\textrm{e} + \textrm{O}_2 \rightarrow \textrm{e} + \textrm{e} + \textrm{O}_2^+` to the reaction set.
We compute

.. math::

   R = kn_{\textrm{e}}n_{\textrm{O}_2^+}

and set

.. math::

   S_{\textrm{e}} = S_{\textrm{O}_2^+} = R.
   
Some caveats when setting the reaction string are:

* Whitespace are separators.
  For example, ``O2+e`` will be interpreted as a species with string identifier ``O2+e``, but ``O2 + e`` will interpreted as a reaction between ``O2`` and ``e``.
* The reaction string *must* contain a left and right hand side separated by ``->``.
  An error will be thrown if this symbol can not be found. 
* The left-hand must consist *only* of neutral or plasma species.
  If the left-hand side consists of species that are not neutral or plasma species, an error will be thrown. 
* The right-hand side can consist of either neutral, plasma species, or photon species.
  Otherwise, an error will be thrown.
* The reaction string will be checked for charge conservation.


Note that if a reaction involves a right-hand side that is not otherwise tracked, the user should omit the species from the right-hand side altogether.
For example, if we have a model which tracks the species :math:`e` and :math:`\textrm{O}_2^+` but we want to include the dissociative recombination reaction :math:`e + \textrm{O}_2^+ \rightarrow O + O`, this reaction should be added to the reaction with an empty right-hand side:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + O2 -> e + e + O2+",
	 "lookup": "constant",
	 "rate": 1E-30
       },     
       {
         "reaction": "e + O2+ -> ",
	 "lookup": "constant",
	 "rate": 1E-30
       }
     ]
   }

.. _CdrPlasma:Wildcards:

Wildcards
_________

Reaction specifiers may include the wildcard ``@`` which is a placeholder for another species.
The wildcards must be specified by including a JSON array ``@`` of the species that the wildcard is replaced by.
For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "N2+ + N2 + @ -> N4+ + @",
	 "@": ["N2", "O2"],
	 "lookup": "constant",
	 "rate": 1E-30
       }
     ]
   }
   
The above code will add two reactions to the reaction set: :math:`N_2 + N_2 + N_2 \rightarrow N_4^+ + N_2` and :math:`N_2 + N_2 + \textrm{O}_2 \rightarrow N_4^+ + \textrm{O}_2`.
It is not possible to set different reaction rates for the two reactions. 


Specifying reaction rates
_________________________

Constant reaction rates
^^^^^^^^^^^^^^^^^^^^^^^

To set a constant reaction rate for a reaction, set the field ``lookup`` to ``"constant"`` and specify the rate.
For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + O2 -> e + e + O2+",
	 "lookup": "constant",
	 "rate": 1E-30
       }
     ]
   }

Single-temperature rates
^^^^^^^^^^^^^^^^^^^^^^^^

* ``functionT A``
  To set a rate dependent on a single species temperature in the form :math:`k(T) = c_1T^{c_2}`, set ``lookup`` to ``functionT A``.
  The user must specify the species from which we compute the temperature :math:`T` by including a field ``T``.
  The constants :math:`c_1` and :math:`c_2` must also be included.

  For example, in order to add a reaction :math:`e + \textrm{O}_2 \rightarrow \varnothing` with rate :math:`k = 1.138\times 10^{-11}T_{\textrm{e}}^{-0.7}` we can add the following:

    .. code-block:: json

     {"plasma reactions":
       [
         {
           "reaction": "e + M+ ->",
	   "lookup": "functionT A",
	   "T": "e",
	   "c1": 1.138
	   "c2": -0.7
         }
       ]
     } 
  

Two-temperature rates
^^^^^^^^^^^^^^^^^^^^^

* ``functionT1T2 A``
  To set a rate dependent on two species temperature in the form :math:`k(T_1, T_2) = c_1\left(T_1/T_2\right)^{c_2}`, set ``lookup`` to ``functionT1T2 A``.
  The user must specify which temperatures are involved by specifying the fields ``T1``, ``T2``, as well as the constants through fields ``c1`` and ``c2``.
  For example, to include the reaction :math:`e + \textrm{O}_2 + \textrm{O}_2 \rightarrow \textrm{O}_2^- + O2` in the set, with this reaction having a rate

  .. math::

     k = 2.4\times 10^{-41}\left(\frac{T_{\textrm{O}_2}}{T_e}\right),

  we add the following:

  .. code-block:: json

     {"plasma reactions":
       [
         {
           "reaction": "e + O2 + O2 -> O2- + O2",
	   "lookup": "functionT1T2 A",
	   "T1": "O2",
	   "T2": "e",
	   "c1": 2.41E-41,
	   "c2": 1
         }
       ]
     }

Townsend ionization and attachment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To set standard Townsend ionization and attachment reactions, set ``lookup`` to ``alpha*v`` and ``eta*v``, respectively.
This will compute the rate constant :math:`k = \alpha \left|\mathbf{v}\right|` where :math:`\mathbf{v}` is the drift velocity of some species.
To specify the species one includes the field ``species``.

For example, to include the reactions :math:`e \rightarrow e + e + M^+` and :math:`e \rightarrow M^-` one can specify the reactions as

.. code-block:: json

     {"plasma reactions":
       [
         {
           "reaction": "e -> e + e + M+",
	   "lookup": "alpha*v",
	   "species": "e"
         },
         {
           "reaction": "e -> M-",
	   "lookup": "eta*v",
	   "species": "e"
         }	 
       ]
     }




Tabulated rates
^^^^^^^^^^^^^^^

To set a tabulated rate with :math:`k = k(E,N)`, set the field ``lookup`` to ``table E/N`` and specify the file, header, and data format to be used.
For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + O2 -> e + e + O2+",
         "lookup": "table E/N",
	 "file": "transport_file.txt",
	 "header": "# O2 ionization (E/N, rate/N)",
	 "E/N ": 0,
	 "rate/N": 1,
	 "min E/N": 10,
	 "max E/N": 1000,
	 "spacing": "exponential",
	 "points": 1000,
	 "plot": true,
	 "dump": "O2_ionization.dat"
       }
     ]
   }

The ``file`` field specifies which field to read the reaction rate from, while ``header`` indicates where in the file the reaction rate is found.
The file parser will read the files below the header line until it reaches an empty line.
The fields ``E/N`` and ``rate/N`` indicate the columns where the reduced electric field and reaction rates are stored.

The final fields ``min E/N``, ``max E/N``, and ``points`` are formatting fields that trim the range of the data input and organizes the data along a table with ``points`` entries.
As with the mobilities (see :ref:`Chap:CdrPlasmaJSONMobility`), the ``spacing`` argument determines whether or not the internal interpolation table uses uniform or exponential grid point spacing.
Finally, the ``dump`` argument will tell ``chombo-discharge`` to dump the table to file, which is useful for debugging or quality assurance of the tabulated data.

Modifying reactions
___________________

Collisional quenching
^^^^^^^^^^^^^^^^^^^^^

To quench a reaction, include a field ``qenching_pressure`` and specify the *quenching pressure* (in atmospheres).
When computing reaction rates, the rate for the reaction will be modified as

.. math::

   k \rightarrow k\frac{p_q}{p_q + p}

where :math:`p^q` is the quenching pressure and :math:`p = p(\mathbf{x})` is the gas pressure.

.. important::

   The quenching pressure should be specified in Pascal. 

For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + N2 -> e + N2 + Y",
         "lookup": "table E/N",
	 "file": "transport_file.txt",
	 "header": "# N2 ionization (E/N, rate/N)",
	 "E/N ": 0,
	 "rate/N": 1,
	 "min E/N": 10,
	 "max E/N": 1000,
	 "points": 1000,
	 "spacing": "exponential",
	 "quenching pressure": 4000
       }
     ]
   }

Reaction efficiencies
^^^^^^^^^^^^^^^^^^^^^

To modify a reaction efficiency, include a field ``efficiency`` and specify it.
This will modify the reaction rate as

.. math::

   k \rightarrow \nu k

where :math:`\nu` is the reaction efficiency.
For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + N2 -> e + N2 + Y",
         "lookup": "table E/N",
	 "file": "transport_file.txt",
	 "header": "# N2 ionization (E/N, rate/N)",
	 "E/N ": 0,
	 "rate/N": 1,
	 "min E/N": 10,
	 "max E/N": 1000,
	 "points": 1000,
	 "spacing": "exponential",
	 "efficiency": 0.6
       }
     ]
   }

Scaling reactions
^^^^^^^^^^^^^^^^^

Reactions can be scaled by including a ``scale`` argument to the reaction.
This works exactly like the ``efficiency`` field outlined above.

Energy correction
^^^^^^^^^^^^^^^^^

Occasionally, it can be necessary to incorporate an energy correction to models, accounting e.g. for electron energy loss near strong gradients.
The JSON interface supports the correction in :cite:t:`Soloviev2009`.
To use it, include an (optional) field ``soloviev`` and specify ``correction`` and ``species``.
For example:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + N2 -> e + N2 + Y",
         "lookup": "table E/N",
	 "file": "transport_file.txt",
	 "header": "# N2 ionization (E/N, rate/N)",
	 "E/N ": 0,
	 "rate/N": 1,
	 "min E/N": 10,
	 "max E/N": 1000,
	 "points": 1000,
	 "spacing": "exponential",
	 "efficiency": 0.6,
	 "soloviev": {
	   "correction": true,
	   "species": "e"
	 }
       }
     ]
   }

When this energy correction is enabled, the rate coefficient is modified as

.. math::

   k \rightarrow k\left(1 + \frac{\mathbf{E}\cdot D_s\nabla n_s}{\mu_s n_s E^2}\right),

where :math:`s` is the species specified in the ``soloviev`` field, :math:`n_s` is the density and :math:`D_s` and :math:`\mu_s` are diffusion and mobility coefficients.
We point out that the correction factor is restricted such that the reaction rate is always non-negative. 
Note that this correction makes sense when rates are dependent only on the electric field, see :cite:t:`Soloviev2009`.

.. note::

   When using the energy correction, the specifies species must be both mobile and diffusive.

Plotting reactions
__________________

It is possible to have ``CdrPlasmaJSON`` include the reaction rates in the HDF5 output files by including a field ``plot`` as follows:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + O2 -> e + e + O2+",
	 "plot": true,
	 "lookup": "constant",
	 "rate": 1E-30,
       }
     ]
   }

Plotting the reaction rate can be useful for debugging or analysis.
Note that it is, by extension, also possible to add useful data to the I/O files from reactions that otherwise do not contribute to the discharge evolution.
For example, if we know the rate :math:`k` for excitation of nitrogen to a specific excited state, but do not otherwise care about tracking the excited state, we can add the reaction as follows:

.. code-block:: json

   {"plasma reactions":
     [
       {
         "reaction": "e + N2 -> e + N2",
	 "plot": true,	 
	 "lookup": "constant",
	 "rate": 1E-30,
       }
     ]
   }

This reaction is a dud in terms of the discharge evolution (the left and right hand sides are the same), but it can be useful for plotting the excitation rate. 

.. note:: This functionality should be used with care because each reaction increases the I/O load.

Warnings and caveats
____________________

Note that the JSON interface *always* computes reactions as if they were specified by the deterministic reaction rate equation

.. math::

   \partial_t n_i = \sum_r k_r n_j n_k n_l\ldots,
   
where the fluid source term for any reaction :math:`r` is :math:`S_r = k_r n_j n_k n_l\ldots`.
Caution should always be exercised when defining a reaction set.

Higher-order reactions
^^^^^^^^^^^^^^^^^^^^^^

Usually, many rate coefficients depend on the output of other software (e.g., BOLSIG+) and the scaling of rate coefficients is not immediately obvious.
This is particularly the case for three-body reactions with BOLSIG+ that may require scaling before running the Boltzmann solver (by scaling the input cross sections), or after running the Boltzmann solver, in which case the rate coefficients themselves might require scaling. In any case the user should investigate the cross-section file that BOLSIG+ uses, and figure out the required scaling. 

.. important::
   
   For two-body reactions, e.g. :math:`A + B\rightarrow \varnothing` the rate coefficient must be specified in units of :math:`\textrm{m}^3\textrm{s}^{-1}`, while for three-body reactions :math:`A + B + C\rightarrow\varnothing` the rate coefficient must have units of :math:`\textrm{m}^6\textrm{s}^{-1}`.

   For three-body reactions the units given by BOLSIG+ in the output file may or may not be incorrect (depending on whether or not the user scaled the cross sections).

Townsend coefficients
^^^^^^^^^^^^^^^^^^^^^

Townsend coefficients are not fundamentally required for specifying the reactions, but as with the higher-order reactions some of the output rates for three-body reactions might be inconsistently represented in the BOLSIG+ output files.
For example, some care might be required when using the Townsend attachment coefficient for air when the reaction :math:`\textrm{e} + \textrm{O}_2 + \textrm{O}_2\rightarrow \textrm{O}_2^- + \textrm{O}_2` is included because the rate constant might require proper scaling after running the Boltzmann solver, but this scaling is invisible to the BOLSIG+'s calculation of the attachment coefficient :math:`\eta/N`.

.. warning::

   The JSON interface *does not guard* against inconsistencies in the user-provided chemistry, and provision of inconsistent :math:`\eta/N` and attachment reaction rates are quite possible.

.. _Chap:PhotoReactionsJSON:

Photo-reactions
---------------

Photo-reactions are reactions between charged/neutral and photons in the form

.. math::

   A + B + \gamma + \ldots \rightarrow C + D + \ldots.

where species :math:`A, B, \ldots` are charged and neutral species and :math:`\gamma` is a photon.
The left hand side can contain only *one* photon species, and the right-hand side can not contain a photon species.
In other words, two-photon absorption is not supported, and photons that are absorbed on the mesh cannot become new photons.
This is not a fundamental limitation, but a restriction imposed by the JSON interface. 

Specifying reactions
____________________

Reactions of the above type are handled by a JSON array ``photo reactions``, with required fields:

* ``reaction`` (string) containing the reaction process.
* ``lookup`` (string) for determining how to compute the reaction rate.

For example:

.. code-block:: json

   {"photo reactions":
     [
       {"reaction": "Y + O2 -> e + O2+"},
     ]
   }

The rules for specifying reaction strings are the same as for the plasma reactions, see :ref:`Chap:CdrPlasmaReactionsJSON`.
Wildcards also apply, see :ref:`CdrPlasma:Wildcards`.

Default behavior
________________

Since the radiative transfer solvers solve for the number of ionizing photons, the CDR solver source terms are incremented by

.. math::

   S \rightarrow S + \frac{\Psi}{\Delta t}.

where :math:`\Psi` is the number of ionizing photons per unit volume (i.e., the solution :math:`\Psi`). 

Helmholtz reconstruction
________________________

When performing a Helmholtz reconstruction the photoionization source term is

.. math::

   S = \frac{c A p_{\textrm{O}_2}}{\sqrt{3}\lambda}\Psi.

To modify the source term for consistency with Helmholtz reconstruction specify the field ``helmholtz`` with variables

* ``A``. the :math:`A` coefficient.
* ``lambda``. the :math:`\lambda` coefficient.
  This value will also be specified in the photon species, but it is not retrieved automatically.
* ``neutral``. The neutral species for which we obtain the partial pressure. 

For example:

.. code-block:: json

   {"photo reactions":
     [
       {
         "reaction": "Y + O2 -> e + O2+",
	 "helmholtz": {
	   "A": 1.1E-4,
	   "lambda": 0.0415,
	   "neutral": "O2"
	 }
       }
     ]
   }

Scaling reactions
_________________

Photo-reactions can be scaled by including a ``scale`` argument.
For example, to completely turn off the photoreaction above:

.. code-block:: json

   {"photo reactions":
     [
       {
         "reaction": "Y + O2 -> e + O2+",
	 "helmholtz": {
	   "A": 1.1E-4,
	   "lambda": 0.0415,
	   "neutral": "O2"
	 },
	 "scale": 0.0
       }
     ]
   }

EB boundary conditions
----------------------

Boundary conditions on the embedded boundary are included by the fields

* ``electrode reactions``, for specifying secondary emission on electrodes.
* ``dielectric reactions``, for specifying secondary emission on dielectrics.

To include secondary emission, the user must specify a reaction string in the form :math:`A \rightarrow B`, and also include an emission rate.
Currently, we only support constant emission rates (i.e., secondary emission coefficients).
This is likely to change in the future. 

The following points furthermore apply:

* By default, standard outflow boundary conditions.
  When ``electrode reactions`` or ``dielectric reactions`` are specified, the user only controls the *inflow* back into the domain.
* Wildcards can appear on the left hand side of the reaction.  
* If one specifies :math:`A + B \rightarrow C` for a surface reaction, this is the same as specifying two reactions :math:`A \rightarrow C` and :math:`B\rightarrow C`.
  The same emission coefficient will be used for both reactions.
* Both photon species and plasma species can appear on the left hand side of the reaction.
* Photon species can not appear on the right-hand side of the reaction; we do not include surface sources for photoionization.
* To scale reactions, include a modifier ``scale``.  

For example, the following specification will set secondary emission efficiencies to :math:`10^{-3}`:

.. code-block:: json

 {"electrode reactions":
   [
     { "reaction": "@ -> e",
       "@": ["N2+", "O2+", "N4+", "O4+", "O2+N2"],
       "lookup": "constant",
       "value": 1E-4
     }
   ],
  "dielectric reactions":
   [
     { "reaction": "@ -> e",
       "@": ["N2+", "O2+", "N4+", "O4+", "O2+N2"],
       "lookup": "constant",
       "value": 1E-3
     }
   ] 		

Adaptive mesh refinement
========================

The AMR functionality for ``CdrPlasmaStepper`` is implemented by subclassing ``CellTagger`` through a series of intermediate classes.
Currently, we mostly use the ``CdrPlasmaStreamerTagger`` class for tagging cells based on the evaluation of the Townsend ionization coefficient as

.. math::

   \left(\alpha-\eta\right)\Delta x > \epsilon,

where :math:`\epsilon` is a refinement threshold.
A similar threshold is used when coarsening grid cells.

.. tip::

   See `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classPhysics_1_1CdrPlasma_1_1CdrPlasmaStreamerTagger.html>`_ for the C++ API for ``CdrPlasmaStreamerTagger``. 

.. _Chap:CdrPlasmaNewProblem:

Setting up a new problem
========================

New problems that use the ``CdrPlasma`` physics model are best set up by using the Python tools provided with the module.
A full description is available in the ``README.md`` file contained in the folder:

.. literalinclude:: ../../../../Physics/CdrPlasma/README.md
   :language: markdown
              
To see the list of available options type

.. code-block:: bash

   cd $DISCHARGE_HOME/Physics/CdrPlasma
   ./setup.py --help
