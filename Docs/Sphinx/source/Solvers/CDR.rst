.. _Chap:CDR:

Convection-Diffusion-Reaction
=============================

Here, we discuss the discretization of the equation 

.. math::
   
   \frac{\partial \phi}{\partial t} + \nabla\cdot\left(\mathbf{v} \phi - D\nabla \phi\right) = S.

We assume that :math:`\phi` is discretized by cell-centered averages (note that cell centers may lie inside solid boundaries), and use finite volume methods to construct fluxes in a cut-cells and regular cells.
Here, :math:`\mathbf{v}` indicates a drift velocity and :math:`D` is the diffusion coefficient. 

.. note::
   
   Using cell-centered versions :math:`\phi` might be problematic for some models since the state is extended outside the valid region.
   Models might have to recenter the state in order compute e.g. physically meaningful reaction terms in cut-cells.

Source code for the convection-diffusion-reaction solvers reside in :file:`$DISCHARGE_HOME/Source/ConvectionDiffusionReaction` and the full API can be found at `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classCdrSolver.html>`_.

.. _Chap:CdrSolver:   

CdrSolver
---------

The ``CdrSolver`` class contains the interface for solving advection-diffusion-reaction problems.
``CdrSolver`` is an abstract class and does not contain any specific advective or diffusive discretization (these are added by implementation classes).

Currently, the implementation layers consist of the following:

#. :ref:`Chap:CdrMultigrid`, which inherits from ``CdrSolver`` and adds a second order accurate discretization for the diffusion operator.
#. :ref:`Chap:CdrCTU` and :ref:`Chap:CdrGodunov` which inherit from ``CdrMultigrid`` and add a second order accurate spatial discretization for the advection operator. 

Currently, we mostly use the ``CTU`` class which contains a second order accurate discretization with slope limiters.
``CdrGodunov`` is a similar operator, but the advection code for this is distributed by the ``Chombo`` team.
The C++ API for these classes can be obtained from `<https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classCdrSolver.html>`_.

The advance methods for ``CdrSolver`` are encapsulated by the following functions:

.. code-block:: c++

   // For advancing the diffusion equation with an implicit Euler method. 
   void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const EBAMRCellData& a_source, const Real a_dt) = 0;

   // For advancing the diffusion equation with a Crank-Nicholson method.
   void advanceCrankNicholson(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const EBAMRCellData& a_source, const Real a_dt) = 0;

   // For computing div(v*phi - D*grad(phi)).
   void computeDivJ(EBAMRCellData& a_divJ, EBAMRCellData& a_phi, const Real a_dt, const bool a_conservativeOnly, const bool a_ebFlux, const bool a_domainFlux) = 0;

   // For computing div(v*phi).
   void computeDivF(EBAMRCellData& a_divJ, EBAMRCellData& a_phi, const Real a_dt, const bool a_conservativeOnly, const bool a_ebFlux, const bool a_domainFlux) = 0;

   // For computing div(D*grad(phi)).
   void computeDivD(EBAMRCellData& a_divJ, EBAMRCellData& a_phi, const Real a_dt, const bool a_conservativeOnly, const bool a_ebFlux, const bool a_domainFlux) = 0;

.. _Chap:CdrDetails:   

Discretization details
----------------------

.. _Chap:ExplicitDivergence:   

Explicit divergences and redistribution
_______________________________________

Computing explicit divergences for equations like

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\mathbf{G} = 0

is problematic because of the arbitarily small volume fractions of cut cells.
In general, we seek a method-of-lines update :math:`\phi^{k+1} = \phi^k - \Delta t \left[\nabla\cdot \mathbf{G}^k\right]` where :math:`\left[\nabla\cdot\mathbf{G}\right]` is a stable numerical approximation based on some finite volume approximation.

Pure finite volume methods use

.. math::
   \phi^{k+1} = \phi^k - \frac{\Delta t}{\kappa \Delta x^{\textrm{DIM}}}\int_V\nabla\cdot\mathbf{G}dV,
   :label: conservativeUpdate
   
where :math:`\kappa` is the volume fraction of a grid cell, :math:`\textrm{DIM}` is the spatial dimension and the volume integral is written as discretized surface integral
   
.. math::
   \int_V\nabla\cdot\mathbf{G}dV =\sum_{f\in f(V)}\left(\mathbf{G}_f\cdot \mathbf{n}_f\right)\alpha_f\Delta x^{\textrm{DIM} -1}.
   
The sum runs over all cell edges (faces in 3D) of the cell where :math:`G_f` is the flux on the edge centroid and :math:`\alpha_f` is the edge (face) aperture.

.. figure:: /_static/figures/CutCell.png
   :width: 40%
   :align: center

   Location of centroid fluxes for cut cells. 

However, taking :math:`[\nabla\cdot\mathbf{G}^k]` to be this sum leads to a time step constraint proportional to :math:`\kappa`, which can be arbitrarily small.
This leads to an unacceptable time step constraint for :eq:`conservativeUpdate`.
We use the Chombo approach and expand the range of influence of the cut cells in order to stabilize the discretization and allow the use of a normal time step constraint.
First, we compute the conservative divergence

.. math::
  \kappa_{\mathbf{i}} D_\mathbf{i}^c =  \sum_f G_f\alpha_f\Delta x^{\textrm{DIM} -1},

where :math:`G_f = \mathbf{G}_f\cdot \mathbf{n}_f`. Next, we compute a non-conservative divergence :math:`D_{\mathbf{i}}^{nc}`

.. math::
   D_\mathbf{i}^{nc} =  \frac{\sum_{\mathbf{j}\in{N}\left(\mathbf{i}\right)}\kappa_{\mathbf{j}}D_\mathbf{i}^c}{\sum_{\mathbf{j}\in{N}\left(\mathbf{i}\right)}\kappa_{\mathbf{j}}}

where :math:`N(\mathbf{i})` indicates some neighborhood of cells around cell :math:`\mathbf{i}`. Next, we compute a hybridization of the divergences, 

.. math::
  D_{\mathbf{i}}^H = \kappa_{\mathbf{i}} D_{\mathbf{i}}^c + (1-\kappa_{\mathbf{i}})D_{\mathbf{i}}^{nc},

and perform an intermediate update
  
.. math::
   \phi_{\mathbf{i}}^{k+1} = \phi_{\mathbf{i}}^k - \Delta tD_{\mathbf{i}}^H.
   
The hybrid divergence update fails to conserve mass by an amount :math:`\delta M_{\mathbf{i}} = \kappa_{\mathbf{i}}\left(1-\kappa_{\mathbf{i}}\right)\left(D_{\mathbf{i}}^c - D_{\mathbf{i}}^{nc}\right)`.
In order to main overall conservation, the excess mass is redistributed into neighboring grid cells.
Let :math:`\delta M_{\mathbf{i}, \mathbf{j}}` be the redistributed mass from :math:`\mathbf{j}` to :math:`\mathbf{i}` where
   
.. math::
   \delta M_{\mathbf{i}} = \sum_{\mathbf{j} \in N(\mathbf{i})}\delta M_{\mathbf{i}, \mathbf{i}}.

This mass is used as a local correction in the vicinity of the cut cells, i.e.
   
.. math::
   \phi_{\mathbf{i}}^{k+1} \rightarrow \phi_{\mathbf{i}}^{k+1} + \delta M_{\mathbf{j}\in N(\mathbf{i}), \mathbf{i}},

where :math:`\delta M_{\mathbf{j}\in N(\mathbf{i}), \mathbf{i}}` is the total mass redistributed to cell :math:`\mathbf{i}` from the other cells.
After these steps, we define
   
.. math::
   \left[\nabla\cdot\mathbf{G}^k\right]_{\mathbf{i}} \equiv \frac{1}{\Delta t}\left(\phi_{\mathbf{i}}^{k+1} - \phi_{\mathbf{i}}^k\right)

Numerically, the above steps for computing a conservative divergence of a one-component flux :math:`\mathbf{G}` are implemented in the convection-diffusion-reaction solvers, which also respects boundary conditions (e.g. charge injection).
The user will need to call the function

.. code-block:: c++
		
   virtual void CdrSolver::computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebG)

where ``a_G`` is the numerical representation of :math:`\mathbf{G}` over the cut-cell AMR hierarchy and must be stored on cell-centered faces, and ``a_ebG`` is the flux on the embedded boundary.
The above steps are performed by interpolating ``a_G`` to face centroids in the cut cells for computing the conservative divergence, and the remaining steps are then performed successively.
The result is put in ``a_divG``.

Note that when refinement boundaries intersect with embedded boundaries, the redistribution process is far more complicated since it needs to account for mass that moves over refinement boundaries.
These additional complicated are taken care of inside ``a_divG``, but are not discussed in detail here. 

.. caution::
   
   Mass redistribution has the effect of not being monotone and thus not TVD, and the discretization order is formally :math:`\mathcal{O}(\Delta x)`. 

.. _Chap:ExplicitAdvection:

Explicit advection
__________________

Scalar advection updates follows the computation of the explicit divergence discussed in :ref:`Chap:ExplicitDivergence`.
The face-centered fluxes :math:`\mathbf{G} = \phi\mathbf{v}` are computed by instantiation classes for the convection-diffusion-reaction solvers.
The function signature for explicit advection is

.. code-block:: c++

   void computeDivF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_dt);

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF`` using the procedure outlined above.
The argument ``a_dt`` is the time step size.
It is not needed in a method-of-lines context, but it is used in e.g. :ref:`Chap:CdrCTU` for computing transverse derivatives in order to expand the stability region (i.e., use larger CFL numbers). 

For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   CdrSolver* solver;

   const Real dt = 1.0;

   EBAMRCellData& phi  = solver->getPhi();     // Cell-centered state
   EBAMRCellData& divF = solver->getScratch(); // Scratch storage in solver
   solver->computeDivF(divF, phi, 0.0);        // Computes divF
   DataOps:incr(phi, divF, -dt);               // makes phi -> phi - dt*divF

.. _Chap:ExplicitDiffusion:
   
Explicit diffusion
__________________

Explicit diffusion is performed in much the same way as implicit advection, with the exception that the general flux :math:`\mathbf{G} = D\nabla\phi` is computed by using centered differences on face centers.
The function signature for explicit diffusion is

.. code-block:: c++

   void computeDivD(EBAMRCellData& a_divF, const EBAMRCellData& a_state);

and we increment in the same way as for explicit advection:

.. code-block:: c++

   CdrSolver* solver;

   const Real dt = 1.0;

   EBAMRCellData& phi  = solver->getPhi();     // Cell-centered state
   EBAMRCellData& divD = solver->getScratch(); // Scratch storage in solver
   solver->computeDivF(divD, phi, 0.0);        // Computes divD
   DataOps:incr(phi, divD, dt);                // makes phi -> phi + dt*divD

.. _Chap:ExplicitAdvectionDiffusion:
   
Explicit advection-diffusion
____________________________

There is also functionality for aggregating explicit advection and diffusion advances.
The reason for this is that the cut-cell overhead is only applied once on the combined flux :math:`\phi\mathbf{v} - D\nabla\phi` rather than on the individual fluxes.
For non-split methods this leads to some performance improvement since the interpolation of fluxes on cut-cell faces only needs to be performed once. 
The signature for this is precisely the same as for explicit advection only:

.. code-block:: c++
		
   void computeDivJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrapDt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF``.
For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   const Real dt = 1.0;

   EBAMRCellData& phi  = solver->getPhi();     // Cell-centered state
   EBAMRCellData& divJ = solver->getScratch(); // Scratch storage in solver
   solver->computeDivJ(divJ, phi, 0.0);        // Computes divD
   DataOps:incr(phi, divJ, -dt);               // makes phi -> phi - dt*divJ


Often, time integrators have the option of using implicit or explicit diffusion.
If the time-evolution is not split (i.e. not using a Strang or Godunov splitting), the integrators will often call ``computeDivJ`` rather separately calling ``computeDivF`` and ``computeDivD``.
If you had a split-step Godunov method, the above procedure for a forward Euler method for both parts would be:

.. code-block:: c++

   CdrSolver* solver;

   const Real dt = 1.0;

   solver->computeDivF(divF, phi, 0.0); // Computes divF = div(n*phi)
   DataOps:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF

   solver->computeDivD(divD, phi);      // Computes divD = div(D*nabla(phi))
   DataOps:incr(phi, divD, dt);         // makes phi -> phi + dt*divD

However, the cut-cell redistribution dance (flux interpolation, hybrid divergence, and redistribution) would be performed twice. 

.. _Chap:ImplicitDiffusion:

Implicit diffusion
__________________

Implicit diffusion can occasionally be necesasry. 
The convection-diffusion-reaction solvers support two basic diffusion solves:
Backward Euler and the Crank-Nicholson methods.
The function signatures for these are

.. code-block:: c++

   void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const EBAMRCellData& a_source, const Real a_dt) = 0;
   void advanceCrankNicholson(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const EBAMRCellData& a_source, const Real a_dt) = 0;		
		
		
where ``a_newPhi`` is the state at the new time :math:`t + \Delta t`, ``a_oldPhi`` is the state at time :math:`t` and ``a_source`` is the source term which strictly speaking should be centered at time :math:`t + \Delta t` for the Euler update and at time :math:`t + \Delta t/2` for the Crank-Nicholson update.
This may or may not be possible for your particular problem. 

For example, performing a split step Godunov method for advection-diffusion is as simple as:

.. code-block:: c++

   // First. Compute phi = phi - dt*div(F)
   solver->computeDivF(divF, phi, 0.0); 
   DataOps:incr(phi, divF, -dt);        

   // Implicit diffusion advance over a time step dt
   DataOps::copy(phiOld, phi);            
   solver->advanceEuler(phi, phiOld, dt); 

.. _Chap:CdrMultigrid:

CdrMultigrid
------------

``CdrMultigrid`` adds second-order accurate implicit diffusion code to ``CdrSolver``, but leaves the advection code unimplemented.
The class can use either implicit or explicit diffusion using second-order cell-centered stencils.
In addition, ``CdrMultigrid`` adds two implicit time-integrators, an implicit Euler method and a Crank-Nicholson method.

The ``CdrMultigrid`` layer uses the Helmholtz discretization discussed in :ref:`Chap:Helmholtz`.
It implements the pure functions required by :ref:`Chap:CdrSolver` but introduces a new pure function

.. code-block::

   virtual void advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_dt) = 0;

The faces states defined by the above function are used when forming a finite-volume approximation to the divergence operators, see :ref:`Chap:CdrDetails`.

.. _Chap:CdrCTU:

CdrCTU
------

``CdrCTU`` is an implementation class that uses the corner transport upwind (CTU) discretization.
The CTU discretization uses information that propagates over corners of grid cells when calculating the face states.
It can combine this with use various limiters:

* No limiter (pure CTU)
* Minmod
* Superbee
* Monotonized central differences

In addition, ``CdrCTU`` can turn off the transverse terms in which case the discretization reduces to the donor cell method.
Our motivation for using the CTU discretization lies in the time step selection for the CTU and donor-cell methods, see :ref:`Chap:CTUStep`.
Typically, we want to achieve a dimensionally independent time step that is the same in 1D, 2D, and 3D, but without directional splitting. 

Face extrapolation
__________________

The finite volume discretization uses an upstream-centered Taylor expansion that extrapolates the cell-centered term to half-edges and half-steps:

.. math::

   \phi_{i+1/2,j,}^{n+1/2} = \phi_{i,j,k}^n + \frac{\Delta x}{2}\frac{\partial \phi}{\partial x} + \frac{\Delta t}{2}\frac{\partial \phi}{\partial t} + \mathcal{O}\left(\Delta t^2\right) + \mathcal{O}\left(\Delta t\Delta x\right)

Note that the truncation order is :math:`\Delta t^2 + \Delta x\Delta t` where the latter term is due to the cross-derivative :math:`\frac{\partial^2\phi}{\partial t\partial x}`. 
The resulting expression in 2D for a velocity field :math:`\mathbf{v} = (u,v)` is

.. math::

   \phi_{i\pm1/2,j}^{n+1/2,+} = \phi_{i,j}^n \pm \frac{1}{2}\min\left[1, 1 \mp \frac{\Delta t}{\Delta x}u_{i,j}^n\right]\left(\Delta^x\phi\right)_{i,j}^n - \frac{\Delta t}{2\Delta x}v_{i,j}^n\left(\Delta^y\phi\right)_{i,j}^n,

Here, :math:`\Delta^x` are the regular (normal) slopes whereas :math:`\Delta^y` are the transverse slopes.
The transverse slopes are given by

.. math::

   (\Delta^y\phi)_{i,j}^n = \begin{cases}
   \phi_{i,j+1}^n - \phi_{i,j}^n, & v_{i,j}^n < 0 \\
   \phi_{i,j}^n - \phi_{i,j-1}^n, & v_{i,j}^n > 0 \\   
   \end{cases}

.. _Chap:CTUSlopes:

Slopes
______

For the normal slopes, the user can choose between the minmod, superbee, and monotonized central difference (MC) slopes.
Let :math:`\Delta_l = \phi_{i,j}^n - \phi_{i-1,j}^n` and :math:`\Delta_r = \phi_{i+1,j}^n - \phi_{i,j}^n`.
The slopes are given by:

.. math::

   \text{minmod: }\left(\Delta^x\phi\right)_{i,j}^n &= \begin{cases}
   \Delta_l & |\Delta_l| < |\Delta_r| \text { and } \Delta_l\Delta_r > 0 \\
   \Delta_r & |\Delta_l| > |\Delta_r| \text { and } \Delta_l\Delta_r > 0 \\
   0 & \text{otherwise}.
   \end{cases} \\[1ex]
   \text{MC: } \left(\Delta^x\phi\right)_{i,j}^n &= \text{sgn}\left(\Delta_l + \Delta_r\right)\min\left(\left|\frac{\Delta_l + \Delta_r}{2}\right|, 2\left|\Delta_l\right|, 2\left|\Delta_r\right|\right),\\[1ex]
   \text{superbee: }\left(\Delta^x\phi\right)_{i,j}^n &= \begin{cases}
   \Delta_1 & |\Delta_1| > |\Delta_2| \text { and } \Delta_1\Delta_2 > 0 \\
   \Delta_2 & |\Delta_1| < |\Delta_2| \text { and } \Delta_1\Delta_2 > 0 \\
   0 & \text{otherwise},
   \end{cases} \\[1ex]   

where for the superbee slope we have :math:`\Delta_1 = \text{minmod}\left(\Delta_l, 2\Delta_r\right)` and :math:`\Delta_2 = \text{minmod}\left(\Delta_r, 2\Delta_l\right)`.

.. note::

   When using slopes, monotonicity is not guaranteed for the CTU discretization.
   If slopes are turned off, however, the scheme is guaranteed to be monotone. 

.. _Chap:CTUStep:

Time step limitation
____________________

The stability region for the donor-cell and corner transport upwind methods are:

.. math::

   \text{Donor-cell}: \Delta t \leq \frac{\Delta x}{|v_x| + |v_y| + |v_z|}

   \text{CTU}: \Delta t \leq \frac{\Delta x}{\max\left(|v_x|,|v_y|,|v_z|\right)}

Note that when the flow is diagonal to the grid, i.e. :math:`|v_x| = |v_y| = |v_z|`, the CTU can use a time step that is three times larger than for the donor-cell method.

Class options
_____________

When running the ``CdrCTU`` solver the user can adjust the advective algorithm by turning on/off slope limiters and the transverse term through the class options:

* ``CdrCTU.slope_limiter``, which must be *none*, *minmod*, *mc*, or *superbee*.
* ``CdrCTU.use_ctu``, which must be *true* or *false*.

If the transverse terms are turned off, ``CdrCTU`` will compute the donor-cell time step.

.. _Chap:CdrGodunov:

CdrGodunov
----------

``CdrGodunov`` inherits from ``CdrMultigrid`` and adds advection code for Godunov methods.
This class borrows from ``Chombo`` internals (specifically, ``EBLevelAdvectIntegrator``) and can do second-order advection with time-extrapolation.

``CdrGodunov`` supplies (almost) the same discretization as ``CdrCTU`` with the exception that the underlying discretization can also be used for the incompressible Navier-Stokes equation.
However, it only supports the monotonized central difference limiter. 

.. caution::
   
   ``CdrGodunov`` will be removed from future versions of ``chombo-discharge``.

Using CdrSolver
----------------

Setting up the solver
_____________________

To set up the ``CdrSolver``, the following commands are usually sufficient: 

.. code-block:: c++

   // Assume m_solver and m_species are pointers to a CdrSolver and CdrSpecies
   auto solver  = RefCountedPtr<CdrSolver>  (new MyCdrSolver());
   auto species = RefCountedPtr<CdrSpecies> (new MyCdrSpecies());

To see an example, the advection-diffusion code in :file:`/Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp` shows how to set up a ``CdrSolver``. 

Filling the solver
__________________

In order to obtain mesh data from the ``CdrSolver``, the user should use the following public member functions:

.. code-block:: c++

   EBAMRCellData& getPhi();                               // Return  phi
   EBAMRCellData& getSource();                            // Returns S   
   EBAMRCellData& getCellCenteredVelocity();              // Get cell-centered velocity
   EBAMRFluxData& getFaceCenteredDiffusionCoefficient();  // Returns D
   EBAMRIVData&   getEbFlux();                            // Returns flux at EB
   EBAMRIFData&   getDomainFlux();                        // Returns flux at domain boundaries

To set the drift velocities, the user will fill the *cell-centered* velocities.
Interpolation to face-centered transport fluxes are done by ``CdrSolver`` during the discretization step.

The general way of setting the velocity is to get a direct handle to the velocity data:

.. code-block:: c++

   CdrSolver solver;
   
   EBAMRCellData& veloCell = solver.getCellCenteredVelocity();

Then, ``veloCell`` can be filled with the cell-centered velocity.
This would typically look something like this:

.. code-block:: c++

   EBAMRCellData& veloCell = m_solver->getCellCenteredVelocity();
   for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
         EBCellFAB& patchVel = (*veloCell[lvl])[dit()];

	 // Do something with patchVel
      }
   }

The same procedure goes for the source terms, diffusion coefficients, boundary conditions and so on.
For example, an explicit Euler discretization for the problem :math:`\partial_t\phi = S` is:

.. code-block:: c++

   CdrSolver* solver;

   const Real dt = 1.0;
   
   EBAMRCellData& phi = solver->getPhi();
   EBAMRCellData& src = solver->getSource();

   DataOps::incr(phi, src, dt);
   

Adjusting output
________________

It is possible to adjust solver output when plotting data.
This is done through the input file for the class that you're using.
For example, for the ``CdrCTU`` implementation the following variables are available: 

.. code-block:: text

   CdrCTU.plt_vars = phi vel src dco ebflux  # Plot variables. Options are 'phi', 'vel', 'dco', 'src'		

Here, you adjust the plotted variables by adding or omitting them from your input script.
E.g. if you only want to plot the cell-centered states you would do:

.. code-block:: bash

   CdrGodunov.plt_vars = phi  # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

.. _Chap:CdrSpecies:

CdrSpecies
----------

The ``CdrSpecies`` class is a supporting class that passes information and initial conditions into ``CdrSolver`` instances.
``CdrSpecies`` specifies whether or not the advect-diffusion solver will use only advection, diffusion, both advection and diffusion, or neither.
It also specifies initial data, and provides a string identifier to the class (e.g. for identifying output in plot files).
However, it does not contain any discretization.

.. note::

   Click `here <https://chombo-discharge.github.io/chombo-discharge/doxygen/html/classCdrSpecies.html>`_ for the ``CdrSpecies`` C++ API. 

The below code block shows an example of how to instantiate a species.
Here, diffusion code is turned off and the initial data is one everywhere. 

.. code-block:: c++

   class MySpecies : public CdrSpecies{
   public:

      MySpecies(){
         m_mobile    = true;
	 m_diffusive = false;
	 m_name      = "mySpecies";
      }

      ~MySpecies() = default;

      Real initialData(const RealVect a_pos, const Real a_time) const override {
         return 1.0;
      }
   }

.. tip::
   
   It is also possible to use computational particles as an initial condition in ``CdrSpecies``. 
   In this case you need to fill ``m_initialParticles``, and these are then deposited with a nearest-grid-point scheme when instantiating the solver.
   See :ref:`Chap:Particles` for further details.   

   
..
   Adding a stochastic flux
   ________________________

   It is possible to add a stochastic flux through the public member functions of ``CdrSolver`` in the odd case that one wants to use fluctuating hydrodynamics (FHD).
   This is done by calling a function that computes the term :math:`\sqrt{2D\phi}\mathbf{Z}`:

   .. code-block:: c++

      void gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi);

   When FHD is used, there is no guarantee that the evolution leads to non-negative values.
   We do our best to ensure that the stochastic flux is turned off when :math:`\phi \Delta V` approaches 0 by computing the face-centered states for the stochastic term using an arithmetic mean that goes to zero as :math:`\phi` approaches 0.

   In the above function, ``a_ransource`` can be used directly in a MOL context, e.g.

   .. code-block:: c++

      solver->computeDivF(divF, phi, 0.0); // Computes divF = div(n*phi)
      DataOps:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF

      solver->gwnDiffusionSource(ransource, phi);  // Compute stochastic flux
      DataOps::copy(phiOld, phi);                  // phiOld = phi - dt*divF
      DataOps::incr(phiOld, ransource, a_dt);      // phiOld = phi - dt*divF + dt*sqrt(2D*phi)Z
      solver->advanceEuler(phi, phiOld, dt);       // Backward Euler diffusion solve.

Example application(s)
----------------------

Example applications that use the ``CdrSolver`` are:

* :ref:`Chap:AdvectionDiffusionModel`.
* :ref:`Chap:CdrPlasmaModel`.    
