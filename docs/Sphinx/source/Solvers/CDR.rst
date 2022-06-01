.. _Chap:CDR:

Convection-Diffusion-Reaction
=============================

Here, we discuss the discretization of the equation 

.. math::
   
   \frac{\partial \phi}{\partial t} + \nabla\cdot\left(\mathbf{v} \phi - D\nabla \phi + \sqrt{2D\phi}\mathbf{Z}\right) = S.

We assume that :math:`\phi` is discretized by cell-centered averages (note that cell centers may lie inside solid boundaries), and use finite volume methods to construct fluxes in a cut-cells and regular cells.
Here, :math:`\mathbf{v}` indicates a drift velocity, :math:`D` is the diffusion coefficient, and the term :math:`\sqrt{2D\phi}\mathbf{Z}` is a stochastic diffusion flux. :math:`S` is the source term.


.. note::
   
   Using cell-centered versions :math:`\phi` might be problematic for some models since the state is extended outside the valid region.
   Models might have to recenter the state in order compute e.g. physically meaningful reaction terms in cut-cells.

Source code for the convection-diffusion-reaction solvers reside in :file:`$DISCHARGE_HOME/Source/ConvectionDiffusionReaction`. 

.. _Chap:CdrSolver:

Design
------

.. _Chap:CdrSpecies:

CdrSpecies
___________

The ``CdrSpecies`` class is a supporting class that passes information and initial conditions into ``CdrSolver`` instances.
``CdrSpecies`` specifies whether or not the advect-diffusion solver will use only advection, diffusion, both advection and diffusion, or neither.
It also specifies initial data, and provides a string identifier to the class (e.g. for identifying output in plot files).

The below code block shows an example of how to instantiate a species.
Here, diffusion code is turned off and the initial data is one everywhere. 

.. code-block:: c++

   class mySpecies : public CdrSpecies{
   public:

      mySpecies(){
         m_mobile    = true;
	 m_diffusive = false;
	 m_name      = "mySpecies";
      }

      ~mySpecies() = default;

      Real initial_data(const RealVect a_pos, const Real a_time) const override {
         return 1.0;
      }
   }

Note that you can also deposit computational particles as an initial condition.
In this case you need to fill ``m_initial_particles``.
By default, these are deposited with a nearest-grid-point scheme. 

CdrSolver
__________

The ``CdrSolver`` class contains the interface for solving advection-diffusion-reaction problems.
The class is abstract but there are currently two specific implementations of this class. 
By design ``CdrSolver`` does not contain any specific advective and diffusive discretization, and these are supposed to be added through inheritance.
For example, ``CdrTGA`` inherits from ``CdrSolver`` and adds a second order diffusion discretization.
It also add multigrid code for performing implicit diffusion. 
Below that, the classes ``CdrGodunov`` and ``CdrMuscl`` inherit everything from ``CdrTGA`` and also adds in the advective discretization.
Thus, adding new advection code is done by inheriting from ``CdrTGA`` and implementing new advection schemes.

Currently, we mostly use the ``CdrGodunov`` class which contains a second order accurate discretization with slope limiters, and the advection code for this is distributed by the ``Chombo`` team. 
The alternative implementation in :file:`/src/CdrMuscl.H(cpp)` contains a MUSCL implementation with van Leer slope limiting (i.e. much the same as the Chombo code), but it does not include extrapolation in time. 

CdrTGA
*******

``CdrTGA`` adds second-order accurate implicit diffusion code to ``CdrSolver``, but leaves the advection code unimplemented.
The class can use either implicit or explicit diffusion using second-order cell-centered stencils.
In addition, ``CdrTGA`` adds two implicit time-integrators, an implicit Euler method and the Twizel-Gumel-Arigu (TGA) method. 

CdrGodunov
**********

``CdrGodunov`` inherits from ``CdrTGA`` and adds advection code for Godunov methods.
This class borrows from ``Chombo`` internals (specifically, ``EBLevelAdvectIntegrator``) and can do second-order advection with time-extrapolation.
For example, when extrapolating cell-centered data to faces, the extrapolation can be done (with Van Leer limiters) in both space and time. 

CdrMuscl
*********

``CdrMuscl`` adds MUSCL advection code to ``CdrTGA``.
It uses the same slope limiters as ``CdrGodunov`` but can not extrapolate in time.

Implementations
_______________

To use a ``CdrSolver``, one must instantiate either ``CdrGodunov`` or ``CdrMuscl`` (which differ only in their treatment of advection).
For example:

.. code-block:: c++

   CdrSpecies* spec  = (CdrSpecies*) mySpecies();
   CdrSolver* solver = (CdrSolver*)  new CdrGodunov();

   solver->set_species(spec);

Instantiating ``CdrSolver`` or ``CdrTGA`` directly will cause compile-time errors.

Note that if you want to add new advection code to ``CdrSolver``, you may inherit from ``CdrTGA`` and implement new advection routines. 

Using CdrSolver
----------------

The ``CdrSolver`` is intended to be used in a method-of-lines context where the user will

1. Fill the solver with relevant data (e.g. velocities, diffusion coefficients, source terms etc.).
2. Call public member functions for explicit advection or diffusion, or for performing implicit diffusion advances.

It is up to the developer to ensure that the solver is filled with appropriate data before calling the public member functions.
This would typically look something like this:

.. code-block:: c++

   EBAMRCellData& vel = m_solver->getCellCenteredVelocity();
   for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
         EBCellFAB& patchVel = (*vel[lvl])[dit()];

	 // Set velocity of some patch
	 callSomeFunction(patchVel);
      }
   }

   // Compute div(v*phi)
   computeDivF(....)

There are no time integration algorithms built into the ``CdrSolver``, and the user will have to supply these through ``TimeStepper``.
More complete code is given in the physics module for advection-diffusion problems in :file:`$DISCHARGE_HOME/Physics/AdvectionDiffusion/`. 
This code is also part of a regression test found in :file:`$DISCHARGE_HOME/Regression/AdvectionDiffusion`.   


Setting up the solver
_____________________

To set up the ``CdrSolver``, the following commands are usually included in ``time_stepper::setup_solvers()``:

.. code-block:: c++

   // Assume m_solver and m_species are pointers to a CdrSolver and CdrSpecies
   m_solver  = RefCountedPtr<CdrSolver>  (new MyCdrSolver());
   m_species = RefCountedPtr<CdrSpecies> (new MyCdrSpecies());

   // Solver setup
   m_solver->setVerbosity(10);
   m_solver->setSpecies(m_species);
   m_solver->parseOptions();
   m_solver->setPhase(phase::gas);
   m_solver->setAmr(m_amr);
   m_solver->setComputational_geometry(m_compgeom);
   m_solver->sanityCheck();
   m_solver->allocateInternals();

To see an example, the advection-diffusion code in :file:`/physics/AdvectionDiffusion/AdvectionDiffusion_stepper` shows how to set up this particular solver. 

Filling the solver
__________________

In order to obtain mesh data from the ``CdrSolver``, the user should use the following public member functions:

.. code-block:: c++

   EBAMRCellData& getPhi();                               // Return  phi
   EBAMRCellData& getSource();                            // Returns S   
   EBAMRCellData& getCellCenteredVelocity();              // Get cell-centered velocity
   EBAMRFluxData& getFaceCenteredDiffusionCoefficient();  // Returns D
   EBAMRIVData& getEbFlux();                              // Returns flux at EB
   EBAMRIFData& getDomainFlux();                          // Returns flux at domain boundaries

To set the drift velocities, the user will fill the *cell-centered* velocities.
Interpolation to face-centered transport fluxes are done by ``CdrSolver`` during the discretization step.

The general way of setting the velocity is to get a direct handle to the velocity data:

.. code-block:: c++

   CdrSolver solver(...);
   
   EBAMRCellData& veloCell = solver.getCellCenteredVelocity();

Then, ``veloCell`` can be filled with the cell-centered velocity.
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
This is done through the input file for the class that you're using (e.g. :file:`/src/CdrSolver/CdrGodunov.options`):

.. code-block:: bash

   CdrGodunov.plt_vars = phi vel src dco ebflux  # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

Here, you adjust the plotted variables by adding or omitting them from your input script.
E.g. if you only want to plot the cell-centered states you would do:

.. code-block:: bash

   CdrGodunov.plt_vars = phi  # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

Discretization details
----------------------

.. _Chap:ExplicitDivergence:   

Computing explicit divergences
______________________________

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
   :width: 480px
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
      

.. note::
   Mass redistribution has the effect of not being monotone and thus not TVD, and the discretization order is formally :math:`\mathcal{O}(\Delta x)`. 
   If negative densities are a problem, the ``CdrSolver`` has an option to use mass-weighted redistribution in order to redistribute mass in the neighborhood of the cut cells.
   The default is false, in which case the redistribution uses volume-weighted redistribution.

.. _Chap:ExplicitAdvection:

Explicit advection
__________________

Scalar advection updates follows the computation of the explicit divergence discussed in :ref:`Chap:ExplicitDivergence`.
The face-centered fluxes :math:`\mathbf{G} = \phi\mathbf{v}` are computed by instantiation classes for the convection-diffusion-reaction solvers.
The function signature for explicit advection is

.. code-block:: c++
		
   void computeDivF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF`` using the procedure outlined above.
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
		
   void computeDivD(EBAMRCellData& a_divF, const EBAMRCellData& a_state)

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

Occasionally, the use of implicit diffusion is necessary.
The convection-diffusion-reaction solvers support two basic diffusion solves:
Backward Euler and the Twizel-Gumel-Arigu (TGA) methods.
The function signatures for these are

.. code-block:: c++
		
   void advanceEuler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
   void advanceTGA(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
		
   void advanceEuler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
   void advanceTGA(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
		
where ``phiNew`` is the state at the new time :math:`t + \Delta t`, ``phiOld`` is the state at time :math:`t` and ``src`` is the source term which strictly speaking should be centered at time :math:`t + \Delta t` for the Euler update and at time :math:`t + \Delta t/2` for the TGA update.
This may or may not be possible for your particular problem. 

For example, performing a split step Godunov method for advection-diffusion is as simple as:

.. code-block:: c++

   solver->computeDivF(divF, phi, 0.0); // Computes divF = div(n*phi)
   DataOps:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells
		
   DataOps::copy(phiOld, phi);            // Copy state
   solver->advanceEuler(phi, phiOld, dt); // Backward Euler diffusion solve

.. note::
   The backward Euler method can easily by turned into a Crank-Nicholson method by modifying the source term and time step. 
   

   
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

Example application
-------------------

An example application of usage of the ``CdrSolver`` is found in :ref:`Chap:AdvectionDiffusionModel`. 
