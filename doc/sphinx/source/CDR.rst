.. _Chap:CDR:

Convection-Diffusion-Reaction
=============================

Here, we discuss the discretization of the equation 

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\left(\mathbf{v} \phi - D\nabla \phi + \sqrt{2D\phi}\mathbf{Z}\right) = S.

We assume that :math:`\phi` is discretized by cell-centered averages (note that cell centers may lie inside solid boundaries), and use finite volume methods to construct fluxes in a cut-cells and regular cells.

Here, :math:`\mathbf{v}` indicates a drift velocity, :math:`D` is the diffusion coefficient, and the term :math:`\sqrt{2D\phi}\mathbf{Z}` is a stochastic diffusion flux. :math:`S` is the source term.

.. _Chap:cdr_solver:

cdr_solver
----------

The ``cdr_solver`` class contains the interface for solving advection-diffusion-reaction problems.
The class is abstract and resides in :file:`/src/cdr_solver/cdr_solver.H(cpp)` together with specific implementations.
By design ``cdr_solver`` does not contain any specific advective and diffusive discretization, and these are supposed to be added through inheritance.
For example, ``cdr_tga`` inherits from ``cdr_solver`` and adds a second order diffusive discretization together with multigrid code for performing implicit diffusion. 
Below that, the classes ``cdr_gdnv`` and ``cdr_muscl`` inherit everything from ``cdr_tga`` and also adds in the advective discretization.
Thus, adding new advection code is done by inheriting from ``cdr_tga`` and implementing new advection schemes.
This is much more lightweight than rewriting all of ``cdr_solver`` (which is several thousand lines of code).

.. graphviz::
   :align: center
	   
   digraph {
      rankdir="LR";
      "cdr_solver" -> "cdr_tga" -> {"cdr_gdnv", "cdr_muscl"};
   }



Currently, we mostly use the implementation given in :file:`/src/cdr_solver/cdr_gdnv.H(cpp)` which contains a second order accurate discretization with slope limiters which was distributed by the Chombo team.
The alternative implementation in :file:`/src/cdr_muscl.H(cpp)` contains a MUSCL implementation with van Leer slope limiting (i.e. much the same as the Chombo code). 

Using cdr_solver
________________

The ``cdr_solver`` is intended to be used in a method-of-lines context where the user will

1. Fill the solver with relevant data (e.g. velocities, diffusion coefficients, source terms etc.).
2. Call public member functions for computing advective or diffusive derivates, or perform implicit diffusion advances.

There are on time integration algorithms built into the ``cdr_solver``, and the user will have to supply these through ``time_stepper``.
It is up to the developer to ensure that the solver is filled with appropriate data before calling the public member functions.
This would typically look something like this:

.. code-block:: c++

   EBAMRCellData& vel = m_solver->get_velo_cell();
   for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
         EBCellFAB& patchVel = (*vel[lvl])[dit()];

	 // Set velocity of some patch
	 callSomeFunction(patchVel);
      }
   }

   // Compute div(v*phi)
   compute_divF(....)

More complete code is given in the physics module for advection-diffusion problems in :file:`/physics/advection_diffusion/advection_diffusion_stepper`.
This code is also part of a regression test found in :file:`/regression/advection_diffusion`. 

Getting data
____________

In order to obtain mesh data from the ``cdr_solver``, the user should use the following public member functions:

.. code-block:: c++

   EBAMRCellData& get_state();        // Return  phi
   EBAMRCellData& get_velo_cell();    // Get cell-centered velocity
   EBAMRFluxData& get_diffco_face();  // Returns D
   EBAMRCellData& get_source();       // Returns S
   EBAMRIVData& get_ebflux();         // Returns flux at EB
   EBAMRIFData& get_domainflux();     // Returns flux at domain boundaries

Adjusting output
________________

It is possible to adjust solver output when plotting data.
This is done through the input file for the class that you're using (e.g. :file:`/src/cdr_solver/cdr_gdnv.options`):

.. code-block:: bash

   cdr_gdnv.plt_vars = phi vel src dco ebflux  # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

Here, you adjust the plotted variables by adding or omitting them from your input script.
E.g. if you only want to plot the cell-centered states you would do:

.. code-block:: bash

   cdr_gdnv.plt_vars = phi  # Plot variables. Options are 'phi', 'vel', 'dco', 'src', 'ebflux'

.. _Chap:cdr_species:

cdr_species
-----------

The ``cdr_species`` class is a supporting class that passes information and initial conditions into ``cdr_solver`` instances. 

.. _Chap:ExplicitDivergence:   

Computing explicit divergences
------------------------------

Computing explicit divergences for equations like

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\mathbf{G} = 0

is problematic because of the arbitarily small volume fractions of cut cells.
In general, we seek to update :math:`\phi^{k+1} = \phi^k - \Delta t \left[\nabla\cdot \mathbf{G}^k\right]` where :math:`\left[\nabla\cdot\mathbf{G}\right]` is a numerical approximation based on some finite volume approximation.
Recall that in finite volume methods we usually seek the update

.. math::
   \phi^{k+1} = \phi^k - \frac{\Delta t}{\kappa \Delta x^{\textrm{DIM}}}\int_V\nabla\cdot\mathbf{G}dV,
   :label: conservativeUpdate
   
where :math:`\kappa` is the volume fraction of a grid cell, :math:`\textrm{DIM}` is the spatial dimension and the volume integral is written as discretized surface integral
   
.. math::
   \int_V\nabla\cdot\mathbf{G}dV =\sum_{f\in f(V)}\left(\mathbf{G}_f\cdot \mathbf{n}_f\right)\alpha_f\Delta x^{\textrm{DIM} -1}.
   
The sum runs over all cell edges (faces in 3D) of the cell where :math:`G_f` is the flux on the edge centroid and :math:`\alpha_f` is the edge (face) aperture.

.. figure:: figures/cutCell.png
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
		
   virtual void cdr_solver::compute_divG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebG)

where ``a_G`` is the numerical representation of :math:`\mathbf{G}` over the cut-cell AMR hierarchy and must be stored on cell-centered faces, and ``a_ebG`` is the flux on the embedded boundary.
The above steps are performed by interpolating ``a_G`` to face centroids in the cut cells for computing the conservative divergence, and the remaining steps are then performed successively.
The result is put in ``a_divG``.

Note that when refinement boundaries intersect with embedded boundaries, the redistribution process is far more complicated since it needs to account for mass that moves over refinement boundaries.
These additional complicated are taken care of inside ``a_divG``, but are not discussed in detail here. 
   
.. _Chap:NonNegative:
      
Maintaining non-negative densities
----------------------------------

Although the redistribution functionality is conservative, the cut-cells represent boundaries that make the evolution non-monotone.
In particular, explicit discretization of divergences in cut-cells do not necessarily lead to non-negative densities in the cut cells themselves.
In some cases, negative values of :math:`\phi` are non-physical and the lack of non-negativeness can lead to serious numerical issues.

In order to handle this case, we support another redistribution step in the cut cells that redistributes mass from regular cells and into the cut cells in order to maintain non-negative densities.

.. code-block:: c++
		
   void redistribute_negative(EBAMRCellData& a_phi);

Note that this *will* give give an :math:`\mathcal{O}(1)` error in the solution.
The alternative of maintaining non-negative densities through mass injection introduces the same error, but also has the side-effect of being non-conservative. 

.. _Chap:ExplicitAdvection:

Explicit advection
------------------

Scalar advective updates follows the computation of the explicit divergence discussed in :ref:`Chap:ExplicitDivergence`.
The face-centered fluxes :math:`\mathbf{G} = \phi\mathbf{v}` are computed by instantiation classes for the convection-diffusion-reaction solvers.
These solvers may compute :math:`\mathbf{G}` in different ways.
There is, for example, support for low-order upwind methods as well as Godunov methods.
The function signature for explicit advection is

.. code-block:: c++
		
   void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF`` using the procedure outlined above.
For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   // Assume that data holders divF and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined velocities. 
   solver->compute_divF(divF, phi, 0.0); // Computes divF
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells

.. _Chap:ExplicitDiffusion:
   
Explicit diffusion
------------------

Explicit diffusion is performed in much the same way as implicit advection, with the exception that the general flux :math:`\mathbf{G} = D\nabla\phi` is computed by using centered differences on face centers.
The function signature for explicit diffusion is

.. code-block:: c++
		
   void compute_divD(EBAMRCellData& a_divF, const EBAMRCellData& a_state)

and we increment in the same way as for explicit advection:

.. code-block:: c++

   // Assume that data holders divD and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined diffusion coefficients
   solver->compute_divD(divD, phi); // Computes divD
   data_ops:incr(phi, divD, dt);    // makes phi -> phi + dt*divD
   solver->redistribute_negative(phi);  // Redist negative mass in cut cells

.. _Chap:ExplicitAdvectionDiffusion:
   
Explicit advection-diffusion
----------------------------

There is also functionality for aggregating explicit advection and diffusion advances.
The reason for this is that the cut-cell overhead is only applied once on the combined flux :math:`\phi\mathbf{v} - D\nabla\phi` rather than on the individual fluxes.
For non-split methods this leads to some performance improvement since the interpolation of fluxes on cut-cell faces only needs to be performed once. 
The signature for this is precisely the same as for explicit advection only:

.. code-block:: c++
		
   void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF``.
For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   // Assume that data holders divJ and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined velocities and
   // diffusion coefficients
   solver->compute_divJ(divJ, phi, 0.0); // Computes divF
   data_ops:incr(phi, divJ, -dt);        // makes phi -> phi - dt*divJ
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells

Often, time integrators have the option of using implicit or explicit diffusion.
If the time-evolution is non-split (i.e. not using a Strang or Godunov splitting), the integrators will often call ``compute_divJ`` rather separately calling ``compute_divF`` and ``compute_divD``.
If you had a split-step Godunov method, the above procedure for a forward Euler method for both parts would be:

.. code-block:: c++

   solver->compute_divF(divF, phi, 0.0); // Computes divF = div(n*phi)
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF

   solver->compute_divD(divD, phi);      // Computes divD = div(D*nabla(phi))
   data_ops:incr(phi, divD, dt);         // makes phi -> phi + dt*divD
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells

However, the cut-cell redistribution dance (flux interpolation, hybrid divergence, and redistribution) would be performed twice. 

.. _Chap:ImplicitDiffusion:

Implicit diffusion
------------------

Occasionally, the use of implicit diffusion is necessary.
The convection-diffusion-reaction solvers support two basic diffusion solves:
Backward Euler and the Twizel-Gumel-Arigu (TGA) methods (it should be straightforward for the user to change the backward Euler method into a Crank-Nicholson scheme).
The function signatures for these are

.. code-block:: c++
		
   void advance_euler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
   void advance_tga(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
		
   void advance_euler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
   void advance_tga(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
		
where ``phiNew`` is the state at the new time :math:`t + \Delta t`, ``phiOld`` is the state at time :math:`t` and ``src`` is the source term which strictly speaking should be centered at time :math:`t + \Delta t` for the Euler update and at time :math:`t + \Delta t/2` for the TGA update.
This may or may not be possible for your particular problem. 

For example, performing a split step Godunov method for advection-diffusion is as simple as:

.. code-block:: c++

   solver->compute_divF(divF, phi, 0.0); // Computes divF = div(n*phi)
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells
		
   data_ops::copy(phiOld, phi);            // Copy state
   solver->advance_euler(phi, phiOld, dt); // Backward Euler diffusion solve

   
Adding a stochastic flux
------------------------

It is possible to add a stochastic flux through the public member functions of ``cdr_solver`` in the odd case that one wants to use fluctuating hydrodynamics (FHD).
This is done by calling a function that computes the term :math:`\sqrt{2D\phi}\mathbf{Z}`:

.. code-block:: c++
		
   void GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states);

When FHD is used, there is no guarantee that the evolution leads to non-negative values.
We do our best to ensure that the stochastic flux is turned off when :math:`\phi \Delta V` approaches 0 by computing the face-centered states for the stochastic term using an arithmetic mean that goes to zero as :math:`\phi` approaches 0.

In the above function, ``a_ransource`` can be used directly in a MOL context, e.g.

.. code-block:: c++

   solver->compute_divF(divF, phi, 0.0); // Computes divF = div(n*phi)
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->redistribute_negative(phi);	 // Redist negative mass in cut cells

   solver->GWN_diffusion_source(ransource, phi); // Compute stochastic flux
   data_ops::copy(phiOld, phi);                  // phiOld = phi - dt*divF
   data_ops::incr(phiOld, ransource, a_dt);      // phiOld = phi - dt*divF + dt*sqrt(2D*phi)Z
   solver->advance_euler(phi, phiOld, dt);       // Backward Euler diffusion solve
