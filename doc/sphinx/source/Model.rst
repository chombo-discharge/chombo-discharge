.. _Chap:Equations:

The `PlasmaC` model
===================

Supported solvers
-----------------

`PlasmaC` aims at being a moderately flexible framework for fluid plasma simulations. There are several abstractions in place that ensure that the code covers non-trivial geometries, multiple time stepping schemes, and fairly general plasma-kinetic couplings. The equation set that `PlasmaC` is (currently) capable of solving is

.. math::
   :nowrap:

   \begin{align}
   &\nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0}, \\[1ex]
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\[1ex]
   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n + \sqrt{2D\phi}\mathbf{Z}\right) = S,
   \end{align}

where :math:`\sqrt{2D\phi}\mathbf{Z}` is a stochastic diffusion flux suitable for fluctuating hydrodynamics models (the user may turn off this flux). The above equations must be supported by additional boundary conditions on electrodes and insulating surfaces (on which the surface charge :math:`\sigma` lives).

Radiative transport is also supported, which is done either in the diffusive approximation or by means of Monte Carlo methods. Diffusive RTE methods involve solving

.. math::
   :nowrap:

   \begin{align}
      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) &= \frac{\eta}{c},
   \end{align}
   
where :math:`\Psi` is the isotropic photon density, :math:`\kappa` is an absorption length and :math:`\eta` is an isotropic source term. The time dependent term can be turned off and the equations can be solved stationary. As an alternative, we also provide discrete photon methods that solve for the photoionization profile on a mesh by sampling discrete photons. Our discrete photon methods are capable of including far more physics; they can easily be adapted to e.g. scattering media and also provide much better qualitative features (like shadows, for example). They are, on the other hand, inherently stochastic which implies that some extra care must be taken when integrating the equations of motion.

.. _Chap:PlasmaInterface:
      
Plasma interface
----------------

Solvers in `PlasmaC` are stand-alone solvers that require a geometry, a mesh, and boundary conditions. Although they can be run on their own, we have prepared an interface where the user can implement his own plasma problem. The coupling that is (currently) available in `PlasmaC` is

.. math::
   :nowrap:

   \begin{align}
      \epsilon_r &= \epsilon_r(\mathbf{x}), (\textrm{can additionally be discontinuous}), \\[1ex]
      \mathbf{v} &= \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      D &= \mathbf{v}\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      S &= S(t, \mathbf{x}, \mathbf{E}, \nabla\mathbf{E}, n, \nabla n, \Psi), \\[1ex]
      \eta &= \eta\left(t, \mathbf{x}, \mathbf{E}, n\right), \\[1ex]
      F &= F(t, \mathbf{x}, \mathbf{E}, n),
   \end{align}

where :math:`F` is the boundary flux on insulators or electrodes (which must be separately implemented).

`PlasmaC` works by embedding the equations above into an abstract C++ framework that the user must implement or reuse existing pieces of, and then compile into a *mini-application*. For most users, this will mostly include implementing a new geometry or a new plasma-kinetic scheme. It is possible to generate entirely new physics interfaces, too. Our goal is that the user does not need to worry about temporal or spatial discretization of these equations, but rather focus on the actual setup of the geometry and physics. 

.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

`PlasmaC` uses structured adaptive mesh refinement (SAMR provided by Chombo :cite:`ebchombo`. SAMR exists in two separate categories, patch-based and tree-based AMR. Patch-based AMR is the more general type and contain tree-based grids as a subset; they can use refinement factors other than 2, as well as accomodate anisotropic resolutions and non-cubic patches. In patch-based AMR the domain is subdivided into a collection of hierarchically nested overlapping patches (or boxes). Each patch is a rectangular block of cells which, in space, exists on a subdomain of the union of patches with a coarser resolution. Patch-based grids generally do not have unique parent-children relations: A fine-level patch may have multiple coarse-level parents. An obvious advantage of a patch-based approach is that entire Cartesian blocks are sent into solvers, and that the patches are not restricted to squares or cubes that align with the coarse-grid boundary. A notable disadvantage is that additional logic is required when updating a coarse grid level from the overlapping region of a finer level. Tree-based AMR use quadtree or octree data structures that describe a hierarchy of unique parent-children relations throughout the AMR levels: Each child has exactly one parent, whereas each parent has multiple children (4 in 2D, 8 in 3D). In `PlasmaC` and Chombo, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.



.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries. However, there is significant overhead with the embedded boundary approach and, furthermore, arbitrarily complex geometries are not possible.

`PlasmaC` offers two algorithm for AMR grid generation. Both algorithms work by taking a set of flagged cells on each grid level and generating new boxes that cover the flags. The first algorithm that we support is the classical Berger-Rigoustous grid algorithm that ships with Chombo, see the figure below. The classical Berger-Rigoustous algorithm is serial-like in the sense that is collects the flagged cells onto each MPI rank and then generates the boxes. The algorithm is typically not used at large scale because of its memory consumption. As an alternative, we also support a tiled algorithm where the grid boxes on each block are generated according to a predefined tiled pattern. If a tile contains a single tag, the entire tile is flagged for refinement. The tiled algorithm produces grids that are similar to octrees, but it is more general since it also supports refinement factors other than 2, and is not restricted to cubic domains. 

.. figure:: figures/amr.png
   :width: 240px
   :align: center

   Classical cartoon of patch-based refinement. Bold lines indicate entire grid blocks. 

.. figure:: figures/tiled.png
   :width: 360px
   :align: center

   Classical cartoon of tiled patch-based refinement. Bold lines indicate entire grid blocks. 
	   
.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for `PlasmaC` follows that of Chombo. In Chombo, the geometries are generated from an implicit function :math:`f(\mathbf{x}) = 0` that describes the level-set surface. In `PlasmaC`, we use the *signed distance* function so that kinetic solvers (like Monte-Carlo photon transport) can query the distance to the closest boundary. A signed distance function is always an implicit function, but not vice versa. 

In `Chombo`, geometry generation is done by first constructing a set of boxes that covers the finest AMR level. If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box. The box is allocated in full, so using a smaller box will reduce the memory consumption. Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function. There are no graphs in boxes that are all-covered or all-regular. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible. Memory consumption generally depends on the complexity of the geometry, and arbitrarily fine grids with cut-cell geometries are not possible. Consider for example a cubic domain of :math:`(16384)^3` cells which is decomposed into :math:`(32)^3` cell size patches. This yields :math:`(512)^3` possible patches in total. Now consider that this domain is cut in half by a plane with normal vector :math:`\mathbf{n} = \hat{\mathbf{x}}`. This surface will require allocation of :math:`512\times512\times 1` patches for the geometry. If each patch is padded with 4 ghost cells, this yields :math:`512^2\times(40)^3 \approx 1.6\times 10^{10}` cells. Inside each cell we must store volume fractions, area fractions, cell centroids positions and so one. Although the surface is simple, the required memory easily ranges in the terabyte range. 

The default load-balancing for geometry generation in `Chombo` is an even division of the uniform finest-level grid among all the available. This is a reasonable approach for porous media where the cut-cells distribute evenly through the computational domain, but the approach is not scalable for geometries that consist of small objects in otherwise large domains. To achieve scalable geometry generation, our computational geometry abstractions also support the concept of *voxels* that describe a single type of material; *inside*, *outside*, or *cut-cell*. Proper use of voxels lead to much better load balancing and usually leads to orders of magnitude improvement in the time it takes to generate a geometry. How to set up geometries is discussed more closely in :ref:`Chap:NewGeometry`.

.. _Chap:CDR:

Convection-Diffusion-Reaction Equations
---------------------------------------

Here, we discuss the discretization of the equation 

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\left(\mathbf{v}\phi - D\nabla\phi\right) = S

We assume that :math:`\phi` is discretized by cell-centered averages (note that cell centers may lie inside solid boundaries), and use finite volume methods to construct fluxes in a cut-cells and regular cells.

.. _Chap:ExplicitDivergence:   

Computing explicit divergences
______________________________

Computing explicit divergences for equations like

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\mathbf{G} = 0

is problematic because of the arbitarily small volume fractions of cut cells. In general, we seek to update :math:`\phi^{k+1} = \phi^k - \Delta t \left[\nabla\cdot \mathbf{G}^k\right]` where :math:`\left[\nabla\cdot\mathbf{G}\right]` is a numerical approximation based on some finite volume approximation. Recall that in finite volume methods we usually seek the update

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

However, taking :math:`[\nabla\cdot\mathbf{G}^k]` to be this sum leads to a time step constraint proportional to :math:`\kappa`, which can be arbitrarily small. This leads to an unacceptable time step constraint for :eq:`conservativeUpdate`. We use the Chombo approach and expand the range of influence of the cut cells in order to stabilize the discretization and allow the use of a normal time step constraint. First, we compute the conservative divergence

.. math::
  \kappa_{\mathbf{i}} D_\mathbf{i}^c =  \sum_f G_f\alpha_f\Delta x^{\textrm{DIM} -1},

where :math:`G_f = \mathbf{G}_f\cdot \mathbf{n}_f`. Next, we compute a non-conservative divergence :math:`D_{\mathbf{i}}^{nc}`

.. math::
   D_\mathbf{i}^{nc} =  \frac{\sum_{\mathbf{j}\in{N}\left(\mathbf{i}\right)}\kappa_{\mathbf{j}}D_\mathbf{i}^c}{\sum_{\mathbf{j}\in{N}\left(\mathbf{i}\right)}\kappa_{\mathbf{j}}}

where :math:`N(\mathbf{i}` indicates some neighborhood of cells around cell :math:`\mathbf{i}`. Next, we compute a hybridization of the divergences, 

.. math::
  D_{\mathbf{i}}^H = \kappa_{\mathbf{i}} D_{\mathbf{i}}^c + (1-\kappa_{\mathbf{i}})D_{\mathbf{i}}^{nc},

and perform an intermediate update
  
.. math::
   \phi_{\mathbf{i}}^{k+1} = \phi_{\mathbf{i}}^k - \Delta tD_{\mathbf{i}}^H.
   
The hybrid divergence update fails to conserve mass by an amount :math:`\delta M_{\mathbf{i}} = \kappa_{\mathbf{i}}\left(1-\kappa_{\mathbf{i}}\right)\left(D_{\mathbf{i}}^c - D_{\mathbf{i}}^{nc}\right)`. In order to main overall conservation, the excess mass is redistributed into neighboring grid cells. Let :math:`\delta M_{\mathbf{i}, \mathbf{j}}` be the redistributed mass from :math:`\mathbf{j}` to :math:`\mathbf{i}` where
   
.. math::
   \delta M_{\mathbf{i}} = \sum_{\mathbf{j} \in N(\mathbf{i})}\delta M_{\mathbf{i}, \mathbf{i}}.

This mass is used as a local correction in the vicinity of the cut cells, i.e.
   
.. math::
   \phi_{\mathbf{i}}^{k+1} \rightarrow \phi_{\mathbf{i}}^{k+1} + \delta M_{\mathbf{j}\in N(\mathbf{i}), \mathbf{i}},

where :math:`\delta M_{\mathbf{j}\in N(\mathbf{i}), \mathbf{i}}` is the total mass redistributed to cell :math:`\mathbf{i}` from the other cells. After these steps, we define
   
.. math::
   \left[\nabla\cdot\mathbf{G}^k\right]_{\mathbf{i}} \equiv \frac{1}{\Delta t}\left(\phi_{\mathbf{i}}^{k+1} - \phi_{\mathbf{i}}^k\right)

Numerically, the above steps for computing a conservative divergence of a one-component flux :math:`\mathbf{G}` are implemented in the convection-diffusion-reaction solvers, which also respects boundary conditions (e.g. charge injection). The user will need to call the function

.. code-block:: c++
		
   virtual void cdr_solver::compute_divG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebG)

where ``a_G`` is the numerical representation of :math:`\mathbf{G}` over the cut-cell AMR hierarchy and must be stored on cell-centered faces, and ``a_ebG`` is the flux on the embedded boundary. The above steps are performed by interpolating ``a_G`` to face centroids in the cut cells for computing the conservative divergence, and the remaining steps are then performed successively. The result is put in ``a_divG``. 
   
.. _Chap:NonNegative:
      
Maintaining non-negative densities
__________________________________

Although the redistribution functionality is conservative, the cut-cells represent boundaries that make the evolution non-monotone. In particular, explicit discretization of divergences in cut-cells do not necessarily lead to non-negative densities in the cut cells themselves. In some cases, negative values of :math:`\phi` are non-physical and the lack of non-negativeness can lead to serious numerical issues.

In order to handle this case, we support another redistribution step in the cut cells that redistributes mass from regular cells and into the cut cells in order to maintain non-negative densities.

.. code-block:: c++
		
   void make_non_negative(EBAMRCellData& a_phi)

Again, the functionality for redistributing negative mass in a conservative way is owned by the convection-diffusion-reaction solvers. 

.. _Chap:ExplicitAdvection:

Explicit advection
__________________

Scalar advective updates follows the computation of the explicit divergence discussed in :ref:`Chap:ExplicitDivergence`. The face-centered fluxes :math:`\mathbf{G} = \phi\mathbf{v}` are computed by instantiation classes for the convection-diffusion-reaction solvers. These solvers may compute :math:`\mathbf{G}` in different ways. There is, for example, support for low-order upwind methods as well as Godunov methods. The function signature for explicit advection is

.. code-block:: c++
		
   void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF`` using the procedure outlined above. For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   // Assume that data holders divF and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined velocities. 
   solver->compute_divF(divF, phi, 0.0); // Computes divF
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->make_non_negative(phi);	 // Redist negative mass in cut cells

.. _Chap:ExplicitDiffusion:
   
Explicit diffusion
__________________

Explicit diffusion is performed in much the same way as implicit advection, with the exception that the general flux :math:`\mathbf{G} = D\nabla\phi` is computed by using centered differences on face centers. The function signature for explicit diffusion is

.. code-block:: c++
		
   void compute_divD(EBAMRCellData& a_divF, const EBAMRCellData& a_state)

and we increment in the same way as for explicit advection:

.. code-block:: c++

   // Assume that data holders divD and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined diffusion coefficients
   solver->compute_divD(divD, phi); // Computes divD
   data_ops:incr(phi, divD, dt);    // makes phi -> phi + dt*divD
   solver->make_non_negative(phi);  // Redist negative mass in cut cells

.. _Chap:ExplicitAdvectionDiffusion:
   
Explicit advection-diffusion
____________________________

There is also functionality for aggregating explicit advection and diffusion advances. The reason for this is that the cut-cell overhead is only applied once on the combined flux :math:`\phi\mathbf{v} - D\nabla\phi` rather than on the individual fluxes. For non-split methods this leads to some performance improvement. The signature for this is precisely the same as for explicit advection only:

.. code-block:: c++
		
   void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt)

where the face-centered fluxes are computed by using the velocities and boundary conditions that reside in the solver, and result is put in ``a_divF``. For example, in order to perform an advective advance over a time step :math:`\Delta t`, one would perform the following:

.. code-block:: c++

   // Assume that data holders divJ and phi are defined, and that 'solver' is
   // a valid convection-diffusion reaction solver with defined velocities and
   // diffusion coefficients
   solver->compute_divJ(divJ, phi, 0.0); // Computes divF
   data_ops:incr(phi, divJ, -dt);        // makes phi -> phi - dt*divJ
   solver->make_non_negative(phi);	 // Redist negative mass in cut cells

Often, time integrators have the option of using implicit or explicit diffusion. If the time-evolution is non-split (i.e. not using a Strang or Godunov splitting), the integrators will often call ``compute_divJ`` rather separately calling ``compute_divF`` and ``compute_divD``. If you had a split-step Godunov method, the above procedure for a forward Euler method for both parts would be:

.. code-block:: c++

   solver->compute_divF(divF, phi, 0.0); // Computes divF = div(n*phi)
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF

   solver->compute_divD(divD, phi);      // Computes divD = div(D*nabla(phi))
   data_ops:incr(phi, divD, dt);         // makes phi -> phi + dt*divD
   solver->make_non_negative(phi);	 // Redist negative mass in cut cells

However, the cut-cell redistribution dance (flux interpolation, hybrid divergence, and redistribution) would be performed twice. 

.. _Chap:ImplicitDiffusion:

Implicit diffusion
__________________

Occasionally, the use of implicit diffusion is necessary. The convection-diffusion-reaction solvers support two basic diffusion solves: Backward Euler and the Twizel-Gumel-Arigu (TGA) methods. The function signatures for these are

.. code-block:: c++
		
   void advance_euler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
   void advance_tga(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const EBAMRCellData& src, const Real dt)
		
   void advance_euler(EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
   void advance_tga(  EBAMRCellData& phiNew, const EBAMRCellData& phiOld, const Real dt)
		
where ``phiNew`` is the state at the new time :math:`t + \Delta t`, ``phiOld`` is the state at time :math:`t` and ``src`` is the source term which strictly speaking should be centered at time :math:`t + \Delta t` for the Euler update and at time :math:`t + \Delta t/2` for the TGA update. This may or may not be possible for your particular problem. 

For example, performing a split step Godunov method for advection-diffusion is as simple as:

.. code-block:: c++

   solver->compute_divF(divF, phi, 0.0); // Computes divF = div(n*phi)
   data_ops:incr(phi, divF, -dt);        // makes phi -> phi - dt*divF
   solver->make_non_negative(phi);	 // Redist negative mass in cut cells
		
   data_ops::copy(phiOld, phi);            // Copy state
   solver->advance_euler(phi, phiOld, dt); // Backward Euler diffusion solve

.. _Chap:FieldSolver:
   
Field solver
------------

The `PlasmaC` field solver has a lot of supporting functionality, but essentially relies on only one critical function: Solving for the potential. This is done by calling a class-specific function

.. code-block:: c++

   bool solve(MFAMRCellData& phi, const MFAMRCellData& rho, const EBAMRIVData& sigma);

where ``phi`` is the resulting potential that was computing with the space charge density ``rho`` and surface charge density ``sigma``.

Currently, only one field solver is implemented and this solver uses a geometric multigrid method for solving for the potential. The solver supports three phases: electrodes, gas, and dielectric. Boundary conditions for the solver must be set by the user through an input script. 

.. _Chap:RadiativeTransfer:

Radiative transfer
------------------

Radiative transfer is supported in the diffusion (i.e. Eddington or Helmholtz) approximation and with Monte Carlo sampling of discrete photons. The solvers share a common interface but since diffusion RTE is deterministic and discrete Monte Carlo photons are stochastic, not all temporal integration methods will support both. The diffusion approximation relies on solving an elliptic equation in the stationary case and a parabolic equation in the time-dependent case, while the Monte-Carlo approach currently only solves for instantaneous photon transport. However, it would be straightforward to include transient photons. 

Diffusion approximation
_______________________

In the diffusion approximation, the radiative transport equation is

.. math::

      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is called the Eddington approximation. The radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`. We do not currently support flux-limited diffusion radiative transfer. In the stationary case this yields a Helmholtz equation

.. math::

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is solved by a geometric multigrid method. The default boundary conditions are of the Robin type. For fully transient radiative transport, we offer discretizations based on the backward Euler and TGA schemes as discussed above. 

Monte Carlo methods
___________________

All types of moment-closed radiative transfer equations contain nonphysical artifacts (which may or may not be acceptable). For example, in the diffusion approximation the radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`, implying that photons can leak around boundaries. I.e. the diffusion approximation does not correctly describe shadows. It is possible to go beyond the diffusion approximation by also solving for higher-order moments like the radiative flux. While such methods can describe shadows, they contain other nonphysical features.

.. figure:: figures/rte_comp.png
   :width: 720px
   :align: center

   Qualitative comparison between predictions made with a diffusion RTE solver and a Monte Carlo RTE solver. Left: Source term: Middle: Solution computed in the diffusion approximation with homogeneous Robin boundary conditions. Right: Solution computed with a Monte Carlo method. 

Monte Carlo methods are offered as an alternative to the diffusion approximation. Currently, we have a fully developed stationary Monte Carlo method and a transient method (which tracks photons in time) is also under development. Neither method currently includes scattering, although this would be comparatively straightforward to incorporate. As with the diffusion approximation, we do not include interaction with the plasma state in the time-of-flight of the photon. That is, we do not support e.g. scattering of a photon off electron densities. The reason for this design choice is that the velocity of a photon is much greater than the velocity of an electron, and we would have to rebin discrete photons in parallel several thousand times for each fluid advance. Thus, once a photon is created, it is invisible for the remaining solvers until it is absorbed at a point in the mesh.

Stationary Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~

The stationary Monte Carlo method proceeds as follows.

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` where :math:`\mathcal{P}` is a Poisson distribution. The user may also choose to use pseudophotons rather than physical photons by modifying photon weights. Each photon is generated in the cell centroid :math:`\mathbf{x}_0` and given a random propagation direction :math:`\mathbf{n}`.

2. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`. The absorbed position of the photon is :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}` or reflect it off symmetry boundaries. 

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).
      

The Monte Carlo methods use computational particles for advancing the photons in exactly the same way a Particle-In-Cell method would use them for advancing electrons. Although a computational photon would normally live on the finest grid level that overlaps its position, this is not practical for all particle deposition kernels. For example, for cloud-in-cell deposition schemes it is useful to have the restrict the interpolation kernels to the grid level where the particle lives. In Chombo-speak, we therefore use a buffer region that extends some cells from a refinement boundary where the photons are not allowed to live. Instead, photons in that buffer region are transferred to a coarser level, and their deposition clouds are first interpolated to the fine level before deposition on the fine level happens. Selecting a deposition scheme and adjusting the buffer region is done through an input script associated with the solver. 
   
Transient Monte Carlo
~~~~~~~~~~~~~~~~~~~~~

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time. The transient method is implemented as follows:

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` as above, and append these to the already existing photons. Each photon is given a uniformly distributed random creation time within :math:`\Delta t`. 
   
2. Each photon is advanced over the time step :math:`\Delta t` by a sequence of :math:`N` substeps (:math:`N` may be different for each photon).

   a. We compute :math:`N` such that we sample :math:`N\Delta \tau = \Delta t` with :math:`c\kappa\Delta\tau < 1`.

   b. A photon at position :math:`\mathbf{x}_0` is moved a distance :math:`\Delta \mathbf{x} = c\mathbf{n}\Delta\tau`. For each step we compute the absorption probability :math:`p = \kappa\left|\Delta\mathbf{x}\right|` where :math:`p\in[0,1]` is a uniform random number. If the photon is absorbed on this interval, draw a new uniform random number :math:`r \in [0,1]` and absorb the photon at the position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`. If the photon is not absorbed, it is moved to position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}`.

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell).
      
.. _Chap:SigmaSolver:

Surface charge solver
---------------------
In order to conserve charge on solid insulators, `PlasmaC` has a solver that is defined on the gas-dielectric interface where the surface charge is updated with the incoming flux

.. math::
   F_\sigma(\phi) = \sum_{\phi}q_\phi F_{\textrm{EB}}(\phi),

where :math:`q_\phi` is the charge of a species :math:`\phi`. This ensures strong conservation on insulating surfaces.

.. bibliography:: references.bib
