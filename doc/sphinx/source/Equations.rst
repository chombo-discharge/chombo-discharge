.. _Chap:Equations:

The `PlasmaC` equation set
============================

`PlasmaC` aims at being a moderately flexible framework for fluid plasma simulations. There are several abstractions in place that ensure that the code covers non-trivial geometries, multiple time stepping schemes, and sophisticated plasma-kinetic couplings. The equation set that `PlasmaC` (currently) solves is

.. math::
   :nowrap:

   \begin{align}
   &\nabla\cdot\left(\epsilon_r\nabla\cdot\Phi\right) = -\frac{\rho}{\epsilon_0}, \\[1ex]
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\[1ex]

   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n\right) = S.
   \end{align}

This must be supported by additional boundary conditions on electrodes and insulating surfaces (on which the surface charge :math:`\sigma` lives). In addition, the user can choose to include radiative transport, which is done either in the diffusive approximation or by means of Monte Carlo methods.

Diffusive RTE methods involve solving

.. math::
   :nowrap:

   \begin{align}
      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) &= \frac{\eta}{c},
   \end{align}
   
where :math:`\Psi` is the isotropic photon density, :math:`\kappa` is an absorption length and :math:`\eta` is an isotropic source term. As an alternative, we also provide discrete photon methods that solve for the photoionization profile on a mesh by sampling discrete photons. Our discrete photon methods are capable of including far more physics; they can easily be adapted to e.g. scattering media and also provide much better qualitative features (like shadows, for example). They are, on the other hand, inherently stochastic which implies that some extra care must be taken when integrating the equations of motion. 

The number of advected species and radiative transport equations is arbitrary; the user can select the coupling through our interface. The coupling that is (currently) available in `PlasmaC` is

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


`PlasmaC` works by embedding the equations above into an abstract C++ framework that the user must implement or reuse existing pieces of, and then compile into a *mini-application*. For most users, this will mostly include implementing a new geometry or a new plasma-kinetic scheme. It is our goal that the user does not need to worry about temporal or spatial discretization of these equations, but rather focus on the actual setup of the geometry and physics. 

.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

`PlasmaC` uses structured adaptive mesh refinement (SAMR provided by Chombo :cite:`ebchombo`. SAMR exists in two separate categories, patch-based and tree-based AMR. Patch-based AMR is the more general type and contain tree-based grids as a subset; they can use refinement factors other than 2, as well as accomodate anisotropic resolutions and non-cubic patches. In patch-based AMR the domain is subdivided into a collection of hierarchically nested overlapping patches (or boxes). Each patch is a rectangular block of cells which, in space, exists on a subdomain of the union of patches with a coarser resolution. Patch-based grids generally do not have unique parent-children relations: A fine-level patch may have multiple coarse-level parents. An obvious advantage of a patch-based approach is that entire Cartesian blocks are sent into solvers, and that the patches are not restricted to squares or cubes. A notable disadvantage is that the overlapping grids inflate memory, and that additional logic is required when updating a coarse grid level from the overlapping region of a finer level. Tree-based AMR use quadtree or octree data structures that describe a hierarchy of unique parent-children relations throughout the AMR levels: Each child has exactly one parent, whereas each parent has multiple children (4 in 2D, 8 in 3D). For CPU cache performance reasons, the leaves of an octree are often cubic patches (e.g. :math:`4^3` or :math:`8^3` boxes), but the mesh can also be refined on a cell-by-cell basis. However, the use of single cell leaves becomes prohibitive at large scale for two reasons. The first is that special case must be taken in order to avoid memory inflation due to a growing tree structure. The second is that such trees, while still being SAMR, use indirect memory referencing, thus adding latency in data accessing and processing. This typically leads to poorer CPU performance since the data defined in neighboring cells may be stored on different cache lines. In `PlasmaC` and Chombo, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries. 

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit. 

.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for `PlasmaC` follows that of Chombo. In Chombo, the geometries are generated from a function :math:`f(\mathbf{x}) = 0` that describes the level-set surface. This is done by first constructing a set of boxes that covers the finest AMR level. If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box. The box is allocated in full, so using a smaller box will reduce the memory consumption. Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function. There are no graphs in boxes that are all-covered or all-regular. Furthermore, geometric data describes by the graph only exists in the cut cells themselves, so that this data is sparse. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible. Arbitrarily with fine grids geometries are not possible. Consider for example a cubic domain of :math:`(16384)^3` cells which is decomposed into :math:`(64)^3` cell size patches. This yields :math:`(256)^3` patches. Now consider that this domain is cut in half along one of the coordinate basis vectors by a planar level set surface. This surface will require allocation of :math:`256\times256\times 1` patches for the geometry. If each patch is padded with 4 ghost cells, this yields :math:`256\times256\times(72)^3 \approx 24\times 10^9` cells. Inside each cell we must store volume fractions, area fractions, cell centroids positions and so one. The required memory easily ranges in the terabyte range. 

.. _Chap:AdvectiveDiscretization:

Advective discretization
------------------------

Here, we discuss the discretization of advective derivates

.. math::
   \frac{\partial \phi}{\partial t} + \nabla\cdot\left(\mathbf{v}\phi\right) = 0

We assume that :math:`\phi` is discretized by cell-centered averages (note that cell centers may lie inside solid boundaries). We use the finite volume method to construct fluxes in a cut cell and discretize the advective derivative as

.. math::
   \int_V\nabla\cdot\left(\mathbf{v}\phi\right)dV =\sum_{f\in f(V)}\left(\mathbf{v}_f\cdot \mathbf{n}_f\right)\phi_f\alpha_f\Delta x^{D -1},
   
where the sum runs over all cell edges (faces in 3D) of the cell, :math:`F_f(\phi) = \left(\mathbf{v}_f\cdot \mathbf{n}_f\right)\phi_f` is the edge (face) centroid flux, :math:`\alpha_f` is the edge (face) aperture, and :math:`D` is the dimension. The evaluation of this expression requires knowledge of the state at the face, which in the current version of `PlasmaC` is given by a Godunov method.  

.. figure:: figures/cutCell.png
   :width: 480px
   :align: center

The possibility of arbitrarily small volume fractions :math:`\kappa` requires modification of the advective discretization in the cut cells. We use the Chombo approach and expand the range of influence of the cut cells. First, we compute the conservative divergence

.. math::
  D_{\mathbf{i}}^c(\phi) =  \sum_fF_f(\phi)\alpha_f\Delta x^{D -1}.

Next, we compute a non-conservative divergence :math:`D_{\mathbf{i}}^{nc}` that uses an extended state on covered cell faces and thereby ignores the presence of the boundaries. The extended states are extrapolated from the interior. We then use a hybrid divergence

.. math::
  D_{\mathbf{i}}^H = \kappa_{\mathbf{i}} D_{\mathbf{i}}^c + (1-\kappa_{\mathbf{i}})D_{\mathbf{i}}^{nc}.

The hybrid divergence fails to conserve mass by an amount :math:`\delta M_{\mathbf{i}} = \kappa_{\mathbf{i}}\left(1-\kappa_{\mathbf{i}}\right)\left(D_{\mathbf{i}}^c - D_{\mathbf{i}}^{nc}\right)`, which is redistributed into neighboring cells that can be reached with a monotone path of radius one. Let :math:`\delta M_{\mathbf{i}, \mathbf{j}}` be the redistributed mass from :math:`\mathbf{i}` to :math:`\mathbf{j}`. The advective discretization of cell :math:`\mathbf{j}` is then

.. math::
   D_{\mathbf{j}} = D_{\mathbf{j}}^H + \delta M_{\mathbf{i}, \mathbf{j}}.

With these definitions, the forward Euler method on :math:`\partial_t\phi = \nabla\cdot\left(\mathbf{v} \phi\right)` can now be written as :math:`\phi_{\mathbf{i}}^{n+1} = \phi_{\mathbf{i}}^n + \Delta t D_{\mathbf{i}}`. 

Charge injection and extraction in `PlasmaC` is currently handled through the advective discretization. In the future, there might exist solvers options to injects this charge though the diffusion operator instead. This would be straightforward to modify in the `PlasmaC` source code. To construct boundary fluxes, the user computes :math:`F_{\textrm{EB}}` through the physics module :ref:`Chap:plasma_kinetics`. This provides a straightforward way of handling charge injection boundary conditions. 

In order to conserve charge on solid insulators, `PlasmaC` always updates the total injection current as

.. math::
   F_\sigma(\phi) = \sum_{\phi}q_\phi F_{\textrm{EB}}(\phi),

where :math:`q_\phi` is the charge of a species :math:`\phi`. This ensures strong conservation on insulating surfaces.

.. _Chap:EllipticDiscretization:

Elliptic discretization
-----------------------

The elliptic discretization in `PlasmaC` follows the Chombo cut-cell approach where cell-centered data is used to construct face centroid centered fluxes. 

Next, we discuss the discretization of the Helmholtz equation

.. math::
   \alpha a(\mathbf{x})\phi + \beta\nabla\cdot\left(b(\mathbf{x})\phi\right) = \rho.
   
For example, the Poisson equation is represented by :math:`\alpha = 0`, :math:`\beta = -\epsilon_0`, :math:`b(\mathbf{x}) = \epsilon_r(\mathbf{x})`. Furthermore temporal discretizations of parabolic equations are also underpinned by a Helmholtz solver. 

We use the finite volume method for the Helmholtz equation. For ease of notation, we restrict the discussion below to the case :math:`a=0` which yields the Poisson equation. Extensions to the full Helmholtz problem is straightforward by adding in another diagonal term. Our implementation of the Helmholtz equation also supports multi-fluids, i.e. cases in which :math:`b(\mathbf{x})` is additionally discontinuous across a level-set surface. The multifluid problem needs additional encapsulation of a quasi-boundary condition on the interface between two materials :math:`p` and :math:`p^\prime`, given by

.. math::
   b_p\frac{\partial \phi}{\partial n_p} +   b_{p^\prime}\frac{\partial \phi}{\partial n_{p^\prime}} = \sigma,

where :math:`\mathbf{n}_p` and :math:`\mathbf{n}_{p^\prime}` are unit normals that point into each fluid, with :math:`\mathbf{n}_{p^\prime} = -\mathbf{n}_p`, and :math:`\sigma` is a surface source term. In integral, the Poisson equation is

.. math::
   \oint_A b(\mathbf{x})\nabla\phi\cdot d\mathbf{A} = \frac{1}{\beta}\int_V\rho d V. 


We consider the cell shown in the figure above. Here, the volume :math:`V_{\mathbf{i}}` is a cut-cell at a domain boundary. Integration of the above integral equation over this cell yields

.. math::
   \oint_A b(\mathbf{x})\nabla\phi\cdot d\mathbf{A} = \left(\alpha_1F_1 + \alpha_2F_2 + \alpha_3F_3 + \alpha_{\textrm{D}}F_{\textrm{D}} + \alpha_{\textrm{EB}}F_{\textrm{EB}}\right)\Delta x,

where the fluxes are centroid-centered on their respective faces and :math:`\alpha_i` are face area fractions. The centroid fluxes are evaluated by constructing second order accurate face-centered fluxes, which are then interpolated to the respective centroids. For example, for the flux through the top face in the figure above we find a standard expression for second order accurate approximations of the first derivative:

.. math::
   F_3 = F_{i,j+\frac{1}{2}} = b_{i, j+\frac{1}{2}}\frac{\phi_{i, j+1} - \phi_{i,j}}{\Delta x},

For fluxes through face centroids we interpolate the face-centered fluxes. For example, the flux :math:`F_2` in the figure above is given by

.. math::
   F_2 = \left[F_{i+\frac{1}{2},j }(1-s) + sF_{i+\frac{1}{2}, j+1}\right],

where :math:`s` is the normalized distance from the face center to the face centroid, and :math:`F_{i+\frac{1}{2},j }` and :math:`F_{i+\frac{1}{2}, j+1}` are face-centered fluxes. 

Flux evaluation on coarse-fine boundaries is slightly more involved. The AMR way of handling this is to reflux the coarse side by setting the flux into the coarse cell to be the sum of fluxes from the abutting finer cells. In Chombo, this is done by precomputing a set of flux registers that hold the face centered fluxes on both sides of the coarse-fine interface. Refluxing is then a matter of subtracting the coarse flux from the divergence computation, and adding in the sum of the fine face fluxes. I.e. let :math:`\{f_{\textrm{f}}(f_{\textrm{c}})\}` be the set of fine faces that are obtained when coarsening of a coarse face :math:`f_{\textrm{c}}`. In the reflux step, the divergence operator in the coarse cell is modified as

.. math::
   \nabla\cdot\mathbf{F} \rightarrow \nabla\cdot\mathbf{F} + \frac{1}{\Delta x}\left(\sum_{f} F_{f} - F_c\right),

where :math:`F_{c}` and :math:`F_{f}` are the coarse and fine-face fluxes, and the sum runs over all the fine faces that abut the coarse face.

.. _Chap:EllipticBoundaryConditions:

Elliptic boundary conditions
----------------------------
Next, we discuss four types of boundary conditions for the Helmholtz equation: Neumann, Dirichlet, Robin, and multifluid type boundary conditions. For Neumann boundary conditions the domain and embedded boundary fluxes are specified directly. For Dirichlet boundary co
nditions the process is more involved. For Dirichlet conditions on domain faces we apply finite differences in order to evaluate the flux through the face. For example, for a constant Dirichlet boundary condition :math:`\phi = \phi_0` the face-centered flux at the bottom face is, to second order

.. math::
  F_{i,j-\frac{1}{2}} = -\frac{b_{i,j-\frac{1}{2}}}{\Delta x}\left(3\phi_{i,j+1} -\frac{1}{3}\phi_{i,j} - \frac{8}{3}\phi_0\right)

As with the flux :math:`F_2` on the interior face, fluxes on domain faces are also interpolated to face centroids. Thus, :math:`F_{\textrm{D}}` becomes

.. math::
  F_{\textrm{D}} = \left[F_{i,j-\frac{1}{2}}(1-t) + tF_{i-1,j-\frac{1}{2}}\right],

where :math:`t` is the distance from the face center to the face centroid.

.. figure:: figures/raycast.png
   :width: 480px
   :align: center

   Ray casting at the EB for obtaining the normal gradient.

The evaluation of Dirichlet boundary conditions on the EB is more complicated because the EB normal does not align with any of the coordinate directions. To evaluate the flux on the boundary we construct ray based or least squares based stencils for evaluating :math:`\partial_n\phi` (see \cite{Johansen1998} or \cite{ebchombo} for details). Regardless of which approach is used, we have

.. math::
  \frac{\partial\phi}{\partial n} = w_0\phi_0 + \sum_{{\mathbf{i}} \in \Psi}w_{{\mathbf{i}}}\phi_{{\mathbf{i}}},

where :math:`\phi_0` is the Dirichlet value on the boundary, :math:`w_0` is a boundary weight and :math:`\Psi` is a stencil that contains only interior points. The weights :math:`w_{{\mathbf{i}}}` are weights for these points. As an example, consider the flux in the figure above. The first order accurate partial derivative on the boundary is given by

.. math::
  \frac{\partial\phi}{\partial n} = \frac{\phi_0 - \overline{\phi}}{l},

where :math:`\overline{\phi}` is the interpolated value at the intersection of the ray and the line that connects :math:`\mathbf{x}_{i-1, j}` and :math:`\mathbf{x}_{i-1, j+1}`. Since :math:`\overline{\phi}` can be linearly interpolated by using these two interior points only, this is clearly in the form of Eq.~\eqref{eq:bndry_stencil}. The boundary derivative stencils are well separated from the boundary (i.e. they do not use the values of the irregular cell itself). For the Poisson equation this is a requirement in order to achieve good conditioning of the discretized system as the volume fraction approaches zero \cite{Johansen1998}. 

Higher-order approximations to the flux are built in a similar way by including more interior cells. In our experience, the best convergence results come from using second order accurate ray-based boundary stencils, which requires 3 ghost cells in the general case. If we cannot find a stencil for computing the normal derivative by ray-casting, which can occur if there aren't enough cells available, we use quadrant-based least squares for computing the normal derivative (again, see \cite{Johansen1998} or \cite{ebchombo}).

We have also implemented Robin boundary conditions of the type

.. math::
  a_1\phi + a_2\frac{\partial \phi}{\partial n} = a_3,

which is an appropriate type of boundary condition for the radiative transfer equation. The normal derivative is given by :math:`\partial_n\phi = (a_3 - a_1\phi)/a_2` so that extrapolation of :math:`\phi` to the boundary is sufficient for imposing the boundary flux. Our way of doing this is simply to extrapolate :math:`\phi` to the boundary by using either least squares or Taylor-based stencils. 

On multifluid boundaries the boundary condition is neither Dirichlet, Neumann, or Robin. Multifluid boundaries are more complex since the state at the boundary is not known, but rather depends on the solution inside both fluids. Our approach follows that of \cite{Crockett2011} where we first compute stencils for the normal derivative on each side of the boundary,

.. math::
  \frac{\partial\phi}{\partial n_q} = w_0^q\phi_B + \sum_{{\mathbf{i}} \in \Psi_q}w_{{\mathbf{i}}}^q\phi_{{\mathbf{i}}},

where :math:`q = p` or :math:`q=p^\prime` and :math:`\phi_B` is the solution on the surface centroid, and the stencil only reaches into one of the fluids. The linear nature of this equation allows one to obtain the surface state :math:`\phi_B` from the matching condition, which can then be eliminated in order to evaluate :math:`\partial\phi/\partial n_p`. 


.. _Chap:GMG:

Geometric multigrid
-------------------

To solve the discretized Helmholtz equation we use the geometric multigrid (GMG) solver template that ships with Chombo :cite:`ebchombo`. GMG involves smoothing of the solutions on progressively coarsened grids and is compatible with AMR. Smoothing on each level involves relaxation (e.g. Jacobi or Gauss-Seidel), which primarily reduces the magnitude of high freqency errors. Removal of low-frequency errors from the solution is much slower. Because of this, multigrid accelerates convergence by projecting the error onto a coarser grid where the error has, from the viewpoint of the grid, a shorter wavelength, making relaxation more efficient. Once a bottom grid level has been reached and an approximate bottom-level solution has been found, the error is prolongated onto a finer grid and relaxation is then re-applied. Geometric multigrid works best when the long wavelength modes of the fine grid operator are well represented as short wavelength modes on the coarse grid operator. For EB applications however, coarsening can result in the removal of finer geometric features so that the relaxation step cannot sufficiently dampen the error modes at which GMG is aimed at. Because of this, geometric multigrid for EB applications usually involve lower convergence rates between each multigrid cycle than it does for geometry-less domains and, moreover, typically involves dropping to the bottom solver sooner. Currently, we only support relaxation solvers as the bottom solver for multi-phase problems, whereas we use the built-in BiCGStab and GMRES solvers in Chombo :cite:`ebchombo` for single-phase elliptic problems. In the future, we would like to use algebraic multigrid from e.g. PETSc as a bottom solver in the V-cycle in order to enhance solver efficiency for very complex geometries. 


Radiative transfer
------------------

Diffusion approximation
_______________________

In the diffusion approximation, the radiative transport equation is

.. math::

      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is called the Eddington approximation. The radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`. In the stationary case the Eddington approximation yields a Helmholtz equation

.. math::

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

which is solved by using the multigrid methods discussed above. For fully transient radiative transport, we offer discretizations based on the backward Euler and TGA schemes as discussed above. 

Monte Carlo methods
___________________

All types of moment-closed radiative transfer equations contain nonphysical artifacts (which may or may not be acceptable). For example, in the diffusion approximation the radiative flux is :math:`F = -\frac{c}{3\kappa}\nabla \Psi`, implying that photons can leak around boundaries. I.e. the diffusion approximation does not correctly describe shadows. It is possible to go beyond the diffusion approximation by also solving for higher-order moments like the radiative flux. While such methods can describe shadows, they contain other nonphysical features.

Monte Carlo methods are offered as an alternative to the diffusion approximation. Currently, we have a fully developed stationary Monte Carlo method and a transient method (which tracks photons in time) is also under development. Neither method currently includes scattering, although this would be comparatively straightforward to incorporate. As with the diffusion approximation, we do not include interaction with the plasma state in the time-of-flight of the photon. That is, we do not support e.g. scattering of a photon off electron densities. The reason for this design choice is that the velocity of a photon is much greater than the velocity of an electron, and we would have to rebin discrete photons in parallel several thousand times for each fluid advance. Thus, once a photon is created, it is invisible for the remaining solvers until it is absorbed at a point in the mesh.

Stationary Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~

The stationary Monte Carlo method proceeds as follows.

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` where :math:`\mathcal{P}` is a Poisson distribution. The user may also choose to use pseudophotons rather than physical photons. Each photon is generated in the cell centroid :math:`\mathbf{x}_0` and given a random propagation direction :math:`\mathbf{n}`.

2. Draw a propagation distance :math:`r` by drawing random numbers from an exponential distribution :math:`p(r) = \kappa \exp\left(-\kappa r\right)`. The absorbed position of the photon is :math:`\mathbf{x} = \mathbf{x}_0 + r\mathbf{n}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}`.

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell). 

Transient Monte Carlo
~~~~~~~~~~~~~~~~~~~~~

The transient Monte Carlo method is almost identical to the stationary method, except that it does not deposit all generated photons on the mesh but tracks them through time. The transient method is implemented as follows:

1. For each cell in the mesh, draw a discrete number of photons :math:`\mathcal{P}\left(\eta \Delta V\Delta t\right)` as above, and append these to the already existing photons. Each photon is given a random creation time in :math:`\Delta t`. 
   
2. Each photon is advanced over the time step :math:`\Delta t` by a sequence of :math:`N` substeps (:math:`N` may be different for each photon).

   a. We compute :math:`N` such that we sample :math:`N\Delta \tau = \Delta t` with :math:`c\kappa\Delta\tau < 1`.

   b. A photon at position :math:`\mathbf{x}_0` is moved a distance :math:`\Delta \mathbf{x} = c\mathbf{n}\Delta\tau`. For each step we compute the absorption probability :math:`p = \kappa\left|\Delta\mathbf{x}\right|` where :math:`p\in[0,1]` is a uniform random number. If the photon is absorbed on this interval, draw a new uniform random number :math:`r \in [0,1]` and absorb the photon at the position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`. If the photon is not absorbed, it is moved to position :math:`\mathbf{x}_0 + r\Delta\mathbf{x}`.

3. Check if the path from :math:`\mathbf{x}_0` to :math:`\mathbf{x}` intersects an internal or domain boundary. If it does, absorb the photon on the boundary. If not, move the photon to :math:`\mathbf{x}`.

4. Rebin the absorbed photons onto the AMR grid. This involves parallel communication. 

5. Compute the resulting photoionization profile. The user may choose between several different deposition schemes (like e.g. cloud-in-cell). 


.. bibliography:: references.bib
