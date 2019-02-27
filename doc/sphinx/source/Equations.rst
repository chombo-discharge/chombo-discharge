.. _Chap:Equations:

The `PlasmaC` equation set
============================

`PlasmaC` aims at being a moderately flexible framework for fluid plasma simulations. There are several abstractions in place that ensure that the code covers non-trivial geometries, multiple time stepping schemes, and sophisticated plasma-kinetic couplings. The equation set that `PlasmaC` (currently) solves is

.. math::
   :nowrap:

   \begin{align}
   &\nabla\cdot\left(\epsilon_r\nabla\cdot\Phi\right) = -\frac{\rho}{\epsilon_0}, \\[1ex]
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\[1ex]
   &\kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},\\[1ex]
   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n\right) = S.
   \end{align}

This must be supported by additional boundary conditions on electrodes and insulating surfaces (on which the surface charge :math:`\sigma` lives). The number of advected species and radiative transport equations is arbitrary; the user can select the coupling through our interfaces: He can even turn off radiative transport altogether.

The coupling that is (currently) available in `PlasmaC` is

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


`PlasmaC` works by embedding the equations above into an abstract C++ framework that the user must implement or reuse existing pieces of, and then compile into a *mini-application*. For most users, this will mostly include implementing a new geometry, a new plasma-kinetic scheme, or new functions for deciding when to coarsen or refine a certain spatial region. It is our goal that the user does not need to worry about temporal or spatial discretization of these equations. 

.. _Chap:SpatialDiscretization:

Spatial discretization
----------------------

`PlasmaC` uses spatial discretization provided by Chombo. In `PlasmaC`, computations occur over a set of levels with different resolutions, where the resolution refinement between levels can be a factor 2 or 4. On each level, the mesh is described by a set of disjoint patches (rectangular box in space), where the patches are distributed among MPI processes.

Embedded boundary applications are supported by additionally describing the mesh with a graph near cut-cells. This allows us to combine the efficiency of patch-based AMR with complex geometries. 

.. figure:: figures/complex_patches.png
   :width: 480px
   :align: center

   Patch-based refinement (factor 4 between levels) of a complex surface. Each color shows a patch, which is a rectangular computational unit. 

.. _Chap:EBMesh:

Geometry generation
___________________

Geometry generation for `PlasmaC` follows that of Chombo. In Chombo, the geometries are generated from a function :math:`f(\mathbf{x})` that describes the level-set surface. This is done by first constructing a set of boxes that covers the finest AMR level. If the function intersects one of these boxes, the box will allocate a *graph* that describes the connectivity of the volume-of-fluid indices in the entire box. The box is allocated in full, so using a smaller box will reduce the memory consumption (but increase run time). Chombo uses sparse storage for the EB mesh information; graphs are only stored in boxes that intersect with the implicit function. There are no graphs in boxes that are all-covered or all-regular. Furthermore, geometric data describes by the graph only exists in the cut cells themselves, so that this data is truly sparse. 

Even with sparse storage of the graph information, the memory overhead associated with the EB graph is not negligible. The finest AMR level always dominates the memory consumption for the EB mesh, and arbitrarily fine grids are not possible. Essentially, one must expect that one dense box is allocated for each box that intersects the graph at the finest AMR level. In other words, the available system memory must be sufficiently so that one can refine the entire embedded boundary at the finest AMR level. For example, if your finest AMR grid contains :math:`10^4` cut-cell boxes each of size :math:`64^3` (both are realistic numbers), the graph will contain roughly :math:`2.6\times 10^9` graph nodes (one node per cell in each cut-cell box). The memory consumption per graph various, but if we estimate this at :math:`512` bytes per node on average, then the graph memory consumption is :math:`1.34\,\textrm{TB}`. 

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
