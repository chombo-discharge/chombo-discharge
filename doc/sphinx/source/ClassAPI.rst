.. _Chap:ClassAPI:

`PlasmaC` full API
==================

Here, we discuss the base classes that makes up the foundation of `PlasmaC`.

`PlasmaC` uses a division-of-labor between physical and numerical modules. For a brief introduction to these, see :ref:`Chap:Design`. The full, internal workings of `PlasmaC` are too complex to provide in detail, but we attempt to provide an overview here, as well as providing a summary of the input variables for the most important base classes.

If you want to view the code in full, please see the :doxy:`Doxygen API <index>`.

Here are the base modules for `PlasmaC`, note that :ref:`Chap:plasma_kinetics`, :ref:`Chap:time_stepper`, :ref:`Chap:computational_geometry`, and :ref:`Chap:cell_tagger` are abstract, and require top-level implementation.

.. _Chap:physical_domain:

physical_domain
---------------

:ref:`Chap:physical_domain` is the simplest class in `PlasmaC`. Internally, it contains only two corners of a rectangular box, and these two corners indicate the low and high ends of your physical domain. The following two options control the origin and extents of your simulation region:

.. literalinclude:: links/physical_domain.options

The :ref:`Chap:physical_domain` class was written very early on in `PlasmaC` and was intended to contain much more functionality than two corners (e.g. boundary conditions). Because the class contains so little information, the class will be integrated with :ref:`Chap:amr_mesh` in the future. 

.. _Chap:plasma_engine:

plasma_engine
-------------

:ref:`Chap:plasma_engine` is the driver class in `PlasmaC`. In order to be instantiated, :ref:`Chap:plasma_kinetics` is fed all the other modules that we discuss below. Broadly speaking, :ref:`Chap:plasma_engine` has the following responsibilities:

* If given a geometry, it will call for generation of the necessary geometric information for the EB grid. This is done by fetching the geometry from :ref:`Chap:computational_geometry`, the physical simulation domain from :ref:`Chap:physical_domain`, and the grid information from :ref:`Chap:amr_mesh`. 
* It is reponsible for setting up simulations for fresh starts, or for restarts (by reading a checkpoint file). See :ref:`Chap:RestartingSimulations` for details. 
* It calls for instantiation of solvers at appropriate times. This is done by calls to :ref:`Chap:time_stepper`, which has the ownership of the solvers.
* It keeps track of the time step, and ensures that the time in all solvers are synchronized.
* The class is responsible for regridding, which occurs through two separate mechanisms: Tagging of irregular cells and tagging of regular cells. The former is only done once, while the latter is controlled through an input parameter. For regular cells, :ref:`Chap:plasma_engine` will perform calls to (derived classes of) :ref:`Chap:cell_tagger`.
* :ref:`Chap:plasma_engine` does all the output and checkpointing.

Here are the options for plasma_engine in the current version of the code:

.. literalinclude:: links/plasma_engine.options

We now discuss these options in turn:

* ``plasma_engine.verbosity`` controls how much info is written in the pout.* files from each MPI rank. Higher numbers means more messages.
* ``plasma_engine.recursive_regrid`` is a special flag that does a recursive type of regrid where not all levels are regridded at once, but are instead regridded at time steps that are factors of ``plasma_engine.regrid_interval``. If this flag is false, all grid levels are regridded simultanously. If this flag is true, only the *finest* AMR level is regridded every ``plasma_engine.regrid_interval`` while the others are regridded at intervals defined by their refinement ratios. 
* ``plasma_engine.plot_interval`` controls how often plot files are written to disk. A negative interval turns off plot I/O. 
* ``plasma_engine.regrid_interval`` controls how often we will call for a regrid. A negative interval turns off regrids. 
* ``plasma_engine.checkpoint_interval`` controls how often checkpoint files are written to disk. A negative interval turns off file checkpointing. 
* ``plasma_engine.initial_regrids`` sets up *initial* regrids, i.e. regrids before the first time step. 
* ``plasma_engine.start_time`` is the simulation start time.
* ``plasma_engine.stop_time`` is the simulation stop time. If the simulation reaches this time, the final time step adjusted so the simulation lands on ``plasma_engine.stop_time`` 
* ``plasma_engine.max_steps`` is the maximum number of time step the simulation is allowed to run. 
* ``plasma_engine.geometry_only`` is a special flag that only plots the geometry. This will skip everything solver-related and only output a plot file of the geometry in the :file:`/geo` subfolder. 
* ``plasma_engine.ebis_memory_load_balance`` is a special flag that load balances geometry generation with respect to memory or not. 
* ``plasma_engine.write_ebis`` is a special flag that writes the EB geometry to an HDF5 file. 
* ``plasma_engine.read_ebis`` is a special flag that reads the EB geometry from an HDF5 file. 
* ``plasma_engine.output_directory`` sets the output directory for the simulation. 
* ``plasma_engine.grow_tags`` species how much to grow *irregular tags* (i.e. boundary tags). 
* ``plasma_engine.output_names`` controls the plot and checkpoint file names. 
* ``plasma_engine.max_plot_depth`` restricts the maximum plot depth. A negative number implies that plot data is taken all the way down to the finest level. 
* ``plasma_engine.max_chk_depth`` restricts the maximum checkpoint depth. A negative number implies that plot data is taken all the way down to the finest level. 
* ``plasma_engine.num_plot_ghost`` controls the number of ghost cells that will be included in plot files. 
* ``plasma_engine.plt_vars`` controls the plot variables for ``plasma_engine``. Currently available plot variables are: `tags`, which plot the cells marked for refinement. `tracer` plots the associated tracer fields from which refinement and coarsening decisions are made. `J` plots the space charge current, and `mpi_rank` plots the MPI rank for the various boxes. 
* ``plasma_engine.restart`` controls the restart step. A value less or equal to zero implies a fresh simulation. Setting this number to larger than zero tells ``plasma_engine`` to look for a checkpoint file with the specified time step, and then start the simulation from that particular time step. 
* ``plasma_engine.restart_mode`` This is a special flag that allows different restart modes (for example, removing all space charge or surface charge). 
* ``plasma_engine.refine_geometry`` controls how far the geometry will be refined. A negative values implies that level-set surfaces are refined down to the maximum allowed level. 
* ``plasma_engine.use_new_io`` is a flag that tells ``plasma_engine`` to use the new I/O routines. As these become standard, this class option will be removed. 

..
   Most options here are self-explanatory. However, we explicitly mention a few that may not be immediately clear. Firstly, ``geometry_only`` is a special option that *only* generates the geometry. It will be written to an HDF5 file whose name depends on ``output_directory`` and ``output_names``. When you develop your applications, it is often convenient to set this one to *true*, since this will skip a bunch of initialization stages (such as solving the Poisson equation, for example). It allows for quick debugging of your geometry. See chapter :ref:`Chap:ControllingOutput` to see where simulation files are placed. Furthermore, since generating the geometric information is a non-negligible work load, you may choose to write the geometric information to an HDF5 file by enable ``write_ebis``. This will write the finest level information to a file, which you can later read using ``read_ebis``. 

   The ``restart`` flag allows you to restart a simulation from a certain checkpoint step. See :ref:`Chap:RestartingSimulations` for details.

   The geometry refinement options specify to which level the geometry will be refined. Internally, this is done by tagging irregular cells down to certain levels. ``refine_geometry`` is a master option, refining all surfaces down to the specified level. The remaining options allow you to individually tag certain surfaces, such as dielectric-gas surfaces, electrode-gas surfaces, electrode-dielectric surfaces and so on. We remark that :ref:`Chap:geo_coarsener` offers a way of leveraging this coarse-grained refinement by removing tags in certain places.

.. _Chap:amr_mesh:

amr_mesh
--------

:ref:`Chap:amr_mesh` handles (almost) all spatial operations in `PlasmaC`. Internally, :ref:`Chap:amr_mesh` contains a bunch of operators that are useful across classes, such as ghost cell interpolation operators, coarsening operators, and stencils for interpolation and extrapolation near the embedded boundaries. :ref:`Chap:amr_mesh` also contains routines for generation and load-balancing of grids based and also contains simple routines for allocation and deallocation of memory. 

:ref:`Chap:amr_mesh` is an integral part of `PlasmaC`, and users will never have the need to modify it unless they are implementing something entirely new. The behavior of :ref:`Chap:amr_mesh` is modified through it's available input parameters, listed below:

.. literalinclude:: links/amr_mesh.options

We now discuss the various ``amr_mesh`` class options.

* ``amr_mesh.verbosity`` controls the verbosity of this class. ``amr_mesh`` can potentially do a lot of output, so it is best to leave this to the default value (-1) unless you are debugging. 
* ``amr_mesh.coarsest_domain`` is the partitioning of the *coarsest* grid level that discretizes your problem domain. The entries in this option must all be integers of ``amr_mesh.max_box_size``.
* ``amr_mesh.blocking_factor`` sets the minimum box size that can be generated by the mesh generation algorithm. We remark that if you are doing particle deposition, ``amr_mesh.blocking_factor`` and ``amr_mesh.max_box_size`` MUST be equal. 
* ``amr_mesh.max_box_size`` sets the maximum box size that can be generated by the mesh generation algorithm. 
* ``amr_mesh.max_ebis_box`` sets the maximum box size that will be used in the geometry generation step. A smaller box will consume less memory, but geometry generation runtime will be longer. 
* ``amr_mesh.max_amr_depth`` defines the largest possible number of grids that can be used. 
* ``amr_mesh.max_sim_depth`` defines the maximum simulation depth for the simulation. This options exists because you may want to run one part of a simulation using a coarser resolution than ``amr_mesh.max_amr_depth``. 
* ``amr_mesh.refine_all_lvl`` is deprecated class option. 
* ``amr_mesh.mg_coarsen`` is a "pre-coarsening" method for multigrid solvers. The deeper multigrids levels are there to facilitate convergence, and it often helps to use fairly large box sizes on some of these levels before aggregating boxes. 
* ``amr_mesh.fill_ratio`` is the fill ratio for the mesh refinement algorithm. This value must be between 0 and 1; a smaller value will result in larger grids. A higher value results in more compact grids, but possibly with more boxes. 
* ``amr_mesh.irreg_growth`` controls how much irregular tags (e.g. boundary tags) are grown before being passed into the mesh refinement algorithm. 
* ``amr_mesh.buffer_size`` is the minimum number of cells between grid levels. 
* ``amr_mesh.ref_rat`` is the refinement factor between levels. Values 2 and 4 are supported, and you may use mixed refinement ratios. The length of this vector must be at least equal to the number of refinement levels. 
* ``amr_mesh.num_ghost`` indicates how many ghost cells to use for all data holders. The typical value is 3 for EB-applications, but non-EB applications might get away with fewer ghost cells. 
* ``amr_mesh.eb_ghost`` controls how ghost cells are used for the EB generation. 
* ``amr_mesh.centroid_sten`` controls which stencil is used for interpolating data to irregular cell centroids. Currently available options are 'pwl' (piecewise linear), 'linear' (bi/trilinear), 'taylor' (higher order Taylor expansion), and 'lsq' which is a least squares fit. Only 'pwl' is guaranteed to have positive weights in the stencil. 
* ``amr_mesh.eb_sten`` controls which stencil is used for interpolation/extrapolation to embedded boundary centroids. We cannot guarantee that the stencils have only positive weights. 
* ``amr_mesh.redist_radius`` is the redistribution radius for hyperbolic redistribution. 
* ``amr_mesh.ghost_interp`` defines the ghost cell interpolation type. Algorithms that require very specific ghost cell interpolation schemes (advection, for example) use their own interpolation method that is outside user control. The available options are 'pwl' (piecewise linear) and 'quad' (quadratic). 
* ``amr_mesh.load_balance`` tells ``amr_mesh`` how to load balance the grids. 
* ``amr_mesh.ebcf`` allows ``amr_mesh`` to turn on certain optimizations when there are **not** crossing between embedded boundaries and grid refinement boundaries. If such crossings exist, and you set this flag to false, `PlasmaC` *will* compute incorrect answers. 

..
   In the above input parameters, ``max_amr_depth`` indicates the maximum AMR depth that we will simulate. ``coarsest_domain``, which is the number of cells on the coarsest AMR level, *must* be divisible by ``blocking_factor``, which is the smallest possible allowed box that will be generated in the grid generation. Likewise, ``max_box_size`` indicates the largest possible box. If you are using any of the particle methods in `PlasmaC`, ``blocking_factor`` and ``max_box_size`` must be equal. 

   ``ref_rat`` indicates the refinement factors between levels; the first entry indicates the refinement between the coarsest AMR level and the next. We only support refinement factors of 2 or 4 (you may use mixed refinement factors). This means that your coarsest spatial resolution will be given by the ``coarsest_domain``, and the finest resolution given by ``coarsest_domain`` multiplied recursively by your refinement ratios. Note that the resolution in each spatial direction *must* be the same.

   ``fill_ratio``, which must be :math:`>0` and :math:`\leq 1` indicates the grid generation effectiveness. For higher values of the fill ratio, grids are smaller and more compact, but more boxes are generated. For lower values, grids tend to be larger, and boxes tend to be more square. Note that if ``max_box_size`` is equal to ``blocking_factor``, the generated grids are essentially octrees.

   To add some flexibility in how refinement levels are handled in different stages of evolution, we've added an option ``max_sim_depth`` which restricts the grid generation to a level equal to or lower than ``max_amr_depth``. This options exists because checkpoint/restart (see :ref:`Chap:RestartingSimulations`) do *not* allow changing the AMR levels since this implies a changed geometry as well.

   Users may also change the stencils for data interpolation across refinement boundaries, as well as stencils for extrapolation and interpolation near the embedded boundaries. Typically, we use linear interpolation for ghost cells but quadratic interpolation is supported as well, by changing ``ghost_interp``. The flags ``stencil_order``, ``stencil_type``, and ``stencil_radius`` control the stencils used for interpolation and extrapolation near the EB. 

   If users use features that imply that the refinement boundaries might cross the EB, he **must** enable the ``ebcf`` flag. This is the case, for example, if one uses :ref:`Chap:geo_coarsener` to remove parts of the EB mesh. Internally, this flag informs :ref:`Chap:amr_mesh` about how to fill ghost and reflux data. If this flag is ``false`` *and* there are crossings between the EB and the refinement boundary, we cannot guarantee stable behavior since ghost cells might not be filled properly. 

   The options ``num_ghost`` and ``eb_ghost`` should not be changed since much of our code requires three ghost cells.

.. _Chap:computational_geometry:

computational_geometry
----------------------

:ref:`Chap:computational_geometry` is the class that implements that geometry. In `PlasmaC`, we use level-set functions for description of surfaces. :ref:`Chap:computational_geometry` is not an abstract class; if you pass in an instance of :ref:`Chap:computational_geometry` (rather than a casted instance), you will get a regular geometry without embedded boundaries. A new :ref:`Chap:computational_geometry` class requires that you set the following class members:

.. code-block:: c++

   Real m_eps0;
   Vector<electrode> m_electrodes;
   Vector<dielectric> m_dielectrics;

Here, ``m_eps0`` is hte gas permittivity, ``m_electrodes`` are the electrodes for the geometry and ``m_dielectrics`` are the dielectrics for the geometry. 

.. _Chap:electrode:

electrode
_________

The :ref:`Chap:electrode` class is responsible for describing an electrode and its boundary conditions. Internally, this class is lightweight and consists only of a tuple that holds a level-set function and an associated boolean value that tells whether or not the level-set function has a live potential or not. The constructor for the electrode class is:

.. code-block:: c++
   
  electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction = 1.0);

where the first argument is the level-set function and the second argument is responsible for setting the potential. The third argument is an optional argument that allows the user to set the potential to a specified fraction of the applied potential.

.. _Chap:dielectric:

dielectric
__________

The :ref:`Chap:dielectric` class is another lightweight class that describes a dielectric and its permittivity. The class has two constructors:


.. code-block:: c++
		
  dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity);

  dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos));

The first constructor defines a dielectric with the first argumennt describing the level-set surface and the second argument describing the permittivity. The second constructor is just like the first one, except that it allows the user to define a spatially dependent permittivity. 

.. _Chap:plasma_kinetics:

plasma_kinetics
---------------

:ref:`Chap:plasma_kinetics` is the physics module of `PlasmaC`. The entire class is an interface, whose implementations runs deep into :ref:`Chap:time_stepper` which contains much of the low-level functionality that couples the schemes. See the :doxy:`Doxygen API <classplasma__kinetics>` for details. By far, :ref:`Chap:plasma_kinetics` represents the most time-consuming tasks of implementing new plasma schemes. The reason for this is that we need to maintain a certain level of abstraction in order to cover a broad spectrum of plasma phenomena. 

There are no default input parameters for :ref:`Chap:plasma_kinetics`, as users must generally implement their own kinetics. A successful implementation of :ref:`Chap:plasma_kinetics` has the following:

* Instantiated :ref:`Chap:species`. These contain metadata for the transport solvers. 
* Instantiated :ref:`Chap:photon`. These contain metadata for the radiative transport solvers. 
* Implemented the core functionality that couple all solvers. 

`PlasmaC` automatically allocates the specified number of convection-diffusion-reaction and radiative transport solvers. For information on how to interface into the CDR solvers, see :ref:`Chap:species`. Likewise, see :ref:`Chap:photon` for how to interface into the RTE solvers.

There are currently no support for existing file formats for describing reactions and so on. If you have a huge list of reactions that need to be implemented, it would probably pay off to write a code-generating interface between :ref:`Chap:plasma_kinetics` and your list of reactions. 

Implementation of the core functionality is comparatively straightforward. In the constructor, the user should fetch his input parameters (if he has any) and **must** instantiate all species and photons in the internal containers. When the class is used by `PlasmaC` later on, all arguments in the core functions follow that ordering. The core functionality is given by the following functions: 

.. code-block:: c++
		
  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const = 0;

  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const = 0;

  virtual Vector<Real> compute_cdr_source_terms(const Real              a_time,
						const RealVect&         a_pos,
						const RealVect&         a_E,
						const RealVect&         a_gradE,
						const Vector<Real>&     a_cdr_densities,
						const Vector<Real>&     a_rte_densities,
						const Vector<RealVect>& a_grad_cdr) const = 0;

  virtual Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
						    const RealVect&     a_pos,
						    const RealVect&     a_normal,
						    const RealVect&     a_E,
						    const Vector<Real>& a_cdr_densities,
						    const Vector<Real>& a_cdr_velocities,
						    const Vector<Real>& a_cdr_gradients,
						    const Vector<Real>& a_rte_fluxes,
						    const Vector<Real>& a_extrap_cdr_fluxes) const = 0;

  virtual Vector<Real> compute_cdr_dielectric_fluxes(const Real&         a_time,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const RealVect&     a_E,
						     const Vector<Real>& a_cdr_densities,
						     const Vector<Real>& a_cdr_velocities,
						     const Vector<Real>& a_cdr_gradients,
						     const Vector<Real>& a_rte_fluxes,
						     const Vector<Real>& a_extrap_cdr_fluxes) const = 0;

  virtual Vector<Real> compute_cdr_domain_fluxes(const Real&           a_time,
						 const RealVect&       a_pos,
						 const int&            a_dir,
						 const Side::LoHiSide& a_side,
						 const RealVect&       a_E,
						 const Vector<Real>&   a_cdr_densities,
						 const Vector<Real>&   a_cdr_velocities,
						 const Vector<Real>&   a_cdr_gradients,
						 const Vector<Real>&   a_rte_fluxes,
						 const Vector<Real>&   a_extrap_cdr_fluxes) const = 0;
						     
  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const = 0;

  virtual Real initial_sigma(const Real      a_time,
			     const RealVect& a_pos) const = 0;

The above code blocks do exactly what their signatures indicate. It is up to the user to implement these. The length of the Vector holding the return values from these functions are expected to be equal to the number of CDR solvers, with the exception of *compute_rte_source_terms* which has the length given by the number of RTE solvers. Note that in all of the above, the ordering of the input vectors are expected to be the same as the ordering of the species vector of :ref:`Chap:plasma_kinetics`. 

For example, if the user has defined only a single advected species, he may implement the constructor as

.. code-block:: c++

		my_kinetics::my_kinetics(){
		   m_num_species = 1;
		   m_num_photons = 1;

		   m_species.resize(m_num_species);
		   m_photons.resize(m_num_photons);

		   m_species[0] = RefCountedPtr<species> (new my_species());
		   m_photons[0] = RefCountedPtr<photon> (new my_photon());
		}

This constructor assumes that *my_species* has already been defined somewhere (for example, as a private class within *my_kinetics*).

Defining drift velocities
_________________________

Next the user may implement the velocity computation function, which sets :math:`\mathbf{v}` in the CDR equations:

.. code-block:: c++
		
		Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
		                                        const RealVect&     a_pos,
							const RealVect&     a_E,
							const Vector<Real>& a_cdr_densities) const {
		   Vector<RealVect> velo(1);
		   velo[0] = a_E;
		   return velo;
		}

This implementation is a full implementation of the velocity coupling of the CDR equations. In this case, the velocity of the advected component is equal to :math:`\mathbf{E}`. For a full plasma simulation, there will also be mobilities involved, which the user is reponsible for obtaining.

Defining diffusion coefficients
_______________________________

In order to define diffusion coefficients, the user implements *compute_cdr_diffusion_coefficients*, which returns the diffusion coefficients for the diffused species. If a species (e.g. positive ions) is not diffusive, it does not matter what diffusion coefficient you set. 

.. code-block:: c++
		Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
		          					const RealVect&     a_pos,
							        const RealVect&     a_E,
							        const Vector<Real>& a_cdr_densities) const {
		   Vector<Real> diffco(2, 0.0);
		   diffco[0] = 1.0;
		   return diffco;
		}

Defining chemistry terms
________________________

The function *compute_cdr_source_terms* is reponsible for computing :math:`S` in the CDR equations. If we want, for example, :math:`S_1 = k n_1n_2`, where :math:`k` is some rate and :math:`n_1` and :math:`n_2` are densities of some species (e.g. electrons and positive ions), 

.. code-block:: c++
		
   Vector<Real> compute_cdr_source_terms(const Real              a_time,
		                         const RealVect&         a_pos,
					 const RealVect&         a_E,
					 const RealVect&         a_gradE,
					 const Vector<Real>&     a_cdr_densities,
					 const Vector<Real>&     a_rte_densities,
					 const Vector<RealVect>& a_grad_cdr) const {
      Vector<Real> source(m_num_species, 0.0);
      source[1] = k*a_cdr_densities[0]*a_cdr_densities[1];
      return source;
   }

In the above function, the user may also implement photoionization: The argument *Vector<Real> a_rte_densities* is the isotropic photon densities, i.e. the number of photons per unit volume. 

Defining photon production terms
________________________________

Reverse coupling between the CDR equations and the RTE equations occur through the *compute_rte_source_terms* function. The return value of this function is the mean number of photons produced per steradian. Often, such functions may be complicated. If we assume, for example, that the RTE source term is :math:`\eta = n/\tau`, where :math:`\tau` is a spontaneous emission lifetime, then we can implement the coupling as

.. code-block:: c++
		
   Vector<Real> compute_rte_source_terms(const Real&         a_time,
  					 const RealVect&     a_pos,
					 const RealVect&     a_E,
					 const Vector<Real>& a_cdr_densities) const {
      Vector<Real> source(1);
      source[0] = a_cdr_densities[0]/tau;
      return source;
   }

   Generally, one wants to ensure consistency in how one handles photon production and excited state relaxation. For the above radiative transfer example, the user should also include a corresponding term in the function that computes the chemistry source terms. 


Setting transport boundary conditions
_____________________________________
Boundary conditions are support through three functions that handle transport through three types of boundaries: domain boundaries, dielectric surfaces, and electrode surfaces. The three functions have (almost) the same signature:

.. code-block:: c++
		
   Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
		                             const RealVect&     a_pos,
					     const RealVect&     a_normal,
					     const RealVect&     a_E,
					     const Vector<Real>& a_cdr_densities,
					     const Vector<Real>& a_cdr_velocities,
					     const Vector<Real>& a_cdr_gradients,
					     const Vector<Real>& a_rte_fluxes,
					     const Vector<Real>& a_extrap_cdr_fluxes) const {
      Vector<Real> fluxes(m_num_species, 0.0);
      return fluxes;
   }

This function expects you to return the transport fluxes at the boundaries - like those occuring in a finite volume context. For domain boundaries, this signature is slightly changed: The argument ``a_normal`` (which is the normal vector *into* the gas volume) is replaced by two arguments that describe the side and direction of the domain wall. The ``a_normal`` is the normal into the gas volume, which is opposite to the convention used in finite volume formulations. The arguments ``a_cdr_velocities`` are the drift velocities projected on the outward normal, and the same convention is used for ``a_cdr_gradients`` (which hold the spatial gradients) and ``a_extrap_cdr_fluxes`` which hold the extrapolated drift fluxes. If you want simple extrapolated boundary conditions you would set one of the fluxes equal to ``a_extrap_cdr_fluxes``. A small caveat: ``a_extrap_cdr_fluxes`` are the extrapolated *drift* fluxes; if you also want the diffusive flux you can use the gradient argument and recompute the diffusion coefficient.



Setting initial surface charge
______________________________

   
Finally, the final function specifies the initial surface charge in the domain. If there is no initial surface charge, then

.. code-block:: c++
		
   Real initial_sigma(const Real      a_time,
		      const RealVect& a_pos) const {
      return 0.0
   }

.. _Chap:species:

species
_______

The :ref:`Chap:species` is a lightweight class used to provide information into convection-diffusion-reaction solvers. This class is mostly used within :ref:`Chap:plasma_kinetics` in order to provide information on how to instantiate CDR solvers. :ref:`Chap:species` is abstract so that the user must implement

.. code-block:: c++

  virtual Real initial_data(const RealVect a_pos, const Real a_time) const = 0;

This function specifies the initial data of the species that is advected. For example, the following implementation sets the initial CDR density value to one:

.. code-block:: c++
		
  Real initial_data(const RealVect a_pos, const Real a_time) const {
     return 1.0;
  }

In addition to this, the user *must* provide information on the charge of the species, and whether or not it is mobile or diffusive. In addition, he should set the name of the species so that it can be identified in output files. In `PlasmaC`, this is done by setting the following four values in the constructor

.. code-block:: c++
		
  std::string m_name; // Solver name
  int m_charge;       // Charge (in units of the elementary charge)
  bool m_diffusive;   // Diffusive species or not
  bool m_mobile;      // Mobile species or not


Usually, these are set through the constructor. The ``m_charge`` unit is in units of the elementary charge. For example, the following is a full implementation of an electron species:


.. code-block:: c++

		class electron : public species {
		  electron() {
		     m_name   = "electrons";
		     m_charge = -1;
		     m_diffusive = true;
		     m_mobile = true;
		  }

		  ~electron(){}

		  Real initial_data(const RealVect a_pos, const Real a_time) const {
		     return 1.0;
		  }
		};



The members ``m_mobile`` and ``m_diffusive`` are used for optimization in `PlasmaC`: If the user specifies that a species is immobile, `PlasmaC` will skip the advection computation. Note that ``m_diffusive`` and ``m_mobile`` override the specifications in :ref:`Chap:plasma_kinetics`. If the user provides a non-zero velocity through :ref:`Chap:plasma_kinetics` function *compute_cdr_velocities*, and sets ``m_mobile`` to ``false``, the species velocity will be zero. Of course, the user will often want to provide additional input information to his species, for example by specifying a seed for the initial conditions. Please see :ref:`Chap:MiniApplications` for how to provide input parameters. 

.. _Chap:photon:

photons
_______

:ref:`Chap:photon` is the class that supplies extra information to the RTE solvers. In those solvers, the source term computation is handled by :ref:`Chap:plasma_kinetics`, so the :ref:`Chap:photon` class is very lightweight. The user must implement a single function which specifies the absorption coefficient at a point in space:

.. code-block:: c++

		virtual const Real get_absorption_coeff(const RealVect& a_pos) const = 0;

In addition, the user should provide a name for the RTE solver so that it can be identified in the output files. This is done by setting a ``m_name`` attribute in the :ref:`Chap:photon` class.

The following is a full implementation of the :ref:`Chap:photon` class:

.. code-block:: c++

		class my_photon : public photon {
		  my_photon() {
		     m_name = "my_photon";
		  }

		  ~my_photon(){}

		  const Real get_absorption_coeff(const RealVect& a_pos) const {
		     return 1.0;
		  }
		};

By default, there are no input parameters available for the :ref:`Chap:photon` class, but the user will often want to include these, for example by modifying the absorption coefficient. Note that you are allowed to use a spatially varying absorption coefficient. Please see :ref:`Chap:MiniApplications` for how to pass input parameters into your classes.


.. _Chap:time_stepper:

time_stepper
------------

The :ref:`Chap:time_stepper` class handles the integration of the plasma equations. :ref:`Chap:time_stepper` is an abstract class for which we have several implementations that use various levels of sophistication. Writing new temporal integrators is usually an extensive task, but the base class :ref:`Chap:time_stepper` contains a lot of basic functionality (such as computing source terms for all grid levels and boxes), and also contains an interface to :ref:`Chap:plasma_kinetics`. Since :ref:`Chap:time_stepper` does not actually contain an advance method for the equations of motion, we recommend that the user refers to the API of the temporal integrator that he uses for a full explanation of the various integrators. 



The input options above are, for the most part, self-explanatory. Mostly, they refer to the handling of the size of the time step, for example by passing the Courant-Friedrichs-Lewy number, or setting a minimum or maximum possible time step. However, all of these options *may* be handled differently by different integrators, since different schemes have different restrictions on stable time steps.

Finally, there is an option to allow radiative transport updates only at certain time steps by modifying (at his own peril) the ``fast_rte`` flag. Yet again, we remark that this flag may be handled differently by different solvers. 

We have various implementation of :ref:`Chap:time_stepper` that allow different temporal integration of the equations of motion. Please see :ref:`Chap:TemporalDiscretization`.  

Typically, time steppers are selected at compile time. However, the user may select time steppers at run-time by modifying his main file in the appropriate way. For each time stepper, there are various options available at run-time through an input script.


.. _Chap:cell_tagger:

cell_tagger
-----------

The :ref:`Chap:cell_tagger` class handles tagging of cells inside the gas phase of the simulation region. Currently, we do not support tagging of cells inside the solid phase, although this would be straightforward to implement. 


In `PlasmaC`, :ref:`Chap:cell_tagger` is an abstract class that the user must implement if he wishes to change how cells are tagged. The user must defined the number of desired tracer fields (through the constructor) and then implement three functions:


.. code-block:: c++

  /*!
    @brief Compute tracer fields
  */
  virtual void compute_tracers() = 0;
		
  /*!
    @brief Coarsen a cell based on a tracer field
  */
  virtual bool coarsen_cell(const RealVect&         a_pos,
			    const Real&             a_time,
			    const Real&             a_dx,
			    const int&              a_lvl,
			    const Vector<Real>&     a_tracer,
			    const Vector<RealVect>& a_grad_tracer) = 0;

  /*!
    @brief Refine a cell based on a tracer field
  */
  virtual bool refine_cell(const RealVect&         a_pos,
			   const Real&             a_time,
			   const Real&             a_dx,
			   const int&              a_lvl,
			   const Vector<Real>&     a_tracer,
			   const Vector<RealVect>& a_grad_tracer) = 0;

Internally, the regrid happens in the following way: In the regrid stage, :ref:`Chap:plasma_engine` will perform a call to :ref:`Chap:cell_tagger` to update a number of *tracer* fields: This is done in the *compute_tracers* function above. The tracer fields are scalar fields that the user defines: A tracer field may for example be :math:`T = |\mathbf{E}|`. Once those fields are updated, :ref:`Chap:cell_tagger` will look through all cells in the domain and perform refinement and coarsening based on the implementations of *coarsen_cell* and *refine_cell*.

For example, if the user has updated a tracer :math:`T = |\mathbf{E}|`, an example implemention of *refine_cell* is

.. code-block:: c++
		
		bool refine_cell(const RealVect&         a_pos,
            		         const Real&             a_time,
				 const Real&             a_dx,
				 const int&              a_lvl,
				 const Vector<Real>&     a_tracer,
				 const Vector<RealVect>& a_grad_tracer) {

		   return (a_grad_tracer[0]*a_dx)/a_tracer[0] > 0.5) ? true : false;
		}

In this implementation, the function will examine the local curvature :math:`\left|\nabla|\mathbf{E}|\right|\Delta x/|\mathbf{E}|` and refine if it is larger than :math:`0.5`. 


For full flexibility, :ref:`Chap:cell_tagger` has been granted solver access so that the user may fetch solution data directly from the solvers. For most users, this is uneccesarily complicated. We therefore have some mid-level implementation classes that simplify this process by fetching some frequently-used data to more well-defined interfaces (as opposed to *compute_tracers*). 

There are options in the :ref:`Chap:cell_tagger` base class that permits the user to restrict tagging only to certain spatial regions by modifying the following parameters:

.. literalinclude:: links/cell_tagger.options

In the above, the user may define an arbitrary number of boxes in which tagging is *allowed*. If you do not specify a box, i.e. if ``num_boxes`` is zero, tagging is allowed everywhere. If specify one or more boxes, the ``boxN_lo`` and ``boxN_hi`` parameters indicate the valid tagging regions. Note that the boxes are not level-specific, since this is controlled through *coarsen_cell* and *refine_cell*, respectively. Again, we remark that :ref:`Chap:cell_tagger` is an abstract class, which means that these options are passed in through the *derived* class. 


.. _Chap:geo_coarsener:

geo_coarsener
-------------

:ref:`Chap:geo_coarsener` is a plug-and-play class for removing "geometric" tags that were fetched by :ref:`Chap:plasma_engine`. This class is very simple: The user defines a number of boxes in space where boundary tags are removed down to a specified level. In this way, the user may remove tags in regions where the solution is well-behaved. The interface to :ref:`Chap:geo_coarsener` is

.. literalinclude:: links/geo_coarsener.options
		    
In the above, the user may modify ``num_boxes`` in order to specify a number of boxes where tags will be removed. As input, the boxes take the low and high corners, and the specified level on which tags are removed. Note that if tags are generated on level :math:`n`, the grid will have depth :math:`n+1`. Thus, if the user wishes to remove *all* boundary tags in a box, he must set ``boxN_lvl`` to zero. 


