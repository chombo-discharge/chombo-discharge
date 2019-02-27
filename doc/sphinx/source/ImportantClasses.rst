.. _Chap:ImportantClasses:

Class API
=========

Here, we discuss the base classes that makes up the foundation of `PlasmaC`.

`PlasmaC` uses a division-of-labor between physical and numerical modules. For a brief introduction to these, see :ref:`Chap:Design`. The full, internal workings of `PlasmaC` are too complex to provide in detail, but we attempt to provide an overview here, as well as providing a summary of the input variables for the most important base classes.

If you want to view the code in full, please see the :doxy:`Doxygen API <index>`.

Here are the base modules for `PlasmaC`, note that :ref:`Chap:plasma_kinetics`, :ref:`Chap:time_stepper`, :ref:`Chap:computational_geometry`, and :ref:`Chap:cell_tagger` are abstract, and require top-level implementation.

.. _Chap:plasma_engine:

plasma_engine
-------------

:ref:`Chap:plasma_engine` is the most important class in `PlasmaC`, and is the object that choreographs a mini-app. In order to be instantiated, :ref:`Chap:plasma_kinetics` must be fed all the other modules that we discuss below. On the whole, :ref:`Chap:plasma_engine` is given the following responsibilities:

* Given a geometry, it will call for generation of the necessary geometric information for the EB grid. This is done by fetching the geometry from :ref:`Chap:computational_geometry`, the physical simulation domain from :ref:`Chap:physical_domain`, and the grid information from :ref:`Chap:amr_mesh`. 
* It is reponsible for setting up simulations for fresh starts, or for restarts (by reading a checkpoint file). See :ref:`Chap:RestartingSimulations` for details. 
* It calls for instantiation of solvers at appropriate times. This is done by calls to :ref:`Chap:time_stepper`, which has the ownership of the solvers.
* It keeps track of the time step, and ensures that the time in all solvers are synchronized.
* The class is responsible for regridding, which occurs through two separate mechanisms: Tagging of irregular cells and tagging of regular cells. The former is only done once, while the latter is controlled through an input parameter. For regular cells, :ref:`Chap:plasma_engine` will perform calls to (derived classes of) :ref:`Chap:cell_tagger`.
* :ref:`Chap:plasma_engine` does all the output and checkpointing.

Here are the options for plasma_engine in the current version of the code:

.. literalinclude:: links/plasma_engine.options

Most options here are self-explanatory. However, we explicitly mention a few that may not be immediately clear. Firstly, ``geometry_only`` is a special option that *only* generates the geometry. It will be written to an HDF5 file whose name depends on ``output_directory`` and ``output_names``. When you develop your applications, it is often convenient to set this one to *true*, since this will skip a bunch of initialization stages (such as solving the Poisson equation, for example). It allows for quick debugging of your geometry. See chapter :ref:`Chap:ControllingOutput` to see where simulation files are placed. Furthermore, since generating the geometric information is a non-negligible work load, you may choose to write the geometric information to an HDF5 file by enable ``write_ebis``. This will write the finest level information to a file, which you can later read using ``read_ebis``. 

The ``restart`` flag allows you to restart a simulation from a certain checkpoint step. See :ref:`Chap:RestartingSimulations` for details.

The geometry refinement options specify to which level the geometry will be refined. Internally, this is done by tagging irregular cells down to certain levels. ``refine_geometry`` is a master option, refining all surfaces down to the specified level. The remaining options allow you to individually tag certain surfaces, such as dielectric-gas surfaces, electrode-gas surfaces, electrode-dielectric surfaces and so on. We remark that :ref:`Chap:geo_coarsener` offers a way of leveraging this coarse-grained refinement by removing tags in certain places.

.. _Chap:plasma_kinetics:

plasma_kinetics
---------------

:ref:`Chap:plasma_kinetics` is the physics module of `PlasmaC`. The entire class is an interface, whose implementations run deep into :ref:`Chap:time_stepper` which contains much of the low-level functionality that couples the schemes. See the :doxy:`Doxygen API <classplasma__kinetics>` for details. By far, :ref:`Chap:plasma_kinetics` is the most time-consuming tasks of implementing new plasma schemes. The reason for this is that we need to maintain a certain level of abstraction in order to cover a broad spectrum of plasma phenomena. 

There are no default input parameters for :ref:`Chap:plasma_kinetics`, as users must generally implement their own kinetics. A successful implementation of :ref:`Chap:plasma_kinetics` has the following:

* Instantiated :ref:`Chap:species`
* Instantiated :ref:`Chap:photon`
* Implemented the core functionality

Based on these three points, `PlasmaC` will automatically allocate the specified numbers of convection-diffusion-reaction and radiative transport solvers. For information on how to interface into the CDR solvers, see :ref:`Chap:species`. Likewise, see :ref:`Chap:photon` for how to interface into the RTE solvers. 

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
						     
  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const = 0;

  virtual Real initial_sigma(const Real      a_time,
			     const RealVect& a_pos) const = 0;

The above code blocks do exactly what their signatures indicate. It is up to the user to implement these. The return values in these functions are expected to be equal to the number of CDR solvers, with the exception of *compute_rte_source_terms* which has the length given by the number of RTE solvers. Note that in all of the above, the ordering of the input vectors are expected to be the same as the ordering of the species vector of :ref:`Chap:plasma_kinetics`. 

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

The remaining functions follow the same principle as the one above. For example, the function *compute_cdr_source_terms* is reponsible for computing :math:`S` in the CDR equations. If we want, for example, :math:`S = -n`, where :math:`n` is the density of the advected species, then we may do

.. code-block:: c++
		
   Vector<Real> compute_cdr_source_terms(const Real              a_time,
		                         const RealVect&         a_pos,
					 const RealVect&         a_E,
					 const RealVect&         a_gradE,
					 const Vector<Real>&     a_cdr_densities,
					 const Vector<Real>&     a_rte_densities,
					 const Vector<RealVect>& a_grad_cdr) const {
      Vector<Real> source(1);
      source[1] = -a_cdr_densities[0];
      return source;
   }

In the above function, the user may also implement photoionization: The argument *Vector<Real> a_rte_densities* is the isotropic photon densities.

Reverse coupling between the CDR equations and the RTE equations occur through the *compute_rte_source_terms* function. Often, such functions are comparatively complicated. If we assume, for example, that the RTE source term is :math:`\eta = \alpha n`, where :math:`\alpha` is an ionization constant defined somewhere, then we can implement the coupling as

.. code-block:: c++
		
   Vector<Real> compute_rte_source_terms(const Real&         a_time,
  					 const RealVect&     a_pos,
					 const RealVect&     a_E,
					 const Vector<Real>& a_cdr_densities) const {
      Vector<Real> source(1);
      source[1] = alpha*a_cdr_densities[0];
      return source;
   }

   
Finally, the final function specifies the initial surface charge in the domain. If there is no initial surface charge, then

.. code-block:: c++
		
   Real initial_sigma(const Real      a_time,
		      const RealVect& a_pos) const {
      return 0.0
		}

.. _Chap:amr_mesh:

amr_mesh
--------

:ref:`Chap:amr_mesh` is the class that handles almost all spatial operations in `PlasmaC`. Internally, :ref:`Chap:amr_mesh` contains a bunch of operators that are useful across classes, such as ghost cell interpolation operators, coarsening operators, and stencils for interpolation and extrapolation near the embedded boundaries. :ref:`Chap:amr_mesh` also contains routines for generation and load-balancing of grids based and also contains simple routines for allocation and deallocation of memory. For details, see the :doxy:`Doxygen API <amr_mesh>`.

:ref:`Chap:amr_mesh` is an integral part of `PlasmaC`, and users will never have the need to modify it unless they are implementing entirely new solvers. The behavior of :ref:`Chap:amr_mesh` is modified through it's available input parameters, listed below:

.. literalinclude:: links/amr_mesh.options     

In the above input parameters, ``max_amr_depth`` indicates the maximum AMR depth that we will simulate. ``coarsest_domain``, which is the number of cells on the coarsest AMR level, *must* be divisible by ``blocking_factor``, which is the smallest possible allowed box that will be generated in the grid generation. Likewise, ``max_box_size`` indicates the largest possible box.

``ref_rat`` indicates the refinement factors between levels; the first entry indicates the refinement between the coarsest AMR level and the next. We only support refinement factors of 2 or 4 (you may use mixed refinement factors). This means that your coarsest spatial resolution will be given by the ``coarsest_domain``, and the finest resolution given by ``coarsest_domain`` multiplied recursively by your refinement ratios. Note that the resolution in each spatial direction *must* be the same.

``fill_ratio``, which must be :math:`>0` and :math:`\leq 1` indicates the grid generation effectiveness. For higher values of the fill ratio, grids are smaller and more compact, but more boxes are generated. For lower values, grids tend to be larger, and boxes tend to be more square. Note that if ``max_box_size`` is equal to ``blocking_factor``, the generated grids are essentially octrees.

To add some flexibility in how refinement levels are handled in different stages of evolution, we've added an option ``max_sim_depth`` which restricts the grid generation to a level equal to or lower than ``max_amr_depth``. This options exists because checkpoint/erestart (see :ref:`Chap:RestartingSimulations`) do *not* allow changing the AMR levels since this implies a changed geometry as well.

Users may also change the stencils for data interpolation across refinement boundaries, as well as stencils for extrapolation and interpolation near the embedded boundaries. Typically, we use linear interpolation for ghost cells but quadratic interpolation is supported as well, by changing ``ghost_interp``. The flags ``stencil_order``, ``stencil_type``, and ``stencil_radius`` control the stencils used for interpolation and extrapolation near the EB. 

If users use features that imply that the refinement boundaries might cross the EB, he must enable the ``ebcf`` flag. This is the case, for example, if one uses :ref:`Chap:geo_coarsener` to remove parts of the EB mesh. Internally, this flag informs :ref:`Chap:amr_mesh` about how to fill ghost and reflux data. If this flag is ``false`` *and* there are crossings between the EB and the refinement boundary, we cannot guarantee stable behavior since ghost cells might not be filled properly. 

The options ``num_ghost`` and ``eb_ghost`` should not be changed since much of our code requires three ghost cells.

.. _Chap:physical_domain:

physical_domain
---------------

:ref:`Chap:physical_domain` is the simplest class in `PlasmaC`. Internally, it contains only two corners of a rectangular box, and these two corners indicate the low and high ends of your physical domain. The following two options control the origin and extents of your simulation region:

.. literalinclude:: links/physical_domain.options

.. _Chap:time_stepper:

time_stepper
------------

The :ref:`Chap:time_stepper` class handles the integration of the plasma equations. :ref:`Chap:time_stepper` is an abstract for which we have several implements (such as a second order Runge Kutta method). Writing new temporal integrators is an extensive task, but the base class :ref:`Chap:time_stepper` contains a lot of short hand functionality for updating equations, and also contains an interface to :ref:`Chap:plasma_kinetics`. Since :ref:`Chap:time_stepper` does not actually contain an advance method for the equations of motion, we recommend that the user refers to the API of the temporal integrator that he uses for a full explanation of the various integrators. However, the following options for :ref:`Chap:time_stepper` are passed into *all* implementation classes:

.. literalinclude:: links/time_stepper.options

The input options above are, for the most part, self-explanatory. Mostly, they refer to the handling of the size of the time step, for example by passing the Courant-Friedrichs-Lewy number, or setting a minimum or maximum possible time step. However, all of these options *may* be handled differently by different integrators, since different schemes have different restrictions on stable time steps.

Finally, there is an option to allow radiative transport updates only at certain time steps by modifying (at his own peril) the ``fast_rte`` flag. Yet again, we remark that this flag may be handled differently by different solvers. 

We have various implementation of :ref:`Chap:time_stepper` that allow different temporal integration of the equations of motion. Please see :ref:`Chap:TemporalDiscretization`.  

Typically, time steppers are selected at compile time. However, the user may select time steppers at run-time by modifying his main file in the appropriate way. For each time stepper, there are various options available at run-time through an input script.


.. _Chap:computational_geometry:

computational_geometry
----------------------

:ref:`Chap:computational_geometry` is the class that implements that geometry. In `PlasmaC`, we use level-set functions for description of surfaces. Please refer to :ref:`Chap:MiniApplications` for descriptions on how to implement new geometries.

:ref:`Chap:computational_geometry` is not an abstract class; if you pass in an instance of :ref:`Chap:computational_geometry` (rather than a casted instance), you will get a regular geometry without embedded boundaries.

     By default, there are no input options available for :ref:`Chap:computational_geometry`, although inherited classes that actual implement a non-regular geometry will typically have many. Please see :ref:`Chap:MiniApplications` for further information.

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

In the above, the user may define an arbitrary number of boxes in which tagging is *allowed*. If you do not specify a box, i.e. if ``num_boxes`` is zero, tagging is allowed everywhere. If specify one or more boxes, the ``boxN_lo`` and ``boxN_hi`` parameters indicate the valid tagging regions. Note that the boxes are not level-specific, since this is controlled through *coarsen_cell* and *refine_cell*, respectively. 


.. _Chap:geo_coarsener:

geo_coarsener
-------------

:ref:`Chap:geo_coarsener` is a plug-and-play class for removing "geometric" tags that were fetched by :ref:`Chap:plasma_engine`. This class is very simple: The user defines a number of boxes in space where boundary tags are removed down to a specified level. In this way, the user may remove tags in regions where the solution is well-behaved. The interface to :ref:`Chap:geo_coarsener` is

.. literalinclude:: links/geo_coarsener.options
		    
In the above, the user may modify ``num_boxes`` in order to specify a number of boxes where tags will be removed. As input, the boxes take the low and high corners, and the specified level on which tags are removed. Note that if tags are generated on level :math:`n`, the grid will have depth :math:`n+1`. Thus, if the user wishes to remove *all* boundary tags in a box, he must set ``boxN_lvl`` to zero. 

.. _Chap:species:

species
-------

The :ref:`Chap:species` is a lightweight class used to provide information into convection-diffusion-reaction solvers. This class is mostly used within :ref:`Chap:plasma_kinetics` in order to provide information on how to instantiate CDR solvers. :ref:`Chap:species` is abstract so that the user must implement

.. code-block:: c++

   /*!
    @brief Initial data. 
    @param[in] a_pos Position
    @param[in] a_time Time
  */
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const = 0;

This function specifies the initial data of the species that is advected. For example, the following implementation sets the initial CDR density value to one:

.. code-block:: c++
		
  Real initial_data(const RealVect a_pos, const Real a_time) const {
     return 1.0;
  }

In addition to this, the user *must* provide information on the charge of the species, and whether or not it is mobile or diffusive. In addition, he should set the name of the species so that it can be identified in output files. In `PlasmaC`, this is done by setting the following four values in the constructor

.. code-block:: c++
		
  /*!
    @brief Species name
  */
  std::string m_name;
  
  /*!
    @brief Charge
  */
  int m_charge;

  /*!
    @brief Diffusive species or not
  */
  bool m_diffusive;

  /*!
    @brief Mobile species or not
  */
  bool m_mobile;


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
-------

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
