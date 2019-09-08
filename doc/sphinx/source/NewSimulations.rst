.. _Chap:NewSimulations:

Setting up simulations
======================

In `PlasmaC`, most of the action occurs through the physics interface :ref:`Chap:plasma_kinetics`.

Setting up a geometry
---------------------


Defining your chemistry
-----------------------

Chemical reactions are defined through our physics interface :ref:`Chap:plasma_kinetics`. There is support for general types of reactions amongst all species through these interfaces, but there is no middleware that translates known formats (e.g. CHEMKIN) to something usable for `PlasmaC`. If the user has a chemical database consisting of hundreds of reactions, it would probably pay off to construct such middleware first.

The must implement a set of :ref:`Chap:species` that describes the various chemical species that will be tracked. These can be coupled with radiative transport through the :ref:`Chap:photon` class. There are currently two kinetic interfaces that are supported. The first is:

.. code-block:: c++

   virtual void advance_reaction_network(Vector<Real>&          a_particle_sources,
		                         Vector<Real>&          a_photon_sources,
					 const Vector<Real>     a_particle_densities,
					 const Vector<RealVect> a_particle_gradients,
					 const Vector<Real>     a_photon_densities,
					 const RealVect         a_E,
					 const RealVect         a_pos,
					 const Real             a_dx,
					 const Real             a_dt,
					 const Real             a_time,
					 const Real             a_kappa) const = 0;

This function is called for all grid cells in a `PlasmaC` simulation. Here, the first two arguments are output arguments that hold the particle and photon sources. The third and fourth argument are input arguments that hold the densities in the grid cell. We have chosen this format since source terms can then be filled using a variety of algorithms. For example, particle source terms can be filled using reaction-rate equations, tau-leaping schemes, or even stochastic simulation algorithms. We would like to remark that the input and output from these functions can be interpreted in different ways by different solvers. For example, the Monte-Carlo radiative transfer solver can take ``a_photon_source`` to be either a number per grid cell, or a volumetric source term. For example, if you use a stochastic simulation algorithm it is natural to describe ``a_photon_sources`` as the number of photons produced in the cell, and the Monte-Carlo solver needs to be informed that its source term contains a number rather than a rate. 

Electrostatic boundary and initial conditions
---------------------------------------------

Setting the electrostatic boundary and initial conditions requires three steps:

1. You must pass a function pointer to :ref:`plasma_engine` that decribes the applied voltage :math:`V(t)`. We will refer to this function as the "live voltage".

   If you use the Python setup tool, this function will automatically be defined for you, and you may manipulate it directly in your main file.
2. Define the boundary conditions on the domain edges (faces in 3D). These have the form:

   .. code-block:: bash
		
		poisson_solver.bc_x_low  = neumann               # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
		poisson_solver.bc_x_high = neumann               # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
		poisson_solver.bc_y_low  = neumann               # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
		poisson_solver.bc_y_high = neumann               # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
		poisson_solver.bc_z_low  = dirichlet_ground      # BC type. "neumann", "dirichlet_ground", "dirichlet_live"
		poisson_solver.bc_z_high = dirichlet_live        # BC type. "neumann", "dirichlet_ground", "dirichlet_live"

3. You must supply the boundary conditions on your electrodes. This is done by defining the electrode as ``live=true`` or ``live=false``, usually through the constructor. However, you *may* apply a fraction of the live voltage :math:`V(t)` to your electrodes by setting the ``m_fraction`` class member. Setting ``m_fraction = 1.0`` will set the potential on the electrode to :math:`V(t)`, setting ``m_fraction = 0.5`` sets the potential to :math:`0.5V(t)` and so on. Please see the :ref:`Chap:electrode` chapter for more details. 

   On dielectric surfaces the electric potential is always computed based on the dielectric boundary condition, and there is not way of setting this directly. 
		

Setting initial conditions
--------------------------

In order to set the initial conditions, the user must provide an implementation of the :ref:`Chap:species` class. This implementation may exist anywhere, but only species defined in :ref:`Chap:plasma_kinetics` will be tracked in the simulation. Through :ref:`Chap:species`, the user may fill CDR solvers with a prescribed volumetric density, *or* may optionally deposit the initial conditions by depositing physical particles onto the grid. For example, the user *must* provide a function

.. code-block:: c++
		
  Real initial_data(const RealVect a_pos, const Real a_time) const {
     return something;
  }

which sets the initial density field. However, :ref:`Chap:species` may deposit particles by providing these to the instantiated object. For example, the following code block is a complete implementation that uses scalar fields *and* particles as an initial condition:

.. code-block:: c++

		class electron : public species {
		  electron() {
		     m_name       = "electrons";
		     m_charge     = -1;
		     m_diffusive  = true;
		     m_mobile     = true;
		     m_deposition = InterpType::CIC;

		     const Real weight  = 1.0;
		     const RealVect pos = RealVect::Zero;
		     m_initial_particles.add(Particle(weight, pos));
		  }

		  ~electron(){}

		  Real initial_data(const RealVect a_pos, const Real a_time) const {
		     return 1.0;
		  }
		};

The ``initial_data`` function sets the density to one everywhere. In addition, we have added a single particle with weight one at the Cartesian coordinates :math:`(x=0, y=0, z=0)`. Note that the two functions are additive. If you only want to use particles as initial data, you could either have ``initial_data`` return zero everywhere, or you can set the :ref:`Chap:species` class member ``m_init_with_function`` to ``false``. 

You may, in principle, add as many particles as you want. However, the particles are shared among all MPI ranks so there *is* a practical limit to how many you can use.

Defining transport boundary conditions
--------------------------------------

Transport boundary conditions are provided through the :ref:`plasma_kinetics` physics interface, please refer to that chapter for additional details.

Setting radiative transport boundary conditions
-----------------------------------------------

Boundary conditions for the radiative transer equations, if available, are set through the implementation classes. For example, for the Monte-Carlo module we have defined the following options:

.. code-block:: bash

   mc_photo.bc_x_low          = outflow       # Boundary condition. 'outflow', 'symmetry', or 'wall'
   mc_photo.bc_x_high         = outflow       # Boundary condition
   mc_photo.bc_y_low          = outflow       # Boundary condition
   mc_photo.bc_y_high         = outflow       # Boundary condition
   mc_photo.bc_z_low          = outflow       # Boundary condition
   mc_photo.bc_z_high         = outflow       # Boundary condition

For the diffusion-limited photon transport module (``eddington_sp1``), boundary conditions are always set through the following options:

.. code-block:: bash
		
   eddington_sp1.bc_x_low            = robin     # Boundary on domain. 'neumann' or 'robin'
   eddington_sp1.bc_x_high           = robin     # Boundary on domain. 'neumann' or 'robin'              
   eddington_sp1.bc_y_low            = robin     # Boundary on domain. 'neumann' or 'robin'
   eddington_sp1.bc_y_high           = robin     # Boundary on domain. 'neumann' or 'robin'
   eddington_sp1.bc_z_low            = robin     # Boundary on domain. 'neumann' or 'robin'
   eddington_sp1.bc_z_high           = robin     # Boundary on domain. 'neumann' or 'robin'		
