.. _Chap:NewSimulations:

Setting up simulations
======================

In `PlasmaC`, most of the action occurs through the physics interface :ref:`Chap:plasma_kinetics`.

.. _Chap:NewGeometry:

Setting up a geometry
---------------------

The creation of a geometry in `PlasmaC` consists of the implementation of a abstraction layer that describes the surface of electrodes and dielectrics, whic is done by implementing the :ref:`Chap:computational_geometry`. This class will ask the user to instantiate objects :ref:`Chap:electrode` and :ref:`Chap:dielectric`. Each of these objects contain an implicit function that describes the level-set surface of the boundary. If the implicit function does not exist in `Chombo` or in `PlasmaC`, the user will have to implement it himself using elements of computational geometry. Often, however, you will be able to re-use consisting implicit functions. These are located in ``src/geometry``. 

Creating a new geometry
_______________________

Creating a new geometry consists of the following steps:

1. Create a new inherited class of :ref:`Chap:computational_geometry`.
2. In the constructor, fill the internal data holders ``m_electrodes`` and ``m_dielectrics``, see :ref:`Chap:computational_geometry` for details.

Often, the geometry will be parametrized and the user will want to add an option file where the inherited class can fetch options at run-time. As an example, the code snippet below describes a complete implementation of an electrode geometry consisting of a single sphere at the center.

Fast geometry generation with pre-voxelization
______________________________________________

Geometry generation is a non-negligible part of `PlasmaC`. Many geometries consist of small object in large domains, and it is generally desireable to make geometry generation scale according to the size of the object rather than the size of the domain. To achieve this, we have introduced the concept of a *voxel*. A voxel consists of a Cartesian grid box described by its two corners, and a flag that describes the nature of that box. The user may provide three types of voxels to his :ref:`Chap:computational_geometry` class; covered voxels, regular voxels, and cut voxels. The default behavior of the pre-voxelization is as follows:

1. All grid boxes that intersect a cut voxel are load balanced among the available MPI ranks and queried for geometric information. 
2. All grid boxes that are completely contained inside a regular (or covered) voxel box is set to regular (or covered) *if* the regular voxel does not intersect a cut voxel. 

The advantage of using a pre-voxelization is that we can supply meta-information to the mesh generator so that all the parts of the domain that do not intersect an object are not actually queried for geometric information. This can lead to orders of magnitude improvements in the geometry generation step.

The user may set the voxels manually when he defines his geometry, or he may use our voxel-generating tool for defining them.

.. _Chap:GridResolutions:

Setting the grid resolution
---------------------------

Setting the grid resolution in `PlasmaC` is done by selecting a base domain for :ref:`Chap:amr_mesh`, a blocking factor (i.e. the minimum box size), the maximum number of refinement levels, and the refinement factor between each level. The *base domains* is the coarsest grid level on which all equations are solved (there are special rules for the multigrid).

Depending on how many cores is used in a simulation, there might be reasons to select quite large base levels. Typically, one should aim to have at least one grid patch on the base level per MPI rank. 

Selecting a base domain
_______________________

Usually, the blocking factor is 8, 16, or 32. Note that `PlasmaC` uses a patch-based formalism and domain sizes are therefore not limited to the conventional :math:`2^N` grids encountered in octrees. In `PlasmaC`, the base level must instead be divisible by the blocking factor along each coordinate direction. For example, you may have a domain which is :math:`112^3` cells with a blocking factor of :math:`16`, but you may *not* use a domain size of :math:`112\times112\times56` since :math:`56` is not (integer) divisible by 16. In this case, the solution would be to use a slightly higher resolution and using a domain size of :math:`128\times128\times64` instead.

The base domain is selected by setting the following run-time parameters:

.. code-block:: bash

		amr_mesh.coarsest_domain = 128 128 128 # Number of cells on coarsest domain
		amr_mesh.max_box_size    = 16          # Maximum allowed box size

Setting the grid refinement
___________________________

The grid refinement is set by selecting a number of refinement levels *and* the refinement factor between each level. Refinement factors of 2 and 4 are supported. Note that using refinement factors of 4 reduce the number of AMR levels, but lead to additional noise in the solutions and, furthermore, usually leads to many more grid cells.

Setting the grid refinement is done through run-time parameters in :ref:`Chap:amr_mesh`:

.. code-block:: bash

		amr_mesh.max_amr_depth   = 4           # Maximum amr depth
		amr_mesh.ref_rat         = 4 2 4 2 4 2 # Refinement ratios


Multigrid coarsening
____________________

The multigrid solvers in `PlasmaC` use grids that are coarser than the base grid in order to facilitate better convergence. In multigrid, the base grid is usually coarsened as far as possible, but the user also has the option to call the bottom solver on a coarser grid level. In `PlasmaC`, this is achieved by first coarsening the base level by a factor of two as far as possible, while keeping the grid coarsenable by the blocking factor. For the example of :math:`128\times128\times64` grid cells with blocking factor 16, we may reach levels down to :math:`4\times4\times2`. For the example of :math:`112\times112\times112` cells the base grid is divided into :math:`7\times7\times7` boxes with a blocking factor of :math:`16`. Coarsening by a factor of two is not possible since the domain :math:`56\times56\times56` is not divisible by a blocking factor of :math:`16`. Instead, the multigrid coarsening routines is allowed to use smaller box sizes by maintaing the :math:`7\times7\times7` box structure and then coarsening the boxes themselves. I.e. the :math:`56\times56\times56` domain is obtained by a :math:`7\times7\times7` box decomposition with a box size of :math:`8`. The next multigrid level uses a box size of :math:`4`, and so on. 

Conflicts may occur if the user attempts to exit multigrid on a coarsened level that does not exist. For the above example of :math:`112\times112\times112` cells the coarsened multigrid levels are :math:`56^3` and :math:`28^3`, so if the user attempts to call the multigrid bottom solver at a coarsening of :math:`32\times32\times32` cells, a run-time error will occur since this level cannot be reached by standard grid coarsening procedures. 
   
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

Implicit plasma chemistry
_________________________

In cases where transport and plasma chemistry is split, implicit treatment of the plasma chemistry terms is possible. However, it is not natively supported. The reason for this is that the plasma chemistry terms can be non-local in space, and even stochastic. If the user wants to use implicit chemistry, he will have to implement it himself.

In `PlasmaC` the plasma chemistry is always advanced through a routine

.. code-block:: c++

  virtual void advance_reaction_network(Vector<Real>&          a_cdr_sources,
					Vector<Real>&          a_photon_sources,
					const Vector<Real>     a_cdr_densities,
					const Vector<RealVect> a_cdr_gradients,
					const Vector<Real>     a_photon_densities,
					const RealVect         a_E,
					const RealVect         a_pos,
					const Real             a_dx,
					const Real             a_dt,
					const Real             a_time,
					const Real             a_kappa) const = 0;		

and the assumption is that this routine will provide source terms for the convection-diffusion-reaction solvers and the radiative transport solvers for advancement over a time step ``a_dt``. In all of `PlasmaC` this routine is used such that the plasma chemistry is *always* implies the advance

.. math::

   \phi^{k+1} = \phi^{k} + \Delta t S,

where :math:`\phi^k` is ``a_cdr_densities`` in the function call above and :math:`S` is the output argument ``a_cdr_sources`` in the ``advance_reaction_network`` routine above. However, we make no assumptions about how :math:`S` is computed. Usually, :math:`S` is computed in some explicit form using tabulated values for ionization coefficients or somesuch, and the above equation becomes a forward Euler method. This is the assumption that we make in e.g. the ``imex_sdc`` class. However, one may certainly perform an implicit advance over the time step ``a_dt`` inside the ``advance_reaction_network`` call, and then set the source term as :math:`S = (\phi^{k+1}-\phi^{k})/\Delta t`. This is perfectly consistent will all the `PlasmaC` integrators and it implies that each plasma chemistry update is done using the internals of ``advance_reaction_network``. 

As an example, consider that one wants to advance :math:`\partial_t\phi = \alpha\phi` implicitly by using the backward Euler method. The solution is :math:`\phi^{k+1} = \phi^k/(1-\alpha\Delta t)` and :math:`S = \frac{\alpha}{1-\alpha \Delta t}\phi^k`, although this latter step would simply be done numerically using :math:`S = \left(\phi^{k+1}-\phi^k\right)/\Delta t`, as implemented below:

.. code-block:: c++
		
  virtual void advance_reaction_network(Vector<Real>&          a_cdr_sources,
					Vector<Real>&          a_photon_sources,
					const Vector<Real>     a_cdr_densities,
					const Vector<RealVect> a_cdr_gradients,
					const Vector<Real>     a_photon_densities,
					const RealVect         a_E,
					const RealVect         a_pos,
					const Real             a_dx,
					const Real             a_dt,
					const Real             a_time,
					const Real             a_kappa) const{

     Real phiOld = a_cdr_densities[0];
     Real phiNew = phiOld/(1-alpha*a_dt);
     a_cdr_sources[0] = (phiNew - phiOld)/a_dt
  }					 

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
