.. _Chap:MinimalPlasmaModel:

Minimal plasma model
====================

The minimal plasma model resides in :file:`/physics/cdr_plasma` and describes plasmas in the drift-diffusion approximation.
This physics model also includes the following subfolders:

* :file:`/physics/cdr_plasma/plasma_models` which contains various implementation of some plasma models that we have used.
* :file:`/physics/cdr_plasma/time_steppers` contains various algorithms for advancing the equations in an EBAMR context.
* :file:`/physics/cdr_plasma/cell_taggers` contains various algorithms for advancing the equations in an EBAMR context.
* :file:`/physics/cdr_plasma/python` contains Python source files for quick setup of "mini-applications".

If users decide to develop new code for the minimal plasma model, e.g. new grid refinement routines, plasma models, integrators, or other types of improvements, they should be put in the folders above. 

Equations of motion
-------------------

In the minimal plasma model, we are solving

.. math::
   :nowrap:

   \begin{align}
   &\nabla\cdot\left(\epsilon_r\nabla\Phi\right) = -\frac{\rho}{\epsilon_0}, \\[1ex]
   &\frac{\partial\sigma}{\partial t} = F_\sigma,\\[1ex]
   &\frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n + \sqrt{2D\phi}\mathbf{Z}\right) = S,
   \end{align}

where :math:`\sqrt{2D\phi}\mathbf{Z}` is a stochastic diffusion flux suitable for fluctuating hydrodynamics models.
By default, this flux is turned off. 

The above equations must be supported by additional boundary conditions on electrodes and insulating surfaces. 

Radiative transport is also supported, which is done either in the diffusive approximation or by means of Monte Carlo methods. Diffusive RTE methods involve solving

.. math::
   :nowrap:

   \begin{align}
      \partial_t\Psi + \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) &= \frac{\eta}{c},
   \end{align}
   
where :math:`\Psi` is the isotropic photon density, :math:`\kappa` is an absorption length and :math:`\eta` is an isotropic source term.
The time dependent term can be turned off and the equations can be solved stationary.
As an alternative, we also provide discrete photon methods that solve for the photoionization profile on a mesh by sampling discrete photons.
Our discrete photon methods are capable of including far more physics; they can easily be adapted to e.g. scattering media and also provide much better qualitative features (like shadows, for example).
They are, on the other hand, inherently stochastic which implies that some extra care must be taken when integrating the equations of motion.

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

`PlasmaC` works by embedding the equations above into an abstract C++ framework that the user must implement or reuse existing pieces of, and then compile into a *mini-application*.

.. _Chap:cdr_plasma_physics:

cdr_plasma_physics
------------------

:ref:`Chap:cdr_plasma_physics` represents the microphysics in :file:`/physics/cdr_plasma`.
The entire class is an interface, whose implementations are used in the time integrators that advance the equations.
As mentioned above, the time integrators are located in :file:`/physics/cdr_plasma/time_steppers`. 

There are no default input parameters for :ref:`Chap:cdr_plasma_physics`, as users must generally implement their own kinetics.
A successful implementation of :ref:`Chap:cdr_plasma_physics` has the following:

* Instantiated :ref:`Chap:species`. These contain metadata for the transport solvers. 
* Instantiated :ref:`Chap:photon`. These contain metadata for the radiative transport solvers. 
* Implemented the core functionality that couple all solvers.

`PlasmaC` automatically allocates the specified number of convection-diffusion-reaction and radiative transport solvers. For information on how to interface into the CDR solvers, see :ref:`Chap:species`. Likewise, see :ref:`Chap:photon` for how to interface into the RTE solvers.

There are currently no support for existing file formats for describing reactions and so on. If you have a huge list of reactions that need to be implemented, it would probably pay off to write a code-generating interface between :ref:`Chap:cdr_plasma_physics` and your list of reactions. 

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

The above code blocks do exactly what their signatures indicate. It is up to the user to implement these. The length of the Vector holding the return values from these functions are expected to be equal to the number of CDR solvers, with the exception of *compute_rte_source_terms* which has the length given by the number of RTE solvers. Note that in all of the above, the ordering of the input vectors are expected to be the same as the ordering of the species vector of :ref:`Chap:cdr_plasma_physics`. 

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

The :ref:`Chap:species` is a lightweight class used to provide information into convection-diffusion-reaction solvers. This class is mostly used within :ref:`Chap:cdr_plasma_physics` in order to provide information on how to instantiate CDR solvers. :ref:`Chap:species` is abstract so that the user must implement

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



The members ``m_mobile`` and ``m_diffusive`` are used for optimization in `PlasmaC`: If the user specifies that a species is immobile, `PlasmaC` will skip the advection computation. Note that ``m_diffusive`` and ``m_mobile`` override the specifications in :ref:`Chap:cdr_plasma_physics`. If the user provides a non-zero velocity through :ref:`Chap:cdr_plasma_physics` function *compute_cdr_velocities*, and sets ``m_mobile`` to ``false``, the species velocity will be zero. Of course, the user will often want to provide additional input information to his species, for example by specifying a seed for the initial conditions. 

.. _Chap:photon:

photons
_______

:ref:`Chap:photon` is the class that supplies extra information to the RTE solvers. In those solvers, the source term computation is handled by :ref:`Chap:cdr_plasma_physics`, so the :ref:`Chap:photon` class is very lightweight. The user must implement a single function which specifies the absorption coefficient at a point in space:

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

By default, there are no input parameters available for the :ref:`Chap:photon` class, but the user will often want to include these, for example by modifying the absorption coefficient. Note that you are allowed to use a spatially varying absorption coefficient. 

For most users, this will mostly include implementing a new geometry or a new plasma-kinetic scheme.
It is possible to generate entirely new physics interfaces, too.
Our goal is that the user does not need to worry about temporal or spatial discretization of these equations, but rather focus on the actual setup of the geometry and physics. 


Temporal discretization
-----------------------

In this chapter we discuss the supported temporal integrators for `PlasmaC`, and discuss their input parameters. These integrators differ in their level of efficiency and accuracy. Currently, none of the integrators can subcycle in time.

Time step limitations
---------------------

Our time integrators have different time step limitations. Fully explicit codes are limited by the advective and diffusive CFL constraints and usually also the dielectric relaxation time. However, some of the `PlasmaC` integrators eliminate the dielectric relaxation time, and all integrators can handle diffusion either implicitly or explicitly. In some cases only the advective CFL constraint is the only time step restriction.

  
Deterministic integrators
-------------------------

.. _Chap:godunov:

godunov
_______
The ``godunov`` temporal integrator is a rather unsophisticated, but very stable, temporal integrator. ``godunov`` uses an operator splitting between charge transport and plasma chemistry in the following way:

1. Advance :math:`\partial_t\phi = -\nabla\cdot\left(\mathbf{v}\phi - D\nabla\phi + \sqrt{2D\phi}\mathbf{Z}\right)` (and the surface charge solver) over a time step :math:`\Delta t`.
2. Compute the electric field 
3. Advance the plasma chemistry over the same time step using the field computed in 2). I.e. advance :math:`\partial_t\phi = S` over a time step :math:`\Delta t`.   
4. Move photons and deposit them on the mesh. 
      
Various integration options for the transport and chemistry steps are available but are discussed elsewhere. Note that the ``godunov`` integrator uses a semi-implicit coupling between the plasma chemistry terms and the electric field, and therefore eliminates the so-called dielectric relaxation time. The formal order of convergence of the ``godunov`` integrator is 1, but the accuracy can be quite good depending on the transport and chemistry schemes that are chosen.

.. _Chap:SISDC:

imex sdc
________
``imex_sdc`` is a semi-implicit spectral deferred correction method for the `PlasmaC` equation set and is an adaptive high-order discretization with implicit diffusion. This method integrates the advection-diffusion-reaction equations in the following way.

Spectral deferred corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given an ordinary differential equation (ODE) as

.. math::
   \frac{\partial u}{\partial t} = F(u,t), \quad u(t_0) = u_0,

the exact solution is
  
.. math::
   u(t) = u_0 + \int_{t_0}^tF\left(u,\tau\right)d\tau.
   
Denote an approximation to this solution by :math:`\widetilde{u}(t)` and the correction by :math:`\delta(t) = u(t) - \widetilde{u}(t)`. The measure of error in :math:`\widetilde{u}(t)` is then

.. math::
   R(\widetilde{u}, t) = u_0 + \int_{t_0}^tF(\widetilde{u}, \tau)d\tau - \widetilde{u}(t).

Equivalently, since :math:`u = \widetilde{u} + \delta`, we can write

.. math::
   \widetilde{u} + \delta = u_0 + \int_{t_0}^t F\left(\widetilde{u}+\delta, \tau\right)d\tau. 

This yields

.. math::
   \delta = \int_{t_0}^t\left[F\left(\widetilde{u}+\delta, \tau\right) - F\left(\widetilde{u}, \tau\right)\right]d\tau + R\left(\widetilde{u},t\right). 

This is called the correction equation. The goal of SDC is to iteratively solve this equation in order to provide a high-order discretization. 

We now discuss the semi-implicit SDC (SISDC) method. First, we apply the method of lines (MOL) such that

.. math::
   :nowrap:
   
   \begin{align}
   \frac{d\phi_{\mathbf{i}}}{dt} &= \mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}\right) + \mathcal{F}_{\textrm{D}}\left(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}\right), \\
   \frac{d\sigma_{\mathbf{i}}}{dt} &= \mathcal{F}_{\sigma}\left(t, \phi_{\mathbf{i}}\right),
   \end{align}

where :math:`\phi_{\mathbf{i}}` denotes a cell-averaged variable, :math:`\mathcal{F}_{\sigma}` is as described in :ref:`Chap:SpatialDiscretization`, :math:`\mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}\right) = -D^c_{\mathbf{i}} + S_{\mathbf{i}}` is the advection-reaction operator , and :math:`\mathcal{F}_{\textrm{D}}(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}) = \frac{1}{\kappa_{\mathbf{i}}}\int_{V_{\mathbf{i}}}\left[\nabla\cdot\left(D\nabla\phi\right)\right]dV_{\mathbf{i}}` is the diffusion operator. Note that the advective operator contains the hybrid divergence discussed in :ref:`Chap:CDR` and :math:`\mathcal{F}_{\textrm{D}}` is parametrically coupled to :math:`\mathbf{E}` through :math:`D = D\left(\mathbf{E}\right)` (we use a semi-colon to indicate this dependence). Strictly speaking, :math:`\mathcal{F}_{\textrm{AR}}` is parametrically coupled in the same way through the  mobilities and boundary conditions, and additionally coupled to :math:`\Psi` through source terms so that the notation :math:`\mathcal{F}_{\textrm{AR}}\left(t, \phi_{\mathbf{i}}; \mathbf{E}_{\mathbf{i}}, \Psi_{\mathbf{i}}\right)` would be appropriate. However, charge injection, advection, and chemistry will be integrated explicitly so this dependence is notationally suppressed. On the other hand, the diffusion part will be solved with the backward Euler method - which yields a Helmholtz equation - and so we need to maintain this dependence for now. Later, we will clarify how this dependence is resolved. The rationale for solving diffusion implicitly is due to the numerical time step constraint of explicit diffusion methods which scales as :math:`\mathcal{O}\left(\Delta x^2\right)`, whereas advection scales more favorably at :math:`\mathcal{O}\left(\Delta x\right)`. We have chosen to integrate the reactive terms explicitly. The reason is that the reactive terms can be non-local, i.e. they can depend on the electron gradient. This is for example the case for fluid models in the local energy approximation where the electron energy source term contains terms that are proportional to the electron diffusion term :math:`D_e\nabla\phi_e`. Implicit discretization of the reactive terms then yield a fully coupled system rather than systems coupled only within individual cells. Charge injection is also handled explicitly. This design choice is mandated by the fact that implicit charge injection through the diffusion terms couples every diffusive species, leading to a system of diffusion equations that are fully coupled through their boundary conditions. Although charge injection could reasonably be performed as a separate step, this leads to numerical instabilities for cut-cell methods since the injected charge must be normalized by the volume fraction of the cell (which can be arbitrarily small). 

SISDC predictor
^^^^^^^^^^^^^^^

Next, we present the SISDC method. In what follows, we suppress the index :math:`{\mathbf{i}}` as it is not explicitly needed. Given an interval :math:`[t_n, t_{n+1}]` on which a solution is sought, SDC methods divide this interval into :math:`p` subintervals :math:`t_n = t_{n,0} < t_{n,1} < \ldots < t_{n,p} = t_{n+1}`. Our discussion, however, pertains only to the interval :math:`[t_n, t_{n+1}]` so we compress the notation to :math:`t_m\equiv t_{n,m}`. We obtain an initial solution :math:`\phi_{m}^0, m=0,1,\ldots,p` as the semi-implicit advance

.. math::
   :nowrap:

   \begin{align}
   \phi_{m+1}^0 &= \phi_m^0 + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m,\phi_m^0\right) + \mathcal{F}_{\textrm{D}}\left(t_{m+1},\phi_{m+1}^0; \mathbf{E}_{m+1}^0\right)\right],\\
   \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m,\phi_m^0\right).
   \end{align}

This defines a Helmholtz problem for :math:`\phi_{m+1}^0` through :math:`\mathcal{F}_{\textrm{D}}`. Generally, the upper subscript denotes an SDC iteration where subscript 0 is the SISDC predictor, and we also have :math:`\phi_0^0 = \phi(t_n)` and :math:`\sigma_0^0 = \sigma(t_n)`. This predictor treats advection and chemistry terms explicitly, and diffusion implicitly. Other types of semi-implicit or multi-implicit couplings are possible :cite:`Bourlioux2003,Layton2004,Nonaka2012`. SDC improves this solution by using deferred corrections: Given a numerical solution :math:`\phi_{m+1}^k`, we compute an error :math:`\delta_{m+1}^k` and obtain the next iterate :math:`\phi_{m+1}^{k+1} = \phi_{m+1}^k + \delta_{m+1}^k`. Each iteration raises the discretization order by one :cite:`Dutt2000,Minion2003`, to maximum order :math:`p+1`. Critical to the success of this approach is the precise evaluation of the numerical quadrature. 

The parametric coupling of the electric field complicates things since the predictor contains :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0\right)`, implying that the Poisson equation and the diffusion advance require concurrent solves for the diffusion update. We simplify this system by using a weak coupling by first computing

.. math::
   :nowrap:
   
   \begin{align}
   \phi_{m+1}^{0,\ast} &= \phi_m^0 + \Delta t_m\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^0\right), \\
   \sigma_{m+1}^0 &= \sigma_m^0 + \Delta t_mF_\sigma\left(t_m, \phi_m^0\right),
   \end{align}

Next, we will approximate :math:`\mathbf{E}_{m+1}^{0}` for use in the predictor. There are two choices for this coupling; one may either use :math:`\mathbf{E}_m^0` for computation of the diffusion coefficients, which we will refer to as the semi-implicit coupling, or one may use fixed-point iteration and compute :math:`\mathbf{E}_{m+1}^{0,\ast} = \mathbf{E}\left(\phi_{m+1}^{0, \ast}, \sigma_{m+1}^0\right)`, followed by the diffusion advance

.. math::
   \phi_{m+1}^{0,\dagger} = \phi_{m+1}^{0,\ast} + \Delta t_m\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^0; \mathbf{E}_{m+1}^\ast\right),

which we will refer to as the implicit coupling. This is e.g. the electric field coupling used in :cite:`Marskar2019`. This approximation can be improved by using more fixed-point iterations that computes :math:`\mathbf{E}_{m+1}^{0,\dagger} = \mathbf{E}\left(\phi_{m+1}^{0,\dagger}, \sigma_{m+1}^0\right)` and then re-solves the predictor equation with :math:`\mathbf{E}_{m+1}^{0,\dagger}` in place of :math:`\mathbf{E}_{m+1}^{0,\ast}`. The process can then be repeated for increased accuracy. Regardless of which coupling is used, we have now calculated :math:`\phi_{m+1}^0`, :math:`\sigma_{m+1}^0`, through which we obtain :math:`\mathbf{E}_{m+1}^0 = \mathbf{E}\left(\phi_{m+1}^0, \sigma_{m+1}^0\right)`, and :math:`\Psi_{m+1}^0 = \Psi\left(\mathbf{E}_{m+1}^0, \phi_{m+1}^0\right)`. Finally, we remark that the SISDC predictor is a sequentially advanced semi-implicit Euler method, which is locally second order accurate and globally first order accurate. Each step of the predictor can be thought of as a Godunov splitting between the advective-reactive and diffusive terms. 

SISDC corrector
^^^^^^^^^^^^^^^
Next, the semi-implicit discretization of the correction equation is

.. math::
   \begin{split}
   \delta_{m+1}^k &= \delta_m^k  + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k + \delta_m^k\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)\right.\\
   &+ \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k + \delta_{m+1}^k; \mathbf{E}_{m+1}^k\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k; \mathbf{E}_{m+1}^k\right)\right] - \left(R_{m+1}^k - R_{m}^k\right).
   \end{split}

We furthermore define

.. math::
   \begin{split}
   R_{m+1}^k - R_m^k &= \int_{t_m}^{t_{m+1}}\left[\mathcal{F}_{\textrm{AR}}\left(\phi^k\right) + \mathcal{F}_{\textrm{D}}\left(\phi^k; \mathbf{E}^k\right)\right]d\tau - \phi_{m+1}^k + \phi_m^k \\
   &\equiv I_m^{m+1}\left(\phi^k\right) - \phi_{m+1}^k + \phi_m^k. 
   \end{split} 

Evaluation of :math:`I_m^{m+1}` yields :math:`p` quadrature rules and we may write

.. math::
   I_m^{m+1}\left(\phi^k\right) = \sum_{l=0}^p q_m^l\left[\mathcal{F}_{\textrm{AR}}\left(t_l, \phi^k_l\right) + \mathcal{F}_{\textrm{D}}\left(t_l, \phi^k_l; \mathbf{E}_l^k\right)\right],

where the weights :math:`q_m^l` are quadrature weights. The final update for :math:`\phi^{k+1}_{m+1}` is then

.. math::
   \begin{split}
   \phi_{m+1}^{k+1} &= \phi_{m}^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) -\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right.\\
   & + \left.\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \phi_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right] + I_{m}^{m+1}\left(\phi^k\right).
   \end{split}

With the exception of :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right)`, all quantities on the right-hand are known and the correction equation is reduced to a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` with error :math:`\delta_{m+1}^k = \phi_{m+1}^{k+1} - \phi_{m+1}^k`. An analogous equation is found for :math:`\sigma_{m+1}^{k+1}`.

The correction step has the same coupling to the electric field as the prediction step in that :math:`\mathbf{E}_{m+1}^{k+1}` appears in the update equation for :math:`\phi_{m+1}^{k+1}`. As for the prediction, we use a weak coupling through which we first compute

.. math::
   :nowrap:
   
   \begin{align}
   \phi_{m+1}^{k+1,\ast} &= \phi_m^{k+1} + \Delta t_m\left[\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k+1}\right) - \mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^{k}\right)\right] + I_m^{m+1}\left(\phi^k\right),\\
   \sigma_{m+1}^{k+1} &= \sigma_m^{k+1} + \Delta t_m\left[F_\sigma\left(t_m, \phi_m^{k+1}\right) - F_\sigma\left(t_m, \phi_m^{k}\right)\right] + \Sigma_m^{m+1}\left(\phi^k\right). 
   \end{align}

The solution for :math:`\sigma_{m+1}^{k+1}` is final since all charge is injected through the advection operator for :math:`\phi`. The term :math:`\Sigma_m^{m+1}` contains the injected charge through :math:`I_m^{m+1}\left(\phi^k\right)`, as was discussed in :ref:`Chap:SpatialDiscretization`. We then solve

.. math::
   \phi_{m+1}^{k+1} = \phi_{m+1}^{k+1, \ast} + \Delta t_m\left[\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k+1}; \mathbf{E}_{m+1}^{k+1}\right) - \mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)\right],

with some approximation for :math:`\mathbf{E}_{m+1}^{k+1}`. As before, this coupling can be made either semi-implicitly or implicitly. The corrector equation defines a Helmholtz equation for :math:`\phi_{m+1}^{k+1}` using :math:`\phi_{m+1}^{k+1,\ast}` as the previous solution and :math:`-\mathcal{F}_{\textrm{D}}\left(\phi_{m+1}^{k}; \mathbf{E}_{m+1}^k\right)` as a source term.

Order, stability, and computational cost
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For consistency with the literature, denote the SISDC method which uses :math:`P` nodes (i.e. :math:`P-1` subintervals) and :math:`K` total iterations (i.e. :math:`K-1` iterations of the correction equation) by :math:`\verb|SISDC|_P^K`. This method will have a global order of accuracy :math:`\min\left(K,P\right)` if the quadrature can be evaluated with appropriate accuracy. Order reductions may occur if the interpolating polynomial in the quadrature suffers from Runge's phenomenon. As we discuss below, uniformly spaced nodes have some computational advantage but is therefore also associated with some risk. Safer choices include Lobatto nodes or Chebyshev nodes (with inclusion of endpoints) to minimize the risk of order reductions. Implications on the choice of quadrature nodes can be found in :cite:`Layton2005`. 

For explicit advection, the deferred correction procedure integrates the correction equation sequentially and therefore does not allow each substep :math:`\Delta t_m` to exceed the CFL-limited time step :math:`\Delta t_{\textrm{cfl}}`, i.e. :math:`\Delta t_m < \Delta t_{\textrm{cfl}} \forall m`. Since we have :math:`\Delta t = \sum_m\Delta t_m`, uniform nodes maximize :math:`\Delta t` subject to the CFL constraint. For example, an :math:`\verb|SISDC|_P^K` method with uniformly spaced nodes has a maximum possible time step :math:`\Delta t < (P-1)\Delta t_{\textrm{cfl}}`. For the same number of function evaluations, the allowed time step with for Lobatto or Chebyshev nodes is smaller. For :math:`P\leq 3`, the uniform nodes, Lobatto nodes, and Chebyshev nodes coincide. Larger time steps are possible with uniform nodes for :math:`P>3`, which has some computational consequence. The table below summarizes the largest possible time steps for the :math:`\verb|SISDC|_P^K` method with the various quadratures. Finally, note that :math:`\Delta t_m < \Delta t_{\textrm{cfl}}` does not guarantee stability since further restrictions may required for stability of the reaction terms.

==========  =================================== ====================================   ================================
 :math:`P`   Lobatto                             Chebyshev                             Uniform
==========  =================================== ====================================   ================================
2           :math:`\Delta t_{\textrm{cfl}}`      :math:`\Delta t_{\textrm{cfl}}`       :math:`\Delta t_{\textrm{cfl}}`
3           :math:`2\Delta t_{\textrm{cfl}}`     :math:`2\Delta t_{\textrm{cfl}}`      :math:`2\Delta t_{\textrm{cfl}}`
4           :math:`2.26\Delta t_{\textrm{cfl}}`  :math:`1.73\Delta t_{\textrm{cfl}}`   :math:`3\Delta t_{\textrm{cfl}}`
5           :math:`3.05\Delta t_{\textrm{cfl}}`  :math:`2.82\Delta t_{\textrm{cfl}}`   :math:`4\Delta t_{\textrm{cfl}}`
6           :math:`3.50\Delta t_{\textrm{cfl}}`  :math:`3.29\Delta t_{\textrm{cfl}}`   :math:`5\Delta t_{\textrm{cfl}}`
7           :math:`4.26\Delta t_{\textrm{cfl}}`  :math:`4.36\Delta t_{\textrm{cfl}}`   :math:`6\Delta t_{\textrm{cfl}}`
==========  =================================== ====================================   ================================

For the predictor step, it is necessary to evaluate :math:`\mathcal{F}_{\textrm{AR}}\left(\phi_m^{k+1}\right)` and thus update the Poisson and radiative transfer equations at each node. In addition, it is necessary to solve the diffusion equation at every node except :math:`m=0` for every diffusive species, which may also require auxiliary updates of the electric field. The corrector step contains extra floating point operator due to the extra terms :math:`\mathcal{F}_{\textrm{AR}}\left(t_m, \phi_m^k\right)` and :math:`\mathcal{F}_{\textrm{D}}\left(t_{m+1}, \phi_{m+1}^k\right)` and the quadrature :math:`I_m^{m+1}`. The computational cost of adding in these terms is small compared to the cost of an Euler update of the advection-reaction equation since one must also computate source terms, drift velocities, and boundary conditions in addition to construction of the hybrid divergence. In short, the computational cost of the predictor and corrector steps are about the same.

Next, we provide some remarks on the extra computational work involved for higher order methods. Broadly speaking, the total amount of floating point operations increases quadratically with the order. Each node requires evaluation of one advection-reaction operator, at least one electric field update, and one radiative transfer update. Likewise, each substep requires one diffusion solve. Thus, :math:`\verb|SISDC|_K^K` requires :math:`K^2` advection-reaction evaluations, :math:`(K-1)^2` diffusion solves, :math:`(K-1)^2` radiative transfer updates, and at least :math:`K^2` electric field updates. In these estimates we have assumed that the diffusion solve couples semi-implicitly to the electric field, thus each corrector iteration requires one electric field update per node, giving a total cost :math:`K^2`. Strictly speaking, the number of advection-reaction evaluations is slightly less since :math:`\mathcal{F}_{\textrm{AR}}\left(t_0, \phi_0^k\right)` does not require re-evaluation in the corrector, and :math:`\mathcal{F}_{\textrm{AR}}\left(t_p,\phi_p^{K-1}\right)` does not need to be computed for the final iteration since the lagged quadrature is not further needed. Nonetheless, the computational work is quadratically increasing, but this is partially compensated by allowance of larger time steps since the :math:`\verb|SISDC|_K^K` has a stability limit of :math:`(K-1)\Delta t_{\textrm{cfl}}` rather than :math:`\Delta t_{\textrm{cfl}}` for uniformly spaced nodes. For comparison with the predictor :math:`\verb|SISDC|_K^1` which is a first order method, the work done for integration over :math:`(K-1)\Delta t_{\textrm{cfl}}` amounts to :math:`K-1` advection-reaction updates, :math:`K-1` diffusion updates, :math:`K-1` radiative transfer updates, and :math:`K` electric field updates. If we take the electric field updates as a reasonable metric for the computational work, the efficiency of the :math:`K` th order method over the first order method is about :math:`K` for integration over the same time interval, i.e. it increases linearly rather than quadratically. However, this estimate is only valid if we do not take accuracy into account. In practice, the predictor does not provide the same accuracy as the corrector over the same integration interval. A fair comparison of the extra computational work involved would require that the accuracy of the two methods be the same after integration over a time :math:`(K-1)\Delta t_{\textrm{cfl}}`, which will generally require more substeps for the first order method. While we do not further pursue this quantification in this paper, the pertinent point is that the extra computational work involved for tolerance-bound higher order discretizations increases sub-linearly rather than quadratically when compared to lower-order equivalents.

We have implemented the SISDC algorithm in the ``sisdc`` class in :file:`/time_steppers/sisdc`. The following class options are available:

.. literalinclude:: links/imex_sdc.options

.. _Chap:MISDC:


Stochastic integrators
----------------------

Deterministic CFD integrators are generally not suitable for stochastic ODEs. For example, the recursive nature of Heun's method or spectral deferred corrections hardly make sense for stochastic ODEs.

For fluctuating hydrodynamics we are preparing several temporal integrators. ``godunov`` is generally useful for FHD but we also provide an Euler-Maruyama integrator (``euler_maruyama``) which is first order accurate in time (although with an accurate advective integrator). Although ``euler_maruyama`` is functionally similar to ``godunov``, it does not eliminate the dielectric relaxation time and is therefore less stable for some simulation cases. 

.. _Chap:euler_maruyama:

euler_maruyama
______________

:ref:`Chap:euler_maruyama` implements the Euler-Maruyama method. This method is based on an Euler method with explicit or implicit diffusion. 
