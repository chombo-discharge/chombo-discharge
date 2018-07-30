.. _Chap:plasma_kinetics:

plasma_kinetics
---------------

:ref:`Chap:plasma_kinetics` is the physics module of PlasmaC. The entire class is an interface, whose implementations run deep into :ref:`Chap:time_stepper` which contains much of the low-level functionality that couples the schemes. See the :doxy:`Doxygen API <plasma_kinetics>` for details. By far, :ref:`Chap:plasma_kinetics` is the most time-consuming tasks of implementing new plasma schemes. The reason for this is that we need to maintain a certain level of abstraction in order to cover a broad spectrum of plasma phenomena. 

There are no default input parameters for :ref:`Chap:plasma_kinetics`, as users must generally implement their own kinetics. A successful implementation of :ref:`Chap:plasma_kinetics` has the following:

* Instantiated :ref:`Chap:species`
* Instantiated :ref:`Chap:photon`
* Implemented the core functionality

Based on these three points, PlasmaC will automatically allocate the specified numbers of convection-diffusion-reaction and radiative transport solvers. For information on how to interface into the CDR solvers, see :ref:`Chap:species`. Likewise, see :ref:`Chap:photon` for how to interface into the RTE solvers. 

Implementation of the core functionality is comparatively straightforward. In the constructor, the user should fetch his input parameters (if he has any) and **must** instantiate all species and photons in the internal containers. When the class is used by PlasmaC later on, all arguments in the core functions follow that ordering. The core functionality is given by the following functions: 

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
