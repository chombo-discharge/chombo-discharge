/*!
  @file plasma_kinetics.H
  @brief Declarationof the plasma kinetics class
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PLASMA_KINETICS_
#define _PLASMA_KINETICS_

#include <RealVect.H>
#include <RefCountedPtr.H>

#include "species.H"
#include "photon_group.H"

/*!
  @brief Abstract base class for specifying plasma kinetics
*/
class plasma_kinetics {
public:

  /*!
    @brief Constructor
  */
  plasma_kinetics(){}

  /*!
    @brief Destructor
  */
  virtual ~plasma_kinetics(){}

  /*!
    @brief Compute velocities for the CDR equations
    @param[in] a_pos Position
    @param[in] a_time Time
    @param[in] a_E Electric field
    @param[in] a_cdr_densites CDR densities
  */
  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const = 0;

    /*!
    @brief Compute the cdr diffusion coefficients
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_E Electric field
    @param[in] a_cdr_densites CDR densities
  */
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const = 0;

  /*!
    @brief Compute source terms for the cdr equations
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_E Electric field
    @param[in] a_gradE \f$\nabla|\mathbf{E}|\f$
    @param[in] a_cdr_densities CDR densities
    @param[in] a_rte_densities RTE densities (isotropic terms)
    @param[in] a_grad_cdr Gradients of cdr densities
  */
  virtual Vector<Real> compute_cdr_source_terms(const Real              a_time,
						const RealVect&         a_pos,
						const RealVect&         a_E,
						const RealVect&         a_gradE,
						const Vector<Real>&     a_cdr_densities,
						const Vector<Real>&     a_rte_densities,
						const Vector<RealVect>& a_grad_cdr) const = 0;
  /*!
    @brief Compute cdr fluxes through electrode-gas interfaces
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_normal Normal vector, points into the gas phase. 
    @param[in] a_E Electric field
    @param[in] a_cdr_velocities CDR velocities. Normal component only. 
    @param[in] a_cdr_densities CDR densities. 
    @param[in] a_rte_fluxes RTE fluxes (normal component only)
    @param[in] a_extrap_cdr_fluxes Extrapolated fluxes from the gas side. 
  */
  virtual Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
						    const RealVect&     a_pos,
						    const RealVect&     a_normal,
						    const RealVect&     a_E,
						    const Vector<Real>& a_cdr_densities,
						    const Vector<Real>& a_cdr_velocities,
						    const Vector<Real>& a_rte_fluxes,
						    const Vector<Real>& a_extrap_cdr_fluxes) const = 0;

  /*!
    @brief Compute cdr fluxes through dielectric-gas interfaces
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_normal Normal vector, points into the gas phase. 
    @param[in] a_E Electric field
    @param[in] a_cdr_velocities CDR velocities. Normal component only. 
    @param[in] a_cdr_densities CDR densities. 
    @param[in] a_rte_fluxes RTE fluxes (normal component only)
    @param[in] a_extrap_cdr_fluxes Extrapolated fluxes from the gas side. 
  */
  virtual Vector<Real> compute_cdr_dielectric_fluxes(const Real&         a_time,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const RealVect&     a_E,
						     const Vector<Real>& a_cdr_densities,
						     const Vector<Real>& a_cdr_velocities,
						     const Vector<Real>& a_rte_fluxes,
						     const Vector<Real>& a_extrap_cdr_fluxes) const = 0;

  /*!
    @brief Compute the isotropic source terms for the RTE
  */
  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const = 0;

  /*!
    @brief Set the initial surface charge
    @param[in] a_time Time
    @param[in] a_pos Position
  */
  virtual Real initial_sigma(const Real      a_time,
			     const RealVect& a_pos) const {};

  
  // ====================================================================================================
  // BELOW HERE IS THE OLD STUFF
  // ====================================================================================================

  virtual Vector<Real> compute_diffusion_coefficients(const RealVect& a_E) const = 0;
  
  
  /*!
    @brief Compute the velocity for an ion
    @param[in] a_E Electric field strength
  */
  virtual Vector<RealVect> compute_velocities(const RealVect& a_E) const = 0;



  /*!
    @brief Compute the source terms for each charge-carrier species
    @param[in] a_densities Various ion and electron densities
    @param[in] a_intensities Various photon intensities
    @param[in] a_E Electric field strength
    @details The output vector should contain the right-hand sides of the various photon solvers in the order as m_ions. 
    The input vectors follow this order.
  */
  virtual Vector<Real> compute_source_terms(const Vector<Real>& a_densities, 
					    const Vector<Real>& a_intensities,
					    const RealVect&     a_E) const = 0;

  /*!
    @brief Compute the source terms (right-hand sides) for RTE equations
    @param[in] a_densities Ion densities
    @param[in] a_E Electric field strength
    @details The output vector should contain the right-hand sides of the various photon solvers in the order as m_photons.
    The input vectors follow this order. 
  */
  virtual Vector<Real> compute_rte_source_terms(const Vector<Real>& a_densities, const RealVect& a_E) const = 0;



  /*!
    @brief Compute ion fluxes through conductor-gas interfaces
    @param[in] a_extrapolated_fluxes Extrapolated fluxes from the gas-side. Useful for inflow/outflow BCs
    @param[in] a_ion_densities Ion densities
    @param[in] a_ion_velocities Ion normal velocities
    @param[in] a_photon_fluxes Photon fluxes
    @param[in] a_E Electric field
    @param[in] a_pos Position
    @param[in] a_normal Normal vector (points into the gas phase)
    @param[in] a_time Time
  */
  virtual Vector<Real> compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						const Vector<Real>& a_ion_densities,
						const Vector<Real>& a_ion_velocities,
						const Vector<Real>& a_photon_fluxes,
						const RealVect&     a_E,
						const RealVect&     a_pos,
						const RealVect&     a_normal,
						const Real&         a_time) const = 0;
  /*!
    @brief Compute ion fluxes through dielectric-gas interfaces
    @param[in] a_extrapolated_fluxes Extrapolated fluxes from the gas-side. Useful for inflow/outflow BCs
    @param[in] a_ion_densities Ion densities
    @param[in] a_ion_velocities Ion normal velocities
    @param[in] a_photon_fluxes Photon fluxes
    @param[in] a_E Electric field
    @param[in] a_pos Position
    @param[in] a_normal Normal vector (points into the gas phase)
    @param[in] a_time Time
  */
  virtual Vector<Real> compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						 const Vector<Real>& a_ion_densities,
						 const Vector<Real>& a_ion_velocities,
						 const Vector<Real>& a_photon_fluxes,
						 const RealVect&     a_E,
						 const RealVect&     a_pos,
						 const RealVect&     a_normal,
						 const Real&         a_time) const = 0;

  /*!
    @brief Set the initial surface charge
    @param[in] a_pos Position
  */
  virtual Real initial_sigma(const RealVect& a_pos) const = 0;

  /*!
    @brief Get all ions
  */
  const Vector<RefCountedPtr<species> >& get_species() const { 
    return m_species; 
  }

  /*!
    @brief Get all photons
  */
  const Vector<RefCountedPtr<photon_group> >& get_photons() const {
    return m_photons;
  }

  /*!
    @brief Return number of ion equations
  */
  const int get_num_species() const{
    return m_species.size();
  }

  /*!
    @brief Return number of RTE equations
  */
  const int get_num_photons() const {
    return m_photons.size();
  }
  
protected:

  /*!
    @brief List of species
  */
  Vector<RefCountedPtr<species> > m_species;

  /*!
    @brief List of optical transitions between species
  */
  Vector<RefCountedPtr<photon_group> > m_photons;

  /*!
    @brief Number of species
  */
  int m_num_species;

  /*!
    @brief Number of RTE groups
  */
  int m_num_photons;
};

#endif
