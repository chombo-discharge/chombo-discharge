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
#include <LoHiSide.H>

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
  plasma_kinetics(){
    m_eed_solve     = false;
    m_eed_index = 0;
  }

  /*!
    @brief Destructor
  */
  virtual ~plasma_kinetics(){}

  /*!
    @brief Routine intended for advancing a full reaction network simultaneous. 
    @details This is the routine that sometimes the reaction network in the FHD description. Things
             that come in/out of this routine are DENSITIES, but you should translate this to number
             of particles/photons per cell. The particle densities are probably not integers, i.e. cells may
	     contain "half an electron". If you're doing nearest-grid point deposition of photons, the
	     photon densities x cell volume should be an integer, although stored as a floating point. 
    @param[out] a_particle_sources   Number of particles produced per unit volume per unit time
    @param[out] a_particle_sources   Number of photons produced per unit volume per unit time
    @param[in]  a_particle_densities Grid-based density for particle species
    @param[in]  a_photon_densities   Grid-based densities for photons
    @param[in]  a_time               Current time
    @param[in]  a_dt                 Advanced time
    @param[in]  a_pos                Position in space
    @param[in]  a_E                  Electric field
  */
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
						const Real              a_kappa,
						const Real              a_dx,
						const RealVect         a_pos,
						const RealVect         a_E,
						const RealVect         a_gradE,
						const Vector<Real>     a_cdr_densities,
						const Vector<Real>     a_rte_densities,
						const Vector<RealVect> a_grad_cdr) const = 0;

  /*!
    @brief Compute the isotropic source terms for the RTE equations
  */
  virtual Vector<Real> compute_rte_source_terms(const Real         a_time,
						const Real         a_kappa,
						const Real         a_dx,
						const RealVect     a_pos,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities) const = 0;

  /*!
    @brief Compute velocities for the CDR equations
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_E Electric field
    @param[in] a_cdr_densities CDR densities
  */
  virtual Vector<RealVect> compute_cdr_velocities(const Real         a_time,
						  const RealVect     a_pos,
						  const RealVect     a_E,
						  const Vector<Real> a_cdr_densities) const = 0;

    /*!
    @brief Compute the cdr diffusion coefficients
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_E Electric field
    @param[in] a_cdr_densities CDR densities
  */
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real         a_time,
							  const RealVect     a_pos,
							  const RealVect     a_E,
							  const Vector<Real> a_cdr_densities) const = 0;


  /*!
    @brief Compute cdr fluxes through electrode-gas interfaces
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_normal Normal vector, points into the gas phase. 
    @param[in] a_E Electric field
    @param[in] a_cdr_velocities CDR velocities. Normal component only. 
    @param[in] a_cdr_densities CDR densities. 
    @param[in] a_cdr_gradients Normal gradients of cdr densities
    @param[in] a_rte_fluxes RTE fluxes (normal component only)
    @param[in] a_extrap_cdr_fluxes Extrapolated fluxes from the gas side. 
  */
  virtual Vector<Real> compute_cdr_electrode_fluxes(const Real         a_time,
						    const RealVect     a_pos,
						    const RealVect     a_normal,
						    const RealVect     a_E,
						    const Vector<Real> a_cdr_densities,
						    const Vector<Real> a_cdr_velocities,
						    const Vector<Real> a_cdr_gradients,
						    const Vector<Real> a_rte_fluxes,
						    const Vector<Real> a_extrap_cdr_fluxes) const = 0;

  /*!
    @brief Compute cdr fluxes through dielectric-gas interfaces
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_normal Normal vector, points into the gas phase. 
    @param[in] a_E Electric field
    @param[in] a_cdr_densities CDR densities. 
    @param[in] a_cdr_velocities CDR velocities. Normal component only. 
    @param[in] a_cdr_gradients Normal gradients of cdr densities
    @param[in] a_rte_fluxes RTE fluxes (normal component only)
    @param[in] a_extrap_cdr_fluxes Extrapolated fluxes from the gas side. 
  */
  virtual Vector<Real> compute_cdr_dielectric_fluxes(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_normal,
						     const RealVect     a_E,
						     const Vector<Real> a_cdr_densities,
						     const Vector<Real> a_cdr_velocities,
						     const Vector<Real> a_cdr_gradients,
						     const Vector<Real> a_rte_fluxes,
						     const Vector<Real> a_extrap_cdr_fluxes) const = 0;

  /*!
    @brief Compute cdr fluxes through domain sides
    @param[in] a_time Time
    @param[in] a_pos Position
    @param[in] a_dir Direction (0 = x, 1=y etc)
    @param[in] a_side Side (low or high side)
    @param[in] a_E Electric field
    @param[in] a_cdr_densities CDR densities. 
    @param[in] a_cdr_velocities CDR velocities. Normal component only. 
    @param[in] a_cdr_gradients Normal gradients of cdr densities
    @param[in] a_rte_fluxes RTE fluxes (normal component only)
    @param[in] a_extrap_cdr_fluxes Extrapolated fluxes from the gas side. 
  */
  virtual Vector<Real> compute_cdr_domain_fluxes(const Real           a_time,
						 const RealVect       a_pos,
						 const int            a_dir,
						 const Side::LoHiSide a_side,
						 const RealVect       a_E,
						 const Vector<Real>   a_cdr_densities,
						 const Vector<Real>   a_cdr_velocities,
						 const Vector<Real>   a_cdr_gradients,
						 const Vector<Real>   a_rte_fluxes,
						 const Vector<Real>   a_extrap_cdr_fluxes) const = 0;


  /*!
    @brief Set the initial surface charge
    @param[in] a_time Time
    @param[in] a_pos Position
  */
  virtual Real initial_sigma(const Real      a_time,
			     const RealVect a_pos) const = 0;

  /*!
    @brief Initialize electron energy density
  */
  virtual Real init_eed(const RealVect a_pos, const Real a_time, const RealVect a_E){
    return 0.0;
  }

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
  int get_num_species() const{
    return m_species.size();
  }

  /*!
    @brief Return number of RTE equations
  */
  int get_num_photons() const {
    return m_photons.size();
  }

  /*!
    @brief Check if solver solves for eed
  */
  bool solve_eed(){
    return m_eed_solve;
  }

  /*!
    @brief Get eed solver index
  */
  int get_eed_index(){
    return m_eed_index;
  }

  void set_dt(const Real a_dt){
    m_dt = a_dt;
  }

  void set_finest_dx(const Real a_finest_dx){
    m_finest_dx = a_finest_dx;
  }

  
protected:
  bool m_eed_solve;     // Doing EED solve or not?
  
  int m_eed_index;      // EED solver index
  
  Real m_dt;            // FHD kinetics need access to time steps
  Real m_finest_dx;     // Finest dx

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
