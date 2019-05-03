/*!
  @file   morrow_lowke.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _MORROW_LOWKE_
#define _MORROW_LOWKE_

#include "plasma_kinetics.H"
#include "perlin_if.H"

/*!
  @brief Class that implements the Morrow-Lowke discharge model for air.
  @details This class uses the Morrow-Lowke kinetic scheme for air and the three-group Bourdon model for the RTE equations. 
  
  Class options
  -------------

  By default, the gas is electrically neutral with initial data supplied through an equal mixture of electrons and positive ions. The class supports uniform distributions, seeds, and random distributions of the initial ionization. 

      morrow_lowke.gas_temperature               = 300                   # Gas temperature (K)
      morrow_lowke.gas_N2_frac                   = 0.8                   # Mixing fraction of nitrogen
      morrow_lowke.gas_O2_frac                   = 0.2                   # Mixing fraction of oxygen
      morrow_lowke.gas_pressure                  = 1.0                   # Gas pressure
      morrow_lowke.gas_quenching_pressure        = 0.03947               # Quenching pressure for photo-emission
      morrow_lowke.positive_species_mobility     = 2.E-4                 # Positive species mobility
      morrow_lowke.negative_species_mobility     = 2.E-4                 # Negative species mobility
      morrow_lowke.excitation_efficiency         = 0.6                   # Impact excitation efficiency
      morrow_lowke.photoionization_efficiency    = 0.1                   # Photo-ionization efficiency
      morrow_lowke.electrode_townsend2           = 1.E-6                 # Townsend coefficient on electrodes
      morrow_lowke.electrode_quantum_efficiency  = 1.E-6                 # Quantum efficiency on electrodes
      morrow_lowke.dielectric_townsend2          = 1.E-6                 # Townsend coefficient on dielectrics
      morrow_lowke.dielectric_quantum_efficiency = 1.E-6                 # Quantum efficiency on dielectrics
      morrow_lowke.dielectric_work_function      = 1.E-6                 # Dielectric work function (ev)
      morrow_lowke.uniform_density               = 1.E10                 # Initial background ionization
      morrow_lowke.seed_density                  = 0.E16                 # Initial seed ionization
      morrow_lowke.seed_radius                   = 1E-3                  # Initial seed radius
      morrow_lowke.seed_position                 = 0. 0. 0.              # Initial seed position
      morrow_lowke.noise_amplitude               = 0.0                   # Initial noise amplitude
      morrow_lowke.noise_octaves                 = 1                     # Initial noise octaves
      morrow_lowke.noise_persistence             = 0.5                   # Reduction factor between each noise octave
      morrow_lowke.noise_frequency               = 1.0 1.0 1.0           # Spatial noise frequency
      morrow_lowke.photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
      morrow_lowke.photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model

*/
class morrow_lowke : public plasma_kinetics {
public:

  /*!
    @brief Constructor
  */
  morrow_lowke();

  /*!
    @brief Destructor
  */
  virtual ~morrow_lowke();

  // ====================================================================================================
// NEW FUNCTION SIGNATURES

  /*!
    @brief Compute diffusion coefficients
  */
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const {
    return this->compute_diffusion_coefficients(a_E);
  }

  /*!
    @brief Compute velocities
  */
  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const {
    return this->compute_velocities(a_E);
  }

  /*!
    @brief Compute source terms
  */
  virtual Vector<Real> compute_cdr_source_terms(const Real              a_time,
						const Real              a_kappa,
						const Real              a_dx,
						const RealVect&         a_pos,
						const RealVect&         a_E,
						const RealVect&         a_gradE,
						const Vector<Real>&     a_cdr_densities,
						const Vector<Real>&     a_rte_densities,
						const Vector<RealVect>& a_grad_cdr) const {
    return this->compute_source_terms(a_cdr_densities, a_rte_densities, a_E);
  }

  /*!
    @brief Compute electrode fluxes
  */
  virtual Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
						    const RealVect&     a_pos,
						    const RealVect&     a_normal,
						    const RealVect&     a_E,
						    const Vector<Real>& a_cdr_densities,
						    const Vector<Real>& a_cdr_velocities,
						    const Vector<Real>& a_cdr_gradients,
						    const Vector<Real>& a_rte_fluxes,
						    const Vector<Real>& a_extrap_cdr_fluxes) const {
    return this->compute_conductor_fluxes(a_extrap_cdr_fluxes,
					  a_cdr_densities,
					  a_cdr_velocities,
					  a_rte_fluxes,
					  a_E,
					  a_pos,
					  a_normal,
					  a_time);
  }

  /*!
    @brief Compute dielectric flxues
  */
  virtual Vector<Real> compute_cdr_dielectric_fluxes(const Real&         a_time,
						     const RealVect&     a_pos,
						     const RealVect&     a_normal,
						     const RealVect&     a_E,
						     const Vector<Real>& a_cdr_densities,
						     const Vector<Real>& a_cdr_velocities,
						     const Vector<Real>& a_cdr_gradients,
						     const Vector<Real>& a_rte_fluxes,
						     const Vector<Real>& a_extrap_cdr_fluxes) const {
    return this->compute_dielectric_fluxes(a_extrap_cdr_fluxes,
					   a_cdr_densities,
					   a_cdr_velocities,
					   a_rte_fluxes,
					   a_E,
					   a_pos,
					   a_normal,
					   a_time);
  }

  /*!
    @brief Compute the isotropic source terms for the RTE
  */
  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const Real&         a_kappa,
						const Real&         a_dx,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const {
    return this->compute_rte_source_terms(a_cdr_densities, a_E);
  }

  /*!
    @brief Compute the velocity for an ion
    @param[in] a_E Electric field strength
  */
  virtual Vector<RealVect> compute_velocities(const RealVect& a_E) const;

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
					    const RealVect&     a_E) const;

  /*!
    @brief Compute the source terms (right-hand sides) for RTE equations
    @param[in] a_densities species densities
    @param[in] a_E Electric field strength
    @details The output vector should contain the right-hand sides of the various photon solvers in the order as m_photons.
    The input vectors follow this order. 
  */
  virtual Vector<Real> compute_rte_source_terms(const Vector<Real>& a_densities, const RealVect& a_E) const;

  /*!
    @brief Compute the ion diffusion coefficient
    @param[in] a_E Electric field strength
  */
  virtual Vector<Real> compute_diffusion_coefficients(const RealVect& a_E) const;

  /*!
    @brief Compute ion fluxes through conductor-gas interfaces
    @param[in] a_extrapolated_fluxes Extrapolated fluxes from the gas-side. Useful for inflow/outflow BCs
    @param[in] a_ion_densities species densities
    @param[in] a_ion_velocities species normal velocities
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
						const Real&         a_time) const;

  /*!
    @brief Compute ion fluxes through dielectric-gas interfaces
    @param[in] a_extrapolated_fluxes Extrapolated fluxes from the gas-side. Useful for inflow/outflow BCs
    @param[in] a_ion_densities species densities
    @param[in] a_ion_velocities species normal velocities
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
						 const Real&         a_time) const;

  /*!
    @brief Set the initial surface charge
  */
  virtual Real initial_sigma(const Real a_time, const RealVect& a_pos) const;


  /*!
    @brief Implementation of electron for this class
  */
  class electron : public species {
  public:
    /*!
      @brief Constructor
    */
    electron();

    /*!
      @brief Destructor
    */
    virtual ~electron();

    /*!
      @brief Implementation of initial data
    */
    virtual Real initial_data(const RealVect a_pos, const Real a_time) const;

    /*!
      @brief Set the noise function
    */
    virtual void set_noise(RefCountedPtr<perlin_if> a_noisePtr);

    /*!
      @brief Noise function
    */
    RefCountedPtr<perlin_if> m_perlin;

    /*!
      @brief Uniform plasma density
    */
    Real m_uniform_density;
    
    /*!
      @brief Plasma seed density
    */
    Real m_seed_density;

    /*!
      @brief Plasma seed radius
    */
    Real m_seed_radius;

    /*!
      @brief Plasma seed radius
    */
    Real m_noise_density;

    /*!
      @brief Plasma seed position
    */
    RealVect m_seed_pos;
  };

  /*!
    @brief Implementation of electron for this class
  */
  class positive_species : public species {
  public:
    /*!
      @brief Constructor
    */
    positive_species();

    /*!
      @brief Destructor
    */
    virtual ~positive_species();

    /*!
      @brief Implementation of initial data
    */
    virtual Real initial_data(const RealVect a_pos, const Real a_time) const;

    /*!
      @brief Set the noise function
    */
    virtual void set_noise(RefCountedPtr<perlin_if> a_noisePtr);

    /*!
      @brief Noise function
    */
    RefCountedPtr<perlin_if> m_perlin;

    /*!
      @brief Uniform plasma density
    */
    Real m_uniform_density;
    
    /*!
      @brief Plasma seed density
    */
    Real m_seed_density;

    /*!
      @brief Plasma seed radius
    */
    Real m_seed_radius;

    /*!
      @brief Plasma seed radius
    */
    Real m_noise_density;

    /*!
      @brief Plasma seed position
    */
    RealVect m_seed_pos;
  };

  /*!
    @brief Implementation of the negative species for the Morrow-Lowke model
  */
  class negative_species : public species {
  public:
    /*!
      @brief Constructor
    */
    negative_species();

    /*!
      @brief Destructor
    */
    ~negative_species();

    /*!
      @brief Implementation of initial data
    */
    virtual Real initial_data(const RealVect a_pos, const Real a_time) const;

    /*!
      @brief Set the noise function
    */
    virtual void set_noise(RefCountedPtr<perlin_if> a_noisePtr);

    /*!
      @brief Noise function
    */
    RefCountedPtr<perlin_if> m_perlin;

    /*!
      @brief Uniform plasma density
    */
    Real m_uniform_density;
    
    /*!
      @brief Plasma seed density
    */
    Real m_seed_density;

    /*!
      @brief Plasma seed radius
    */
    Real m_seed_radius;

    /*!
      @brief Plasma seed radius
    */
    Real m_noise_density;

    /*!
      @brief Plasma seed position
    */
    RealVect m_seed_pos;
  };

  /*!
    @brief Implementation of a single photon transition
  */
  class photon_one : public photon_group {
  public:
    /*!
      @brief Constructor
    */
    photon_one();

    /*!
      @brief Destructor
    */
    ~photon_one();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real get_kappa(const RealVect a_pos) const;

    /*!
      @brief Get the lambda coefficient
    */
    Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    Real get_pO2() const {
      return m_pO2;
    }

  protected:
    /*!
      @brief Lambda-coefficient
    */
    Real m_lambda;

    /*!
      @brief A-coefficient
    */
    Real m_A;

    /*!
      @brief O2-fraction
    */
    Real m_pO2;
  };

  /*!
    @brief Implementation of a single photon transition
  */
  class photon_two : public photon_group {
  public:
    /*!
      @brief Constructor
    */
    photon_two();

    /*!
      @brief Destructor
    */
    ~photon_two();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real get_kappa(const RealVect a_pos) const;

    /*!
      @brief Get the lambda coefficient
    */
    Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    Real get_pO2() const {
      return m_pO2;
    }

  protected:
    /*!
      @brief Lambda-coefficient
    */
    Real m_lambda;

    /*!
      @brief A-coefficient
    */
    Real m_A;

    /*!
      @brief O2-fraction
    */
    Real m_pO2;
  };

  /*!
    @brief Implementation of a single photon transition
  */
  class photon_three : public photon_group {
  public:
    /*!
      @brief Constructor
    */
    photon_three();

    /*!
      @brief Destructor
    */
    ~photon_three();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real get_kappa(const RealVect a_pos) const;

    /*!
      @brief Get the lambda coefficient
    */
    Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    Real get_pO2() const {
      return m_pO2;
    }

  protected:
    /*!
      @brief Lambda-coefficient
    */
    Real m_lambda;

    /*!
      @brief A-coefficient
    */
    Real m_A;

    /*!
      @brief O2-fraction
    */
    Real m_pO2;
  };

  /*!
    @brief Perlin noise function
  */
  RefCountedPtr<perlin_if> m_perlin;

  /*!
    @brief electron index 
  */
  int m_nelec_idx;

  /*!
    @brief Positive species index
  */
  int m_nplus_idx;

  /*!
    @brief Negative species index
  */
  int m_nminu_idx;

  /*!
    @brief photon_one index
  */
  int m_photon1_idx;

  /*!
    @brief photon_two index
  */
  int m_photon2_idx;

  /*!
    @brief photon_three index
  */
  int m_photon3_idx;

  /*!
    @brief Number of octaves for initial data noise
  */
  int m_noise_octaves;

  /*!
    @brief Background ionization rate
  */
  Real m_background_rate;

  /*!
    @brief Second Townsend coefficient for conductors
  */
  Real m_townsend2_conductor;

  /*!
    @brief Second Townsend coefficient for dielectrics
  */
  Real m_townsend2_dielectric;

  /*!
    @brief Photo-emission yield on electrodes
  */
  Real m_electrode_yield;

  /*!
    @brief Photo-emission yield on dielectrics
  */
  Real m_dielectric_yield;

  /*!
    @brief Work function for the dielectric
  */
  Real m_dielectric_work;

  /*!
    @brief Real Uniform background density
  */
  Real m_uniform_denisty;

  /*!
    @brief Neutral number density
  */
  Real m_N;

  /*!
    @brief Oxygen fraction
  */
  Real m_fracO2;

  /*!
    @brief Nitrogen fraction
  */
  Real m_fracN2;

  /*!
    @brief Quenching pressure
  */
  Real m_pq;

  /*!
    @brief Pressure
  */
  Real m_p;

  /*!
    @brief Excitation efficiency
  */
  Real m_exc_eff;

  /*!
    @brief Photoionization efficiency
  */
  Real m_photo_eff;

  /*!
    @brief Gas temperature
  */
  Real m_temp;

  /*!
    @brief Noise amplitude for initial data
  */
  Real m_noise_amp;

  /*!
    @brief Reduction factor between each octave of noise
  */
  Real m_noise_persist;

    /*!
    @brief Perlin noise frequency for initial data
  */
  RealVect m_noise_freq;

  /*!
    @brief Compute the flux on the cathode
  */
  virtual Vector<Real> compute_cathode_flux(const Vector<Real>& a_extrapolated_fluxes,
					    const Vector<Real>& a_ion_densities,
					    const Vector<Real>& a_ion_velocities,
					    const Vector<Real>& a_photon_fluxes,
					    const RealVect&     a_E,
					    const RealVect&     a_pos,
					    const RealVect&     a_normal,
					    const Real&         a_time) const;

  /*!
    @brief Compute the flux on the anode
  */
  virtual Vector<Real> compute_anode_flux(const Vector<Real>& a_extrapolated_fluxes,
					  const Vector<Real>& a_ion_densities,
					  const Vector<Real>& a_ion_velocities,
					  const Vector<Real>& a_photon_fluxes,
					  const RealVect&     a_E,
					  const RealVect&     a_pos,
					  const RealVect&     a_normal,
					  const Real&         a_time) const;


  /*!
    @brief Compute \f$\alpha\f$-coefficient
  */
  virtual Real compute_alpha(const RealVect& a_E) const;

  /*!
    @brief Compute \f$\eta\f$ coefficient
  */
  virtual Real compute_eta(const RealVect& a_E) const;

  /*!
    @brief Compute \f$eta_2\f$-coefficient
  */
  virtual Real compute_eta2(const RealVect& a_E) const;

  /*!
    @brief Compute \f$eta_3\f$-coefficient
  */
  virtual Real compute_eta3(const RealVect& a_E) const;

  /*!
    @brief Compute beta-coefficient
  */
  virtual Real compute_beta(const RealVect& a_E) const;

  /*!
    @brief Compute electron diffusion coefficient
  */
  virtual Real compute_De(const RealVect& a_E) const;

  /*!
    @brief Compute electron velocity
  */
  virtual RealVect compute_ve(const RealVect& a_E) const;

  /*!
    @brief Compute positive species velocity
  */
  virtual RealVect compute_vp(const RealVect& a_E) const;

  /*!
    @brief Compute negative species velocity
  */
  virtual RealVect compute_vn(const RealVect& a_E) const;

};

#endif
