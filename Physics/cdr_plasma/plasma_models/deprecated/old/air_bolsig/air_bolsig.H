/*!
  @file   air_bolsig.H
  @brief  Declaration of a BOLSIG based transport model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _AIR_BOLSIG_
#define _AIR_BOLSIG_

#include "plasma_kinetics.H"
#include "lookup_table.H"
#include "perlin_if.H"

/*!
  @class air_bolsig 
  @brief BOLSIG-based kinetic model for air. Transport coefficients for electrons are taken by calls to BOLSIG+. 
  @details 

  Class options
  -------------

  By default, the gas is electrically neutral with initial data supplied through an equal mixture of electrons and positive ions. The class supports uniform distributions, seeds, and random distributions of the initial ionization. Some kinetic things that are not accessible through BOLSIG can be controlled directly here. The following options are available for this class:

  air_bolsig.transport_file                = transport_data.txt    # Transport data file
  air_bolsig.gas_temperature               = 300                   # Gas temperature (K)
  air_bolsig.gas_N2_frac                   = 0.8                   # Mixing fraction of nitrogen
  air_bolsig.gas_O2_frac                   = 0.2                   # Mixing fraction of oxygen
  air_bolsig.gas_pressure                  = 1.0                   # Gas pressure
  air_bolsig.gas_quenching_pressure        = 0.03947               # Quenching pressure for photo-emission
  air_bolsig.ionization_rate               = 1.E9                  # Background ionization (events/s)
  air_bolsig.electron_diffusion            = true                  # Electron diffusion on. Options are true/false
  air_bolsig.electron_recombination        = 5.E-14                # Electron recombination rate
  air_bolsig.electron_detachment           = 1.E-18                # Electron detachment rate
  air_bolsig.ion_recombination             = 2.07E-12              # Ion recombination rate
  air_bolsig.positive_species_mobility     = 2.E-4                 # Positive species mobility
  air_bolsig.negative_species_mobility     = 2.E-4                 # Negative species mobility
  air_bolsig.excitation_efficiency         = 0.6                   # Impact excitation efficiency
  air_bolsig.photoionization_efficiency    = 0.1                   # Photo-ionization efficiency
  air_bolsig.electrode_townsend2           = 1.E-6                 # Townsend coefficient on electrodes
  air_bolsig.electrode_quantum_efficiency  = 1.E-6                 # Quantum efficiency on electrodes
  air_bolsig.dielectric_townsend2          = 1.E-6                 # Townsend coefficient on dielectrics
  air_bolsig.dielectric_quantum_efficiency = 1.E-6                 # Quantum efficiency on dielectrics
  air_bolsig.dielectric_work_function      = 2.0                   # Dielectric work function (ev)
  air_bolsig.uniform_density               = 1.E10                 # Initial background ionization
  air_bolsig.seed_density                  = 0.E16                 # Initial seed ionization
  air_bolsig.seed_radius                   = 1E-3                  # Initial seed radius
  air_bolsig.seed_position                 = 0. 0. 0.              # Initial seed position
  air_bolsig.noise_amplitude               = 0.0                   # Initial noise amplitude
  air_bolsig.noise_octaves                 = 1                     # Initial noise octaves
  air_bolsig.noise_persistence             = 0.5                   # Reduction factor between each noise octave
  air_bolsig.noise_frequency               = 1.0 1.0 1.0           # Spatial noise frequency
  air_bolsig.Photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
  air_bolsig.Photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
  air_bolsig.Photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
  air_bolsig.Photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
  air_bolsig.Photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
  air_bolsig.Photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model
      
*/
class air_bolsig : public plasma_kinetics {
public:
  /*!
    @brief Constructor
  */
  air_bolsig();

  /*!
    @brief Destructor
  */
  virtual ~air_bolsig();

  // ====================================================================================================
  // NEW FUNCTION SIGNATURES

  /*!
    @brief Compute cdr diffusion coefficients
  */
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const {
    return this->compute_diffusion_coefficients(a_E);
  }

  /*!
    @brief Compute cdr velocities
  */
  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const {
    return this->compute_velocities(a_E);
  }

  /*!
    @brief Compute cdr source terms
  */
  virtual Vector<Real> compute_cdr_source_terms(const Real              a_time,
						const RealVect&         a_pos,
						const RealVect&         a_E,
						const RealVect&         a_gradE,
						const Vector<Real>&     a_cdr_densities,
						const Vector<Real>&     a_rte_densities,
						const Vector<RealVect>& a_grad_cdr) const;

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
    @brief Compute dielectric fluxes
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
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const {
    return this->compute_rte_source_terms(a_cdr_densities, a_E);
  }


  // NEW FUNCTION SIGNATURES END
  // ====================================================================================================

  /*!
    @brief Compute the velocity for an ion
  */
  virtual Vector<RealVect> compute_velocities(const RealVect& a_E) const;

  /*!
    @brief Compute the source terms (right-hand sides) for RTE equations
  */
  virtual Vector<Real> compute_rte_source_terms(const Vector<Real>& a_densities, const RealVect& a_E) const;

  /*!
    @brief Compute the ion diffusion coefficient
  */
  virtual Vector<Real> compute_diffusion_coefficients(const RealVect& a_E) const;

  /*!
    @brief Compute ion fluxes through conductor-gas interfaces
  */
  virtual Vector<Real> compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						const Vector<Real>& a_ion_densities,
						const Vector<Real>& a_ion_velocities,
						const Vector<Real>& a_Photon_fluxes,
						const RealVect&     a_E,
						const RealVect&     a_pos,
						const RealVect&     a_normal,
						const Real&         a_time) const;

  /*!
    @brief Compute ion fluxes through dielectric-gas interfaces
  */
  virtual Vector<Real> compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
						 const Vector<Real>& a_ion_densities,
						 const Vector<Real>& a_ion_velocities,
						 const Vector<Real>& a_Photon_fluxes,
						 const RealVect&     a_E,
						 const RealVect&     a_pos,
						 const RealVect&     a_normal,
						 const Real&         a_time) const;

  /*!
    @brief Set the initial surface charge
  */
  virtual Real initial_sigma(const Real a_time, const RealVect& a_pos) const; 

protected:

  /*!
    @brief String that we use to identify ionization coefficient in bolsig output
  */
  static std::string s_bolsig_alpha;

  /*!
    @brief String that we use to identify attachment coefficient in bolsig output
  */
  static std::string s_bolsig_eta;

  /*!
    @brief String that we use to identify electron mobility in bolsig output
  */
  static std::string s_bolsig_mob;

  /*!
    @brief String that we use to identify electron diffusion in bolsig output
  */
  static std::string s_bolsig_diff;

  /*!
    @brief Safety in case BOLSIG file was written with fitted functions
  */
  static std::string s_skip_fit;

  /*!
    @brief Name of temporary data file
  */
  std::string m_transport_file;

  /*!
    @brief Noise function
  */
  RefCountedPtr<perlin_if> m_perlin;
  
  /*!
    @brief Electron mobility
  */
  lookup_table m_electron_mobility;

  /*!
    @brief Electron mobility
  */
  lookup_table m_electron_diffusion;

  /*!
    @brief Impact ionization coefficient
  */
  lookup_table m_alpha;

  /*!
    @brief Impact attachment coefficient
  */
  lookup_table m_eta;

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
    @brief Photon_one index
  */
  int m_Photon1_idx;

  /*!
    @brief Photon_two index
  */
  int m_Photon2_idx;

  /*!
    @brief Photon_three index
  */
  int m_Photon3_idx;

  /*!
    @brief Number of octaves to use for noise function
  */
  int m_noise_octaves;

  /*!
    @brief 
  */
  Real m_gas_temp;

  /*!
    @brief O2 fraction
  */
  Real m_frac_O2;

  /*!
    @brief N2 fraction
  */
  Real m_frac_N2;

  /*!
    @brief Pressure
  */
  Real m_p;

  /*!
    @brief Quenching pressure
  */
  Real m_pq;

  /*!
    @brief Neutral number density
  */
  Real m_N;

  /*!
    @brief Background ionization rate
  */
  Real m_background_rate;

  /*!
    @brief Second Townsend coefficient for conductors
  */
  Real m_townsend2_electrode;

  /*!
    @brief Second Townsend coefficient for dielectrics
  */
  Real m_townsend2_dielectric;

  /*!
    @brief Quantum efficiency on electrodes
  */
  Real m_electrode_quantum_efficiency;

  /*!
    @brief Quantum efficiency on dielectric surfaces
  */
  Real m_dielectric_quantum_efficiency;

  /*!
    @brief Uniform background density
  */
  Real m_uniform_density;

  /*!
    @brief Excitation efficiency
  */
  Real m_excitation_efficiency;

  /*!
    @brief Photoionization efficiency
  */
  Real m_photoionization_efficiency;

  /*!
    @brief Ion mobilities
  */
  Real m_negative_species_mobility;

  /*!
    @brief Ion mobilities
  */
  Real m_positive_species_mobility;

  /*!
    @brief Electron-ion recombination rate
  */
  Real m_electron_recombination;

  /*!
    @brief ion-ion recombination rate
  */
  Real m_ion_recombination;

  /*!
    @brief Electron detachment coefficient
  */
  Real m_electron_detachment;

  /*!
    @brief Noise amplitude for initial data
  */
  Real m_noise_amplitude;

  /*!
    @brief Reduction factor between each octave of noise
  */
  Real m_noise_persistence;

  /*!
    @brief Perlin noise frequency for initial data
  */
  RealVect m_noise_frequency;

  /*!
    @brief Call BOLSIG in order to read in electron transport data
  */
  virtual void compute_transport_coefficients();

  /*!
    @brief Compute the flux on the cathode
  */
  virtual Vector<Real> compute_cathode_flux(const Vector<Real>& a_extrapolated_fluxes,
					    const Vector<Real>& a_ion_densities,
					    const Vector<Real>& a_ion_velocities,
					    const Vector<Real>& a_Photon_fluxes,
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
					  const Vector<Real>& a_Photon_fluxes,
					  const RealVect&     a_E,
					  const RealVect&     a_pos,
					  const RealVect&     a_normal,
					  const Real&         a_time) const;


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
    virtual Real initialData(const RealVect a_pos, const Real a_time) const;

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
    virtual Real initialData(const RealVect a_pos, const Real a_time) const;

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
    virtual Real initialData(const RealVect a_pos, const Real a_time) const;

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
    @brief Implementation of a single Photon transition
  */
  class Photon_one : public Photon_group {
  public:
    /*!
      @brief Constructor
    */
    Photon_one();

    /*!
      @brief Destructor
    */
    ~Photon_one();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real getKappa(const RealVect a_pos) const;

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
    @brief Implementation of a single Photon transition
  */
  class Photon_two : public Photon_group {
  public:
    /*!
      @brief Constructor
    */
    Photon_two();

    /*!
      @brief Destructor
    */
    ~Photon_two();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real getKappa(const RealVect a_pos) const;

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
    @brief Implementation of a single Photon transition
  */
  class Photon_three : public Photon_group {
  public:
    /*!
      @brief Constructor
    */
    Photon_three();

    /*!
      @brief Destructor
    */
    ~Photon_three();

    /*!
      @brief Get the absorption coefficient
    */
    virtual Real getKappa(const RealVect a_pos) const;

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

#include "CD_NamespaceFooter.H"
#endif
