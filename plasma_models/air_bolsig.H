/*!
  @file   air_bolsig.H
  @brief  Declaration of a BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _AIR_BOLSIG_
#define _AIR_BOLSIG_

#include "plasma_kinetics.H"
#include "lookup_table.H"
#include "perlin_if.H"

/*!
  @brief BOLSIG-based kinetic model for air. Transport coefficients for electrons are taken by calls to BOLSIG+. 
  @details 

  Class options
  -------------

  By default, the gas is electrically neutral with initial data supplied through an equal mixture of electrons and positive ions. The class supports uniform distributions, seeds, and random distributions of the initial ionization. Some kinetic things that are not accessible through BOLSIG can be controlled directly here. The following options are available for this class:

      air_bolsig.bolsig_path                   = ./plasma_models/bolsigplus032016-linux
      air_bolsig.lxcat_path                    = ./plasma_models/bolsigplus032016-linux
      air_bolsig.lxcat_file                    = LXCat-June2013.txt    # Where to find database file (LXCat)
      air_bolsig.grid_points                   = 1000                  # Number of grid points for BOLSIG call
      air_bolsig.min_townsend                  = 0.0                   # Transport coefficient computation limit (low)
      air_bolsig.max_townsend                  = 2000                  # Transport coefficient computation limit (low)
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
      air_bolsig.photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
      air_bolsig.photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
      air_bolsig.photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
      air_bolsig.photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
      air_bolsig.photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
      air_bolsig.photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model
      
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

  /*!
    @brief Compute the velocity for an ion
  */
  virtual Vector<RealVect> compute_velocities(const RealVect& a_E) const;

  /*!
    @brief Compute the source terms for each charge-carrier species
  */
  virtual Vector<Real> compute_source_terms(const Vector<Real>& a_densities, 
					    const Vector<Real>& a_intensities,
					    const RealVect&     a_E) const;

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
						const Vector<Real>& a_photon_fluxes,
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
						 const Vector<Real>& a_photon_fluxes,
						 const RealVect&     a_E,
						 const RealVect&     a_pos,
						 const RealVect&     a_normal,
						 const Real&         a_time) const;

  /*!
    @brief Set the initial surface charge
  */
  virtual Real initial_sigma(const RealVect& a_pos) const; 

protected:

  /*!
    @brief Name of temporary script file
  */
  static std::string s_script_file;

  /*!
    @brief Name of temporary data file
  */
  static std::string s_data_file;

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
    @brief String that we use to identify electron excitation rate e + N2 -> N2_ex(12.25eV) + e
  */
  static std::string s_bolsig_kex23;

  /*!
    @brief String that we use to identify electron excitation rate e + N2 -> N2_ex(13eV) + e
  */
  static std::string s_bolsig_kex24;

  /*!
    @brief Path to BOLSIG
  */
  std::string m_bolsig_path;

  /*!
    @brief Path to LXCAT data file
  */
  std::string m_lxcat_path;

  /*!
    @brief LXCAT file name
  */
  std::string m_lxcat_file;

  /*!
    @brief Name of temporary lxcat file
  */
  std::string m_local_lxcat_file;
  
  /*!
    @brief Name of temporary script file
  */
  std::string m_local_script_file;

  /*!
    @brief Name of temporary data file
  */
  std::string m_local_data_file;

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
    @brief Electron excitation rate e + N2 -> N2^*(12.25eV) + e
  */
  lookup_table m_kex23;

  /*!
    @brief Electron excitation rate e + N2 -> N2^*(13eV) + e
  */
  lookup_table m_kex24;

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
    @brief Number of grid points to use for BOLSIG computation
  */
  int m_grid_points;

  /*!
    @brief Number of octaves to use for noise function
  */
  int m_noise_octaves;

  /*!
    @brief Minimum electric field to compute transport coefficients from
  */
  Real m_min_townsend;
  
  /*!
    @brief Maximum electric field to compute transport coefficients from
  */
  Real m_max_townsend;

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
    @brief Work function for dielectric
  */
  Real m_dielectric_work;

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
    @brief Build the BOLSIG input script
  */
  virtual void build_bolsig_script();

  /*!
    @brief Let BOLSIG build the input script
  */
  virtual void call_bolsig();

  /*!
    @brief Extract necessary data from BOLSIG output file
  */
  virtual void extract_bolsig_data();

  /*!
    @brief Delete the BOLSIG output file
  */
  virtual void delete_bolsig_data();

  /*!
    @brief Delete the BOLSIG script
  */
  virtual void delete_bolsig_script();

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
    virtual const Real initial_data(const RealVect a_pos, const Real a_time) const;

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
    virtual const Real initial_data(const RealVect a_pos, const Real a_time) const;

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
    virtual const Real initial_data(const RealVect a_pos, const Real a_time) const;

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
    const Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    const Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    const Real get_pO2() const {
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
    const Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    const Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    const Real get_pO2() const {
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
    const Real get_lambda() const {
      return m_lambda;
    }

    /*!
      @brief Get the A-coefficient
    */
    const Real get_A() const {
      return m_A;
    }

    /*!
      @brief Get the partial oxygen pressure
    */
    const Real get_pO2() const {
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

};
#endif
