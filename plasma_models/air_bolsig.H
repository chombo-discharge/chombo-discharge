/*!
  @file   air_bolsig.H
  @brief  Declaration of a BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _AIR_BOLSIG_
#define _AIR_BOLSIG_

#include "plasma_kinetics.H"

/*!
  @brief BOLSIG-based kinetic model for air. Transport coefficients for electrons are taken by calls to BOLSIG+. 
  @details 

  Class options
  -------------

  By default, the gas is electrically neutral with initial data supplied through an equal mixture of electrons and positive ions. The class supports uniform distributions, seeds, and random distributions of the initial ionization. Some kinetic things that are not accessible through BOLSIG can be controlled directly here. The following options are available for this class:

      air_bolsig.grid_points                  = 1000           # Number of grid points for BOLSIG call
      air_bolsig.min_townsend                 = 0.0            # Transport coefficient computation limit (low)
      air_bolsig.max_townsend                 = 2000           # Transport coefficient computation limit (low)
      air_bolsig.gas_temperature              = 300            # Gas temperature (K)
      air_bolsig.gas_N2_frac                  = 0.8            # Mixing fraction of nitrogen
      air_bolsig.gas_O2_frac                  = 0.2            # Mixing fraction of oxygen
      air_bolsig.gas_pressure                 = 1.0            # Gas pressure
      air_bolsig.quenching_pressure           = 0.03947        # Quenching pressure for photo-emission
      air_bolsig.ionization_rate              = 1.E9           # Background ionization (events/s)
      air_bolsig.electron_recombination       = 5.E-14         # Electron recombination rate
      air_bolsig.electron_detachment          = 1.E-18         # Electron detachment rate
      air_bolsig.ion_recombination            = 2.07E-12       # Ion recombination rate
      air_bolsig.positive_species_mobility    = 2.E-4          # Positive species mobility
      air_bolsig.negative_species_mobility    = 2.E-4          # Negative species mobility
      air_bolsig.excitation_eff               = 0.6            # Impact excitation efficiency
      air_bolsig.photo_eff                    = 0.1            # Photo-ionization efficiency
      air_bolsig.electrode_townend2           = 1.E-6          # Townsend coefficient on electrodes
      air_bolsig.electrode_quantum_eff        = 1.E-6          # Quantum efficiency on electrodes
      air_bolsig.dielectric_townsend2         = 1.E-6          # Townsend coefficient on dielectrics
      air_bolsig.dielectric_quantum_eff       = 1.E-6          # Quantum efficiency on dielectrics
      air_bolsig.uniform_density              = 1.E10          # Initial background ionization
      air_bolsig.seed_density                 = 1.E16          # Initial seed ionization
      air_bolsig.seed_radius                  = 1E-3           # Initial seed radius
      air_bolsig.seed_position                = 0. 0. 0.       # Initial seed position
      air_bolsig.noise_amplitude              = 0.0            # Initial noise amplitude
      air_bolsig.noise_octaves                = 1              # Initial noise octaves
      air_bolsig.noise_persist                = 0.5            # Reduction factor between each noise octave
      air_bolsig.noise_frequency              = 1.0 1.0 1.0    # Spatial noise frequency
      
*/
class air_bolsig : public plasma_kinetics {
public:

  class electron;
  class positive_species;
  class negative_species;

  class photon_one;
  class photon_two;
  class photon_three;

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
    @brief Path to BOLSIG
  */
  static std::string s_bolsig_path;

  /*!
    @brief Path to LXCAT data file
  */
  static std::string s_lxcat_path;

  /*!
    @brief LXCAT file name
  */
  static std::string s_lxcat_file;

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
    @brief Minimum electric field to compute transport coefficients from
  */
  Real m_bolsig_min_townsend;
  
  /*!
    @brief Maximum electric field to compute transport coefficients from
  */
  Real m_bolsig_max_townsend;

  /*!
    @brief Number of grid points to use for BOLSIG computation
  */
  int m_bolsig_grid_points;

  /*!
    @brief Call BOLSIG in order to read in electron transport data
  */
  virtual void call_bolsig();

  /*!
    @brief Build the BOLSIG input script
  */
  virtual void build_script();

  /*!
    @brief Let BOLSIG build the input script
  */
  virtual void build_data();

  /*!
    @brief Extract necessary data from BOLSIG output file
  */
  virtual void extract_data();

  /*!
    @brief Delete the BOLSIG output file
  */
  virtual void delete_data();

  /*!
    @brief Delete the BOLSIG script
  */
  virtual void delete_script();

};
#endif
