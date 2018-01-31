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

  /*!
    @brief Compute the velocity for an ion
    @param[in] a_E Electric field strength
  */
  virtual Vector<RealVect> compute_velocities(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_velocities");
  }

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
					    const RealVect&     a_E) const{
    MayDay::Abort("morrow_lowke::compute_source_terms");
  }

  /*!
    @brief Compute the source terms (right-hand sides) for RTE equations
    @param[in] a_densities species densities
    @param[in] a_E Electric field strength
    @details The output vector should contain the right-hand sides of the various photon solvers in the order as m_photons.
    The input vectors follow this order. 
  */
  virtual Vector<Real> compute_rte_source_terms(const Vector<Real>& a_densities, const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_rte_source_terms");
  }

  /*!
    @brief Compute the ion diffusion coefficient
    @param[in] a_E Electric field strength
  */
  virtual Vector<Real> compute_diffusion_coefficients(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_diffusion_coefficients");
  }

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
						const Real&         a_time) const{
    MayDay::Abort("morrow_lowke::compute_conductor_fluxes");
  }
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
						 const Real&         a_time) const{
    MayDay::Abort("morrow_lowke::compute_dielectric_fluxes");
  }

  /*!
    @brief Set the initial surface charge
    @param[in] a_pos Position
  */
  virtual Real initial_sigma(const RealVect& a_pos) const;

protected:

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
					    const Real&         a_time) const{
    MayDay::Abort("morrow_lowke::compute_cathode_flux");
  }

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
					  const Real&         a_time) const{
    MayDay::Abort("morrow_lowke::compute_anode_flux");
  }


  /*!
    @brief Compute \f$\alpha\f$-coefficient
  */
  virtual Real compute_alpha(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_alpha");
  }

  /*!
    @brief Compute \f$\eta\f$ coefficient
  */
  virtual Real compute_eta(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_eta");
  }

  /*!
    @brief Compute \f$eta_2\f$-coefficient
  */
  virtual Real compute_eta2(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_eta2");
  }

  /*!
    @brief Compute \f$eta_3\f$-coefficient
  */
  virtual Real compute_eta3(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_eta3");
  }

  /*!
    @brief Compute beta-coefficient
  */
  virtual Real compute_beta(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_beta");
  }

  /*!
    @brief Compute electron diffusion coefficient
  */
  virtual Real compute_De(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_De");
  }

  /*!
    @brief Compute electron velocity
  */
  virtual RealVect compute_Ve(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_Ve");
  }

  /*!
    @brief Compute positive ion velocity
  */
  virtual RealVect compute_Vp(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_Vp");
  }

  /*!
    @brief Compute negative ion velocity
  */
  virtual RealVect compute_Vn(const RealVect& a_E) const{
    MayDay::Abort("morrow_lowke::compute_Vn");
  }

};

#endif
