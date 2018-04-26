/*!
  @file   air7.H
  @brief  Declaration of a 3-species BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR7_
#define _AIR7_

#include "plasma_kinetics.H"
#include "perlin_if.H"

/*!
  @brief 7-species BOLSIG-based transport model for air at 1atm and 300K. 
  temperature in the 

  @details This model includes variations due to the electron temperature in the LFA

  Class options
  -------------

      air7.electron_diffusion            = true                  # Turn off electron diffusion (debugging only)
      air7.excitation_efficiency         = 0.6                   # Impact excitation efficiency
      air7.thermal_outflow               = false                 # Thermal outflow on/off
      air7.diffusive_outflow             = false                 # Diffusive outflow on/off
      air7.photoionization_efficiency    = 0.1                   # Photo-ionization efficiency
      air7.electrode_townsend2           = 1.E-2                 # Townsend coefficient on electrodes
      air7.electrode_quantum_efficiency  = 1.E-2                 # Quantum efficiency on electrodes
      air7.dielectric_townsend2          = 1.E-1                 # Townsend coefficient on dielectrics
      air7.dielectric_quantum_efficiency = 1.E-1                 # Quantum efficiency on dielectrics
      air7.initial_ionization            = 1.E10                 # Initial background ionization
      air7.quenching_pressure            = 0.03947               # Quenching pressure for photo-emission
      air7.photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
      air7.photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
      air7.photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
      air7.photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
      air7.photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
      air7.photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model
      air7.simple_cathode_bc             = false                 # Use simpiflied BC on electrode cathodes

*/
class air7 : public plasma_kinetics {
public:

  class electron;
  class N2plus;
  class O2plus;
  class N4plus;
  class O4plus;
  class O2plusN2;
  class O2minus;

  class photon_one;
  class photon_two;
  class photon_three;

  /*!
    @brief Get gas temperature, pressure, density, and molecular fractions
  */
  static void get_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac);

  /*!
    @brief Constructor
  */
  air7();

  /*!
    @brief Destructor
  */
  virtual ~air7();

  /*!
    @brief Compute diffusion coefficients
  */
  virtual Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
							  const RealVect&     a_pos,
							  const RealVect&     a_E,
							  const Vector<Real>& a_cdr_densities) const;

  /*!
    @brief Compute drift velocities
  */
  virtual Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const;

  
  /*!
    @brief Compute source terms
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
						    const Vector<Real>& a_extrap_cdr_fluxes) const;

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
						     const Vector<Real>& a_extrap_cdr_fluxes) const;
  /*!
    @brief Compute the isotropic source terms for the RTEs
  */
  virtual Vector<Real> compute_rte_source_terms(const Real&         a_time,
						const RealVect&     a_pos,
						const RealVect&     a_E,
						const Vector<Real>& a_cdr_densities) const;

  /*!
    @brief Set the initial surface charge
  */
  virtual Real initial_sigma(const Real a_time, const RealVect& a_pos) const; 

protected:

  /*!
    @brief Noise function. Used for random initial conditions.
  */
  RefCountedPtr<perlin_if> m_perlin;
  Real m_noise_amp;      // Noise amplitude for Perlin noise
  RealVect m_noise_freq; // Noise frequency for Perlin noise
  Real m_noise_persist;  // Persistence for Perlin noise
  int m_noise_octaves;   // Number of octaves for Perlin noise

  int m_electron_idx;
  int m_N2plus_idx;
  int m_O2plus_idx;
  int m_N4plus_idx;
  int m_O4plus_idx;
  int m_O2plusN2_idx;
  int m_O2minus_idx;

  int m_photon1_idx;
  int m_photon2_idx;
  int m_photon3_idx;

  /*!
    @brief Simple cathode (electrode) boundary condition
  */
  bool m_simple_cathode;

  /*!
    @brief Thermal outflow on/off
  */
  bool m_thermal_outflow;

  /*!
    @brief Diffusive outflow on/off
  */
  bool m_diffusive_outflow;

  /*!
    @brief Background gas temperature
  */
  Real m_Tg;

  /*!
    @brief Background gas pressure
  */
  Real m_p;

  /*!
    @brief Quenching pressure
  */
  Real m_pq;

  /*!
    @brief Background gas density
  */
  Real m_N;

  /*!
    @brief O2 fraction of gas
  */
  Real m_O2frac;

  /*!
    @brief N2 fraction of gas
  */
  Real m_N2frac;

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
    @brief Photoionization efficiency
  */
  Real m_photoionization_efficiency;

  /*!
    @brief Excitation efficiency
  */
  Real m_excitation_efficiency;

  /*!
    @brief Ion mobility
  */
  Real m_ion_mobility;

  /*!
    @brief Electron-positive ion recombination
  */
  Real m_electron_recombination;

  /*!
    @brief Ion-ion recombination
  */
  Real m_ion_recombination;

  /*!
    @brief Electron detachmetn
  */
  Real m_electron_detachment;

  /*!
    @brief Compute the electron temperature
  */
  Real compute_electron_temperature(const Real a_EbyN) const;

  /*!
    @brief Compute electron mobility
  */
  Real compute_electron_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute electron diffusion
  */
  Real compute_electron_diffusion(const Real a_EbyN) const;

  /*!
    @brief Compute Townsend ionization coefficient
  */
  Real compute_townsend_ionization_N2(const Real a_EbyN) const;

  /*!
    @brief Compute Townsend ionization coefficient
  */
  Real compute_townsend_ionization_O2(const Real a_EbyN) const;

  /*!
    @brief Compute rate constant for \f$N_2^+ + N_2 + M \rightarrow N_4^+ + M\f$
  */
  Real compute_N2plus_N2_M_to_N4plus_M() const;

  /*!
    @brief Compute rate constant for \f$N_4^+ + O_2 \rightarrow O_2^+ + N_2 + N_2\f$
  */
  Real compute_N4plus_O2_to_O2plus_2N2() const;

  /*!
    @brief Compute rate constant for \f$N_2^+ + O_2 \rightarrow O_2^+ N_2\f$
  */
  Real compute_N2plus_O2_to_O2plus_N2() const;
  
  /*!
    @brief Compute rate constant for \f$O_2^+ + 2N_2 \rightarrow O_2^+N_2 N_2\f$
  */
  Real compute_O2plus_2N2_to_O2plusN2_N2() const;

  /*!
    @brief Compute rate constant for \f$O_2^+N2 + N_2 \rightarrow O_2^+ + 2N_2
  */
  Real compute_O2plusN2_N2_to_O2plus_2N2() const;

  /*!
    @brief Compute rate constant for \f$O_2^+N2 + O_2 \rightarrow O_4^+ + N_2
  */
  Real compute_O2plusN2_O2_to_O4plus_N2() const;

  /*!
    @brief Compute rate constant for \f$O_2^+ + O_2 + M\rightarrow O_4^+ + M\f$
  */
  Real compute_O2plus_O2_M_to_O4plus_M() const;

  /*!
    @brief Compute rate constant for the dissociative recombination \f$e + O_4^+ \rightarrow 2O_2\f$
  */
  Real compute_e_O4plus_to_2O2(const Real a_Te) const;

  /*!
    @brief Compute rate constant for the dissociative recombination \f$e + O_2^+ \rightarrow 2O\f$
  */
  Real compute_e_O2plus_to_2O(const Real a_Te) const;

  /*!
    @brief Compute rate constant for attachment \f$e + 2O_2 \rightarrow O_2^- + O_2\f$
  */
  Real compute_e_2O2_to_O2minus_O2(const Real a_Te) const;

  /*!
    @brief Compute rate constant for ion recombination \f$O_2^- + O_4^+ \rightarrow 3O_2\f$
  */
  Real compute_O2minus_O4plus_to_3O2() const;

  /*!
    @brief Compute rate constant for ion recombination \f$O_2^- + O_4^+ M \rightarrow 3O_2 + M\f$
  */
  Real compute_O2minus_O4plus_M_to_3O2_M() const;

  /*!
    @brief Compute rate constant for ion recombination \f$O_2^- + O_4^+ M \rightarrow 3O_2 + M\f$
  */
  Real compute_O2minus_O2plus_M_to_2O2_M() const;

  /*!
    @brief Common interface for cdr fluxes
  */
  Vector<Real> compute_cdr_fluxes(const Real&         a_time,
				  const RealVect&     a_pos,
				  const RealVect&     a_normal,
				  const RealVect&     a_E,
				  const Vector<Real>& a_cdr_densities,
				  const Vector<Real>& a_cdr_velocities,
				  const Vector<Real>& a_cdr_gradients,
				  const Vector<Real>& a_rte_fluxes,
				  const Vector<Real>& a_extrap_cdr_fluxes,
				  const Real&         a_townsend,
				  const Real&         a_quantum_efficiency,
				  const bool&         a_extrap_electrons) const;
};
#endif
