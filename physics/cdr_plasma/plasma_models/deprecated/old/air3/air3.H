/*!
  @file   air3.H
  @brief  Declaration of a 3-species BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR3_
#define _AIR3_

#include "plasma_kinetics.H"

/*!
  @brief 3-species BOLSIG-based transport model for air at 1atm and 300K

  @details

  Class options
  -------------

  air3.quenching_pressure            = 0.03947               # Quenching pressure for photo-emission
  air3.electron_diffusion            = true                  # Turn off electron diffusion (debugging only)
  air3.ion_mobility                  = 2.E-4                 # Ion mobility
  air3.excitation_efficiency         = 0.6                   # Impact excitation efficiency
  air3.photoionization_efficiency    = 0.1                   # Photo-ionization efficiency
  air3.electrode_townsend2           = 1.E-2                 # Townsend coefficient on electrodes
  air3.electrode_quantum_efficiency  = 1.E-2                 # Quantum efficiency on electrodes
  air3.dielectric_townsend2          = 1.E-1                 # Townsend coefficient on dielectrics
  air3.dielectric_quantum_efficiency = 1.E-1                 # Quantum efficiency on dielectrics
  air3.initial_ionization            = 1.E10                 # Initial background ionization
  air3.photon1_A_coeff               = 1.12E-4               # Parameters from Bourdon et. al photoionization model
  air3.photon1_lambda_coeff          = 4.15E-2               # Parameters from Bourdon et. al photoionization model
  air3.photon2_A_coeff               = 2.88E-2               # Parameters from Bourdon et. al photoionization model
  air3.photon2_lambda_coeff          = 1.09E-1               # Parameters from Bourdon et. al photoionization model
  air3.photon3_A_coeff               = 2.76E-1               # Parameters from Bourdon et. al photoionization model
  air3.photon3_lambda_coeff          = 6.69E-1               # Parameters from Bourdon et. al photoionization model
  air3.electron_recombination        = 5.E-14                # Electron-ion recombination rate
  air3.ion_recombination             = 2.E-12                # Ion-ion recombination rate
  air3.electron_detachment           = 1.E-18                # Electron detachment rate

*/
class air3 : public plasma_kinetics {
public:

  class electron;
  class positive_species;
  class negative_species;

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
  air3();

  /*!
    @brief Destructor
  */
  virtual ~air3();

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

  int m_electron_idx;
  int m_positive_idx;
  int m_negative_idx;
  int m_photon1_idx;
  int m_photon2_idx;
  int m_photon3_idx;

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
  Real compute_townsend_ionization(const Real a_EbyN) const;

  /*!
    @brief Compute Townsend ionization coefficient
  */
  Real compute_townsend_attachment(const Real a_EbyN) const;

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
				  const Real&         a_quantum_efficiency) const;
#include "CD_NamespaceFooter.H"
#endif
