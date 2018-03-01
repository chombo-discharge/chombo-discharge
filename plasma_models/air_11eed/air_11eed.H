/*!
  @file   air_11eed.H
  @brief  Declaration of a BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR_11EED_
#define _AIR_11EED_

#include "plasma_kinetics.H"

/*!
  @brief 11-species plasma model for air (79/21 N2/O2 ratio) at 300K

  @details 

  This class uses functions from BOLSIG+ to compute ionization and transport coefficients for air at room temperature. BOLSIG+ has been pre-run for an electron energy range of \f$[10^{-3}\textrm{eV}, 20\textrm{eV}]\f$ which corresponds to electron temperatures in the range of 1K to 170000K. If you need transport data for other mixing ratios and gas temperatures, you should really, really re-run BOLSIG before running the kinetics in this program. 

  Class options
  -------------

  @todo If anyone finds the time, please improve the ion mobilities. 

*/
class air_11eed : public plasma_kinetics {
public:

  // Forward declaration of cdr species
  class eed;      // Electron energy distribution
  class electron;
  class N2;       
  class O2;       
  class N2plus;   
  class N4plus;   
  class O2plus;
  class O4plus;
  class O2plusN2;
  class O2minus;
  class Ominus;
  class O;


  // Forward declaration of photon excitations
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
  air_11eed();

  /*!
    @brief Destructor
  */
  virtual ~air_11eed();


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

  int m_eed_idx;
  int m_electron_idx;
  int m_N2_idx;
  int m_O2_idx;
  int m_N2plus_idx;
  int m_N4plus_idx;
  int m_O2plus_idx;
  int m_O4plus_idx;
  int m_O2plusN2_idx;
  int m_O2minus_idx;
  int m_Ominus_idx;
  int m_O_idx;
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
    @brief Compute electron mobility
  */
  Real compute_eed_mobility(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute electron mobility
  */
  Real compute_electron_mobility(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute N2plus mobility
  */
  Real compute_N2plus_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute N4plus mobility
  */
  Real compute_N4plus_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute O2plus mobility
  */
  Real compute_O2plus_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute O4plus mobility
  */
  Real compute_O4plus_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute O2plusN2 mobility
  */
  Real compute_O2plusN2_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute O2minus mobility
  */
  Real compute_O2minus_mobility(const Real a_EbyN) const;

  /*!
    @brief Compute Ominus mobility
  */
  Real compute_Ominus_mobility(const Real a_EbyN) const;

    /*!
    @brief Compute electron diffusion coefficient
  */
  Real compute_eed_diffco(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute electron diffusion coefficient
  */
  Real compute_electron_diffco(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute N2 diffusion coefficient
  */
  Real compute_N2_diffco() const;

  /*!
    @brief Compute O2 diffusion coefficient
  */
  Real compute_O2_diffco() const;

  /*!
    @brief Compute N2plus diffusion coefficient
  */
  Real compute_N2plus_diffco() const;

  /*!
    @brief Compute N4plus diffusion coefficient
  */
  Real compute_N4plus_diffco() const;

    /*!
    @brief Compute O2plus diffco
  */
  Real compute_O2plus_diffco() const;

  /*!
    @brief Compute O4plus diffco
  */
  Real compute_O4plus_diffco() const;

  /*!
    @brief Compute O2plusN2 diffco
  */
  Real compute_O2plusN2_diffco() const;

  /*!
    @brief Compute O2minus diffco
  */
  Real compute_O2minus_diffco() const;

  /*!
    @brief Compute Ominus diffco
  */
  Real compute_Ominus_diffco() const;

  /*!
    @brief Compute rate constant \f$e + N_2 \rightarrow 2e + N_2^+\f$
  */
  Real compute_electron_N2_impact_ionization(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute rate constant \f$e + O_2 \rightarrow 2e + O_2^+\f$
  */
  Real compute_electron_O2_impact_ionization(const Real a_energy, const Real a_N) const;

  /*!
    @brief Rate constant for \f$N_2^+ + N_2 + M \rightarrow N_4^+ + M\f$
  */
  Real compute_N2plus_N2_M_to_N4plus_M() const;

  /*!
    @brief Rate constant for \f$N_4^+ + O_2 \rightarrow O_2^+ + 2N_2\f$
  */
  Real compute_N4plus_O2_to_O2_2N2() const;

  /*!
    @brief Rate konstant for \f$N2^+ + O_2 \rightarrow O_2^+ + N2\f$
    @param[in] a_Tg Gas temperature
  */
  Real compute_N2plus_O2_to_O2plus_N2(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$O_2^+ + 2N_2 \rightarrow O_2^+N_2 + N_2\f$
    @param[in] a_Tg Gas temperature
  */
  Real compute_O2plus_2N2_to_O2plusN2_N2(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$O_2^+N_2 + N_2 \rightarrow O_2^+ + 2_N2\f$
  */
  Real compute_O2plusN2_N2_to_O2plus_2N2(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$O_2^+N_2 + N_2 \rightarrow O_2^+ + 2_N2\f$
  */
  Real compute_O2plusN2_O2_to_O4plus_N2() const;

  /*!
    @brief Rate constant for \f$O_2^+ O_2 + M \rightarrow O_4^+ + M\f$.
    @param[in] a_Tg Gas temperature
  */
  Real compute_O2plus_O2_M_to_O4plus_M(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$e + O4plus \rightarrow 2O_2\f$
    @param[in] a_Tg Gas temperature
  */
  Real compute_e_O4plus_to_2O2(const Real a_Te) const;

  /*!
    @brief Rate constant for \f$e + O_2^+ \rightarrow O_2\f$
    @param[in] a_Tg Gas temperature
  */
  Real compute_e_O2plus_to_O2(const Real a_Te) const;

  /*!
    @brief Rate constant for \f$e + 2O_2 \rightarrow O_2^- + O_2\f$
    @param[in] a_Tg Gas temperature
  */
  Real compute_e_2O2_to_O2minus_O2(const Real a_Te) const;

  /*!
    @brief Rate constant for \f$O_2^- + O_4^+ \rightarrow 3O_2\f$
   */
  Real compute_O2minus_O4plus_to_3O2() const;

  /*!
    @brief Rate constant for \f$O_2^- + O_4^+ + M \rightarrow 3O_2 + M\f$
  */
  Real compute_O2minus_O4plus_M_to_3O2_M(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$O_2^- + O_2^+ + M \rightarrow 2O_2 + M\f$
  */
  Real compute_O2minus_O2plus_M_to_2O2_M(const Real a_Tg) const;

  /*!
    @brief Rate constant for \f$e + O_2 \rightarrow e + 2O\f$. 
    @details This is the contribution at 5.58eV
  */
  Real compute_e_O2_to_e_2O_c1(const Real a_energy, const Real a_N) const;

  /*!
    @brief Rate constant for \f$e + O_2 \rightarrow e + 2O\f$
    @details This is the contribution at 8.4eV
  */
  Real compute_e_O2_to_e_2O_c2(const Real a_energy, const Real a_N) const;

  /*!
    @brief Rate constant for \f$e + O_2 \rightarrow O^- + O\f$
  */
  Real compute_e_O2_to_Ominus_O(const Real a_energy, const Real a_N) const;

  /*!
    @brief Rate constant for \f$O^- + O_2^+ \rightarrow O + O_2\f$
  */
  Real compute_Oplus_O2_to_O_O2(const Real a_Tg) const;

  /*!
    @brief Compute the energy loss \f$e + N_2 \rightarrow e + N_2\f$
  */
  Real compute_e_N2_to_e_N2(const Real a_energy, const Real a_N) const;

  /*!
    @brief Compute the energy loss \f$e + N_2 \rightarrow e + N_2\f$
  */
  Real compute_e_O2_to_e_O2(const Real a_energy, const Real a_N) const;

  /*!
    @brief Energy loss for impact ionization of nitrogen. Returns 15.6 eV
  */
  Real compute_e_N2_ionization_loss() const;

  /*!
    @brief Energy loss for impact ionization of oxygen. Returns 12.07 eV
  */
  Real compute_e_O2_ionization_loss() const;

  /*!
    @brief Compute energy loss for the \f$e + O_2\rightarrow e + 2O\f$ dissociation.
    @details This is the one at 5.58eV
  */
  Real compute_e_O2_dissociation_loss_c1() const;

    /*!
    @brief Compute energy loss for the \f$e + O_2\rightarrow e + 2O\f$ dissociation.
    @details This is the one at 8.4eV
  */
  Real compute_e_O2_dissociation_loss_c2() const;

  /*!
    @brief Compute the energy loss for the \f$e + O_2 \rightarrow O^- + O\f$ reaction. 
  */
  Real compute_e_O2_dissociative_attachment_loss() const;

  /*!
    @brief Compute the energy loss term for the \f$e + O_2\f$ inelastic collisisons
  */
  Real compute_e_O2_scattering_loss() const;

  /*!
    @brief Compute the energy loss term for the \f$e + N_2\f$ inelastic collisisons
  */
  Real compute_e_N2_scattering_loss() const;

  /*!
    @brief Compute fluxes at the anode
  */
  Vector<Real> compute_anode_fluxes(const Real&         a_time,
				    const RealVect&     a_pos,
				    const RealVect&     a_normal,
				    const RealVect&     a_E,
				    const Vector<Real>& a_cdr_densities,
				    const Vector<Real>& a_cdr_velocities,
				    const Vector<Real>& a_rte_fluxes,
				    const Vector<Real>& a_extrap_cdr_fluxes) const;
};

#endif
