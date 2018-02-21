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
  @brief 11-species plasma model for air. 

  @details 

  Class options
  -------------

  Optional options:


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

};

#endif
