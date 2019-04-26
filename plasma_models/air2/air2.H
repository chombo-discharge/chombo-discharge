/*!
  @file   air2.H
  @brief  Declaration of the simplified two-species air model
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR2_
#define _AIR2_

#include "plasma_kinetics.H"
#include "lookup_table.H"
#include <random>

/*!
  @brief 2-species plasma model for 80/20 air at standard conditions
*/
class air2 : public plasma_kinetics {
public:

  // Forward declaration of cdr species
  class electron;
  class positive_ion;   

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
  air2();

  /*!
    @brief Destructor
  */
  virtual ~air2();


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
						const Real              a_kappa,
						const Real              a_dx,
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
  
  int m_photon1_idx;
  int m_photon2_idx;
  int m_photon3_idx;

  Real m_townsend2_electrode;            // 2nd Townsend on electrodes
  Real m_townsend2_dielectric;           // 2nd Townsend on dielectrics
  Real m_electrode_quantum_efficiency;   // Quantum efficiency on electrodes
  Real m_dielectric_quantum_efficiency;  // Quantum efficiency on dielectrics
  Real m_photoionization_efficiency;     // Photoionization efficiency
  Real m_excitation_efficiency;          // Excitation efficiency

  Real m_p;
  Real m_pq;

  // RNG stuff
  std::mt19937_64* m_rng;     // Engine
  int  m_seed;                // RNG seed
  Real m_poiss_exp_swap;      // Swap from Poisson to exponential distribution


private:

  inline Real get_alpha(const Real a_E) const;
  inline Real get_mobility(const Real a_E) const;

  inline Real stochastic_reaction(const Real a_S, const Real a_vol, const Real a_dt) const;
  
  /*!
    @brief Common function for dielectric and electrode fluxes, only difference is the secondary emission coefficients
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
				  const Real&         a_townsend2,
				  const Real&         a_quantum_efficiency) const;
};

#endif
