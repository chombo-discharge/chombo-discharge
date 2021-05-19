/*!
  @file   morrow_fhd.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _MORROW_FHD_
#define _MORROW_FHD_

#include "simple_kinetics.H"

#include <random>

/*!
  @brief Class that implements the fluctuating Morrow-Lowke discharge model for air.
*/
class morrow_fhd : public simple_kinetics {
public:

  morrow_fhd();
  ~morrow_fhd();

  Vector<Real> compute_cdr_diffusion_coefficients(const Real&         a_time,
						  const RealVect&     a_pos,
						  const RealVect&     a_E,
						  const Vector<Real>& a_cdr_densities) const;
  
  Vector<RealVect> compute_cdr_velocities(const Real&         a_time,
					  const RealVect&     a_pos,
					  const RealVect&     a_E,
					  const Vector<Real>& a_cdr_densities) const;
  
  Vector<Real> compute_cdr_source_terms(const Real              a_time,
					const Real              a_kappa,
					const Real              a_dx,
					const RealVect&         a_pos,
					const RealVect&         a_E,
					const RealVect&         a_gradE,
					const Vector<Real>&     a_cdr_densities,
					const Vector<Real>&     a_rte_densities,
					const Vector<RealVect>& a_grad_cdr) const;

  Vector<Real> compute_cdr_domain_fluxes(const Real&           a_time,
					 const RealVect&       a_pos,
					 const int&            a_dir,
					 const Side::LoHiSide& a_side,
					 const RealVect&       a_E,
					 const Vector<Real>&   a_cdr_densities,
					 const Vector<Real>&   a_cdr_velocities,
					 const Vector<Real>&   a_cdr_gradients,
					 const Vector<Real>&   a_rte_fluxes,
					 const Vector<Real>&   a_extrap_cdr_fluxes) const;
  
  Vector<Real> compute_cdr_electrode_fluxes(const Real&         a_time,
					    const RealVect&     a_pos,
					    const RealVect&     a_normal,
					    const RealVect&     a_E,
					    const Vector<Real>& a_cdr_densities,
					    const Vector<Real>& a_cdr_velocities,
					    const Vector<Real>& a_cdr_gradients,
					    const Vector<Real>& a_rte_fluxes,
					    const Vector<Real>& a_extrap_cdr_fluxes) const;

  Vector<Real> compute_cdr_dielectric_fluxes(const Real&         a_time,
					     const RealVect&     a_pos,
					     const RealVect&     a_normal,
					     const RealVect&     a_E,
					     const Vector<Real>& a_cdr_densities,
					     const Vector<Real>& a_cdr_velocities,
					     const Vector<Real>& a_cdr_gradients,
					     const Vector<Real>& a_rte_fluxes,
					     const Vector<Real>& a_extrap_cdr_fluxes) const;

  Vector<Real> compute_rte_source_terms(const Real&         a_time,
					const Real&         a_kappa,
					const Real&         a_dx,
					const RealVect&     a_pos,
					const RealVect&     a_E,
					const Vector<Real>& a_cdr_densities) const;

  Vector<Real> compute_conductor_fluxes(const Vector<Real>& a_extrapolated_fluxes,
					const Vector<Real>& a_ion_densities,
					const Vector<Real>& a_ion_velocities,
					const Vector<Real>& a_Photon_fluxes,
					const RealVect&     a_E,
					const RealVect&     a_pos,
					const RealVect&     a_normal,
					const Real&         a_time) const;

  Vector<Real> compute_dielectric_fluxes(const Vector<Real>& a_extrapolated_fluxes,
					 const Vector<Real>& a_ion_densities,
					 const Vector<Real>& a_ion_velocities,
					 const Vector<Real>& a_Photon_fluxes,
					 const RealVect&     a_E,
					 const RealVect&     a_pos,
					 const RealVect&     a_normal,
					 const Real&         a_time) const;

  Real initial_sigma(const Real a_time, const RealVect& a_pos) const;


  Vector<int> m_wallBc; // Wall boundary conditions

  bool m_fhd; // Using FHD or not
  
  Real m_seed;
  Real m_cutoff_poisson;
  Real m_poiss_exp_swap;
  std::mt19937_64* m_rng; 

  int m_nelec_idx;      // Electron index
  int m_nplus_idx;      // Positive ion index
  int m_nminu_idx;      // Negative ion index
  int m_Photon1_idx;    // UV Photon index

  Real m_background_rate;               // Background ionization
  Real m_townsend2_electrode;           // 2nd Townsend coefficient on electrodes
  Real m_townsend2_dielectric;          // 2nd Townsend coefficient on dielectrics
  Real m_electrode_quantum_efficiency;  // Quantum yield on electrodes
  Real m_dielectric_quantum_efficiency; // Quantum yield on dielectircs
  
  Real m_N;      // Neutral number density
  Real m_p;      // Pressure
  Real m_T;      // Temperature
  Real m_fracO2; // O2 fraction
  Real m_fracN2; // N2 fraction
  Real m_pO2;    // 

  Real m_pq;         // Quenching pressure
  Real m_exc_eff;    // Excitation efficiency
  Real m_photoi_eff; // Photoionization efficiency

  // Common function for computing CDR fluxes
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

  inline Real stochastic_reaction(const Real a_S, const Real a_vol, const Real a_dt) const;
  
  inline Real compute_alpha(const RealVect& a_E) const;
  inline Real compute_eta(const RealVect& a_E)   const;
  inline Real compute_eta2(const RealVect& a_E)  const;
  inline Real compute_eta3(const RealVect& a_E)  const;
  inline Real compute_beta(const RealVect& a_E)  const;

  inline Real compute_De(const RealVect& a_E) const;

  inline RealVect compute_ve(const RealVect& a_E) const;
  inline RealVect compute_vp(const RealVect& a_E) const;
  inline RealVect compute_vn(const RealVect& a_E) const;

  class electron : public species {
  public:
    electron();
    ~electron();

    Real initialData(const RealVect a_pos, const Real a_time) const;

    Real m_uniform_density;
    Real m_seed_density;
    Real m_seed_radius;
    
    RealVect m_seed_pos;
  };

  class positive_species : public species {
  public:
    positive_species();
    ~positive_species();

    Real initialData(const RealVect a_pos, const Real a_time) const;

    Real m_uniform_density;
    Real m_seed_density;
    Real m_seed_radius;
    
    RealVect m_seed_pos;
  };

  class negative_species : public species {
  public:
    negative_species();
    ~negative_species();
    
    Real initialData(const RealVect a_pos, const Real a_time) const;
    Real m_uniform_density;
    
    Real m_seed_density;
    Real m_seed_radius;

    RealVect m_seed_pos;
  };

  class uv_Photon : public Photon_group {
  public:
    uv_Photon();
    ~uv_Photon();


    Real getKappa(const RealVect a_pos) const;
    Real get_random_kappa() const;

    inline Real get_pO2() const { return m_pO2; }

  protected:
    Real m_pO2;
    Real m_f1;
    Real m_f2;
    Real m_K1;
    Real m_K2;

    Real m_seed;
    std::mt19937_64* m_rng; 
    std::uniform_real_distribution<Real>* m_udist01;
  };

#include "CD_NamespaceFooter.H"

#endif
