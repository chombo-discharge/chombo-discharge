/*!
  @file   air6_mc8_bourdon.H
  @brief  6-species (3/3 charged/excited) and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR6_MC8_
#ifndef _AIR6_MC8_

#include "network_kinetics.H"
#include "lookup_table.H"

class air6_mc8 : public network_kinetics {
public:

  class electron;
  class M_plus;
  class M_minus;
  class N2_c4_v0;
  class N2_c4_v1;
  class N2_b1_v1;

  class phot_c4v0_X1v0;
  class phot_c4v0_X1v1;
  class phot_c4v1_X1v0;
  class phot_c4v1_X1v1;
  class phot_c4v1_X1v2;
  class phot_c4v1_X1v3;
  class phot_b1v1_X1v0;
  class phot_b1v1_X1v1;

  void advance_reaction_network(Vector<Real>&           a_particle_sources,
				Vector<Real>&           a_photon_sources,
				const Vector<Real>     a_particle_densities,
				const Vector<RealVect> a_particle_gradients,
				const Vector<Real>     a_photon_densities,
				const RealVect         a_E,
				const RealVect         a_pos,
				const Real             a_dx,
				const Real             a_dt,
				const Real             a_time,
				const Real             a_kappa) const;


  Vector<Real> compute_cdr_diffusion_coefficients(const Real         a_time,
						  const RealVect     a_pos,
						  const RealVect     a_E,
						  const Vector<Real> a_cdr_densities) const;
  
  Vector<RealVect> compute_cdr_velocities(const Real         a_time,
					  const RealVect     a_pos,
					  const RealVect     a_E,
					  const Vector<Real> a_cdr_densities) const;
  
  Vector<Real> compute_cdr_domain_fluxes(const Real           a_time,
					 const RealVect       a_pos,
					 const int            a_dir,
					 const Side::LoHiSide a_side,
					 const RealVect       a_E,
					 const Vector<Real>   a_cdr_densities,
					 const Vector<Real>   a_cdr_velocities,
					 const Vector<Real>   a_cdr_gradients,
					 const Vector<Real>   a_rte_fluxes,
					 const Vector<Real>   a_extrap_cdr_fluxes) const;
  
  Vector<Real> compute_cdr_electrode_fluxes(const Real         a_time,
					    const RealVect     a_pos,
					    const RealVect     a_normal,
					    const RealVect     a_E,
					    const Vector<Real> a_cdr_densities,
					    const Vector<Real> a_cdr_velocities,
					    const Vector<Real> a_cdr_gradients,
					    const Vector<Real> a_rte_fluxes,
					    const Vector<Real> a_extrap_cdr_fluxes) const;

  Vector<Real> compute_cdr_dielectric_fluxes(const Real         a_time,
					     const RealVect     a_pos,
					     const RealVect     a_normal,
					     const RealVect     a_E,
					     const Vector<Real> a_cdr_densities,
					     const Vector<Real> a_cdr_velocities,
					     const Vector<Real> a_cdr_gradients,
					     const Vector<Real> a_rte_fluxes,
					     const Vector<Real> a_extrap_cdr_fluxes) const;

  Vector<Real> compute_dielectric_fluxes(const Vector<Real> a_extrapolated_fluxes,
					 const Vector<Real> a_ion_densities,
					 const Vector<Real> a_ion_velocities,
					 const Vector<Real> a_photon_fluxes,
					 const RealVect     a_E,
					 const RealVect     a_pos,
					 const RealVect     a_normal,
					 const Real         a_time) const;

  Real initial_sigma(const Real a_time, const RealVect a_pos) const;

protected:

  static std::string s_bolsig_mobility;
  static std::string s_bolsig_diffco;
  static std::string s_bolsig_alpha;
  static std::string s_bolsig_eta;

  enum source_comp{
    ssa = 0,
    tau = 1,
    rre = 2
  };
  
  source_comp m_scomp;

  // Tables built from BOLSIG+
  lookup_table m_e_mobility;
  lookup_table m_e_alpha;
  lookup_table m_e_eta;

  Vector<int> m_wallbc; // Wall boundary conditions

  bool m_mobile_ions;         // Mobile ions or not
  bool m_diffusive_ions;      // Diffusive ions or not
  bool m_diffusive_electrons; // Diffusive electrons or not
  
  Real m_seed;
  Real m_cutoff_poisson;
  Real m_poiss_exp_swap;

  // Random number generator and some distributions
  std::mt19937_64*                      m_rng; 
  std::uniform_real_distribution<Real>* m_udist01;
  std::uniform_real_distribution<Real>* m_udist11;
  std::normal_distribution<Real>        m_gauss;

  int m_uniform_entries; // # of entries in uniform-made BOLSIG tables

  // Indices for various species
  int m_elec_idx; // Electrons
  int m_plus_idx; // Positive ions
  int m_minu_idx; // Negative ions
  int m_c4v0_idx; // N2 c4(v=0) excited state
  int m_c4v1_idx; // N2 c4(v=1) excited state
  int m_b1v1_idx; // N2 b1(v=0) excited state

  // Indicates for various photons
  int m_c4v0_X1v0_idx; // Transition c4(v=0) -> X1(v=0)
  int m_c4v0_X1v1_idx; // Transition c4(v=0) -> X1(v=0)
  int m_c4v1_X1v0_idx; // Transition c4(v=1) -> X1(v=0)
  int m_c4v1_X1v1_idx; // Transition c4(v=1) -> X1(v=1)
  int m_c4v1_X1v2_idx; // Transition c4(v=1) -> X1(v=2)
  int m_c4v1_X1v3_idx; // Transition c4(v=1) -> X1(v=3)
  int m_b1v1_X1v0_idx; // Transition b1(v=1) -> X1(v=0)
  int m_b1v1_X1v1_idx; // Transition b1(v=1) -> X1(v=1)

  Real m_townsend2_electrode;           // 2nd Townsend coefficient on electrodes
  Real m_townsend2_dielectric;          // 2nd Townsend coefficient on dielectrics
  Real m_electrode_quantum_efficiency;  // Quantum yield on electrodes
  Real m_dielectric_quantum_efficiency; // Quantum yield on dielectircs
  
  Real m_N;      // Neutral number density
  Real m_p;      // Pressure
  Real m_T;      // Temperature
  Real m_fracO2; // O2 fraction
  Real m_fracN2; // N2 fraction
  Real m_pO2;    // Partial O2 pressure

  // Lifetimes for various transitions
  Real m_c4v0_X1v0_tau_r; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v0_X1v1_tau_r; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v1_X1v0_tau_r; // Transition c4(v=1) -> X1(v=0)
  Real m_c4v1_X1v1_tau_r; // Transition c4(v=1) -> X1(v=1)
  Real m_c4v1_X1v2_tau_r; // Transition c4(v=1) -> X1(v=2)
  Real m_c4v1_X1v3_tau_r; // Transition c4(v=1) -> X1(v=3)
  Real m_b1v1_X1v0_tau_r; // Transition b1(v=1) -> X1(v=0)
  Real m_b1v1_X1v1_tau_r; // Transition b1(v=1) -> X1(v=1)

  // Quenching lifetimes for various transitions
  Real m_c4v0_X1v0_tau_q; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v0_X1v1_tau_q; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v1_X1v0_tau_q; // Transition c4(v=1) -> X1(v=0)
  Real m_c4v1_X1v1_tau_q; // Transition c4(v=1) -> X1(v=1)
  Real m_c4v1_X1v2_tau_q; // Transition c4(v=1) -> X1(v=2)
  Real m_c4v1_X1v3_tau_q; // Transition c4(v=1) -> X1(v=3)
  Real m_b1v1_X1v0_tau_q; // Transition b1(v=1) -> X1(v=0)
  Real m_b1v1_X1v1_tau_q; // Transition b1(v=1) -> X1(v=1)

  // Photoionization efficiency for various transitions
  Real m_c4v0_X1v0_eff; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v0_X1v1_eff; // Transition c4(v=0) -> X1(v=0)
  Real m_c4v1_X1v0_eff; // Transition c4(v=1) -> X1(v=0)
  Real m_c4v1_X1v1_eff; // Transition c4(v=1) -> X1(v=1)
  Real m_c4v1_X1v2_eff; // Transition c4(v=1) -> X1(v=2)
  Real m_c4v1_X1v3_eff; // Transition c4(v=1) -> X1(v=3)
  Real m_b1v1_X1v0_eff; // Transition b1(v=1) -> X1(v=0)
  Real m_b1v1_X1v1_eff; // Transition b1(v=1) -> X1(v=1)

  Real m_pq;         // Quenching pressure
  Real m_exc_eff;    // Excitation efficiency
  Real m_photoi_eff; // Photoionization efficiency
  Real m_particle_w; // Particle weights
  Real m_exc_tau;    // Lifetime for excited state
  Real m_quench_tau; // Lifetime for excited state

  Real m_ion_mobility;  // Ion mobility
  Real m_ion_diffusion; // Ion diffusion coefficient

  void read_file_entries(lookup_table& a_table, const std::string a_string);
  void parse_transport_file();
  void parse_transport();
  void read_gas_params();
  void parse_electron_mobility();
  void parse_alpha();
  void parse_eta();
  void parse_see();
  void parse_domain_bc();
  void parse_reaction_settings();
  void parse_initial_particles();
  void instantiate_species();
  void init_rng();
  void add_uniform_pairs(List<Particle>& a_particles);
  void add_gaussian_pairs(List<Particle>& a_particles);
  void add_uniform_excited(List<Particle>& a_particles);
  void add_gaussian_excited(List<Particle>& a_particles);

  void parse_photoi();
  void parse_dissociation_rates();
  void parse_radiative_rates();
  void parse_photoi_eff();

private:

  void read_file_entries(lookup_table& a_table, const std::string a_string);
  
  /*!
    @brief Common function for dielectric and electrode fluxes, only difference is the secondary emission coefficients
  */
  Vector<Real> compute_cdr_fluxes(const Real         a_time,
				  const RealVect     a_pos,
				  const RealVect     a_normal,
				  const RealVect     a_E,
				  const Vector<Real> a_cdr_densities,
				  const Vector<Real> a_cdr_velocities,
				  const Vector<Real> a_cdr_gradients,
				  const Vector<Real> a_rte_fluxes,
				  const Vector<Real> a_extrap_cdr_fluxes,
				  const Real         a_townsend2,
				  const Real         a_quantum_efficiency) const;
};

#endif
