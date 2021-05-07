/*!
  @file   air3_zheleznyak.H
  @brief  3-species and 8-photon model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR3_ZHELEZNYAK_
#define _AIR3_ZHELEZNYAK_

#include "cdr_plasma_physics.H"
#include "lookup_table.H"

#include <random>

namespace ChomboDischarge {
  namespace physics {
    namespace cdr_plasma {

      class air3_zheleznyak : public cdr_plasma_physics {
      public:

	class electron;
	class M_plus;
	class M_minus;

	class uv_photon;


	air3_zheleznyak();
	~air3_zheleznyak();
  

	void advance_reaction_network(Vector<Real>&          a_particle_sources,
				      Vector<Real>&          a_photon_sources,
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

	Real initial_sigma(const Real a_time, const RealVect a_pos) const;

	Real compute_alpha(const RealVect a_E) const;

      protected:

	static std::string s_bolsig_mobility;
	static std::string s_bolsig_diffco;
	static std::string s_bolsig_alpha;
	static std::string s_bolsig_eta;
	static std::string s_bolsig_b1_exc;
	static std::string s_bolsig_c4_exc;

	std::string m_transport_file;

	enum chemistryAlgorithm {
	  euler = 0,
	  rk2   = 1,
	  rk4   = 2
	};

	chemistryAlgorithm m_chemistryAlgorithm;

	// Tables built from BOLSIG+
	lookup_table m_e_mobility;
	lookup_table m_e_diffco;
	lookup_table m_e_alpha;
	lookup_table m_e_eta;
	lookup_table m_b1_exc;
	lookup_table m_c4_exc;

	Vector<int> m_wallbc; // Wall boundary conditions

	bool m_alpha_corr;          // Alpha correction
	bool m_mobile_electrons;    // Mobile electrons or not
	bool m_mobile_ions;         // Mobile ions or not
	bool m_diffusive_ions;      // Diffusive ions or not
	bool m_diffusive_electrons; // Diffusive electrons or not
  
	Real m_rng_seed;
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
	int m_phot_idx; // Photon idx

	Real m_townsend2_electrode;           // 2nd Townsend coefficient on electrodes
	Real m_townsend2_dielectric;          // 2nd Townsend coefficient on dielectrics
	Real m_electrode_quantum_efficiency;  // Quantum yield on electrodes
	Real m_dielectric_quantum_efficiency; // Quantum yield on dielectircs
  
	Real m_N;      // Neutral number density
	Real m_p;      // Pressure
	Real m_pq;     // Quenching pressure
	Real m_T;      // Temperature
	Real m_O2frac; // O2 fraction
	Real m_N2frac; // N2 fraction
	Real m_pO2;    // Partial O2 pressure

	Real m_factor;

	Real m_chemistry_dt;

	Real m_ion_mobility;  // Ion mobility
	Real m_ion_diffusion; // Ion diffusion coefficient

	// Functions for creating random Gaussian distributions
	RealVect random_gaussian();
	RealVect random_direction();
#if CH_SPACEDIM==2
	RealVect random_direction2D();
#elif CH_SPACEDIM==3
	RealVect random_direction3D();
#endif

	int poisson_reaction(const Real a_propensity, const Real a_dt) const;

	void read_file_entries(lookup_table& a_table, const std::string a_string);
	void parse_transport_file();
	void parse_transport();
	void parse_gas_params();
	void parse_electron_mobility();
	void parse_electron_diffco();
	void parse_alpha();
	void parse_eta();
	void parse_see();
	void parse_domain_bc();
	void parse_initial_particles();
	void parse_chemistry();
	void instantiate_species();
	void init_rng();
	void add_uniform_particles(List<Particle>& a_particles, const int a_num, const Real a_weight);
	void add_gaussian_particles(List<Particle>& a_particles,
				    const int       a_num,
				    const Real      a_weight,
				    const Real      a_rad,
				    const RealVect  a_center);


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

	void advance_chemistry_euler(Vector<Real>&          a_particle_sources,
				     Vector<Real>&          a_photon_sources,
				     Vector<Real>&          a_particle_densities,
				     const Vector<RealVect> a_particle_gradients,
				     const Vector<Real>     a_photon_densities,
				     const RealVect         a_E,
				     const RealVect         a_pos,
				     const Real             a_dx,
				     const Real             a_dt,
				     const Real             a_time,
				     const Real             a_kappa) const;

	inline Real excitation_rates(const Real a_E) const;
	inline Real sergey_factor(const Real a_O2frac) const;
      };
    }
  }
}
#endif
