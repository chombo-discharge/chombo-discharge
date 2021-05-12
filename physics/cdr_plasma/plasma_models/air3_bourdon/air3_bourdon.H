#ifndef _AIR3_BOURDON_
#define _AIR3_BOURDON_

#include "cdr_plasma_physics.H"
#include "lookup_table.H"

#include <random>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {

    class air3_bourdon : public cdr_plasma_physics {
    public:

      class electron;
      class M_plus;
      class M_minus;

      class photon_one;
      class photon_two;
      class photon_three;

      air3_bourdon();
      ~air3_bourdon();
  

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
				    const Real             a_kappa) const override;


      Vector<Real> compute_cdr_diffusion_coefficients(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const override;
  
      Vector<RealVect> compute_cdr_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const override;
  
      Vector<Real> compute_cdr_domain_fluxes(const Real           a_time,
					     const RealVect       a_pos,
					     const int            a_dir,
					     const Side::LoHiSide a_side,
					     const RealVect       a_E,
					     const Vector<Real>   a_cdr_densities,
					     const Vector<Real>   a_cdr_velocities,
					     const Vector<Real>   a_cdr_gradients,
					     const Vector<Real>   a_rte_fluxes,
					     const Vector<Real>   a_extrap_cdr_fluxes) const override;
  
      Vector<Real> compute_cdr_electrode_fluxes(const Real         a_time,
						const RealVect     a_pos,
						const RealVect     a_normal,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities,
						const Vector<Real> a_cdr_velocities,
						const Vector<Real> a_cdr_gradients,
						const Vector<Real> a_rte_fluxes,
						const Vector<Real> a_extrap_cdr_fluxes) const override;

      Vector<Real> compute_cdr_dielectric_fluxes(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_normal,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities,
						 const Vector<Real> a_cdr_velocities,
						 const Vector<Real> a_cdr_gradients,
						 const Vector<Real> a_rte_fluxes,
						 const Vector<Real> a_extrap_cdr_fluxes) const override;

      Real initial_sigma(const Real a_time, const RealVect a_pos) const override;

      Real compute_alpha(const RealVect a_E) const override;

    protected:

      static std::string s_bolsig_mobility;
      static std::string s_bolsig_diffco;
      static std::string s_bolsig_alpha;
      static std::string s_bolsig_eta;

      std::string m_transport_file;

      // Tables built from BOLSIG+
      lookup_table m_e_mobility;
      lookup_table m_e_diffco;
      lookup_table m_e_alpha;
      lookup_table m_e_eta;

      Vector<int> m_wallbc; // Wall boundary conditions

      bool m_alpha_corr;          // Alpha correction
      bool m_mobile_electrons;    // Mobile electrons or not
      bool m_mobile_ions;         // Mobile ions or not
      bool m_diffusive_ions;      // Diffusive ions or not
      bool m_diffusive_electrons; // Diffusive electrons or not
  
      int m_uniform_entries; // # of entries in uniform-made BOLSIG tables

      // Indices for various species
      int m_elec_idx; // Electrons
      int m_plus_idx; // Positive ions
      int m_minu_idx; // Negative ions
      int m_pho1_idx; // Photon idx
      int m_pho2_idx; // Photon idx
      int m_pho3_idx; // Photon idx

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

      Real m_ion_mobility;  // Ion mobility
      Real m_ion_diffusion; // Ion diffusion coefficient

      Real m_photoexc_eff; // Photoexcitation efficiency
      Real m_photoi_eff;    // Photoionizationn efficiency


      void read_file_entries(lookup_table& a_table, const std::string a_string);
      void parseDomain_bc();
      void instantiate_species();


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
  }
}
#include "CD_NamespaceFooter.H"
#endif
