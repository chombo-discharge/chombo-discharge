/*!
  @file   air7_zheleznyak.H
  @brief  7-species model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR7_ZHELEZNYAK_
#define _AIR7_ZHELEZNYAK_

#include <CD_CdrPlasmaPhysics.H>
#include "CD_LookupTable.H"

#include <random>

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace CdrPlasma {
  
    class air7_zheleznyak : public CdrPlasmaPhysics {
    public:

      class electron;
      class N2plus;
      class O2plus;
      class N4plus;
      class O4plus;
      class O2plusN2;
      class O2minus;

      class uv_Photon;


      air7_zheleznyak();
      ~air7_zheleznyak();
  

      void advanceReactionNetwork(Vector<Real>&          a_particle_sources,
				    Vector<Real>&          a_Photon_sources,
				    const Vector<Real>     a_particle_densities,
				    const Vector<RealVect> a_particle_gradients,
				    const Vector<Real>     a_Photon_densities,
				    const RealVect         a_E,
				    const RealVect         a_pos,
				    const Real             a_dx,
				    const Real             a_dt,
				    const Real             a_time,
				    const Real             a_kappa) const;


      Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const;
  
      Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const;
  
      Vector<Real> computeCdrDomainFluxes(const Real           a_time,
					     const RealVect       a_pos,
					     const int            a_dir,
					     const Side::LoHiSide a_side,
					     const RealVect       a_E,
					     const Vector<Real>   a_cdr_densities,
					     const Vector<Real>   a_cdr_velocities,
					     const Vector<Real>   a_cdr_gradients,
					     const Vector<Real>   a_rte_fluxes,
					     const Vector<Real>   a_extrap_cdr_fluxes) const;
  
      Vector<Real> computeCdrElectrodeFluxes(const Real         a_time,
						const RealVect     a_pos,
						const RealVect     a_normal,
						const RealVect     a_E,
						const Vector<Real> a_cdr_densities,
						const Vector<Real> a_cdr_velocities,
						const Vector<Real> a_cdr_gradients,
						const Vector<Real> a_rte_fluxes,
						const Vector<Real> a_extrap_cdr_fluxes) const;

      Vector<Real> computeCdrDielectricFluxes(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_normal,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities,
						 const Vector<Real> a_cdr_velocities,
						 const Vector<Real> a_cdr_gradients,
						 const Vector<Real> a_rte_fluxes,
						 const Vector<Real> a_extrap_cdr_fluxes) const;

      Real initialSigma(const Real a_time, const RealVect a_pos) const;

      Real computeAlpha(const RealVect a_E) const;

    protected:

      static std::string s_bolsig_mobility;
      static std::string s_bolsig_diffco;
      static std::string s_bolsig_alpha;
      static std::string s_bolsig_eta;
      static std::string s_bolsig_alphaN2;
      static std::string s_bolsig_alphaO2;

      static std::string s_bolsig_energy;
  
      std::string m_transport_file;

      enum chemistryAlgorithm {
	euler = 0,
	rk2   = 1,
	rk4   = 2
      };

      chemistryAlgorithm m_chemistryAlgorithm;

      // Tables built from BOLSIG+
      LookupTable m_e_mobility;    // Mobility
      LookupTable m_e_diffco;      // Diffusion coefficient
      LookupTable m_e_alphaN2;     // Ionization rate
      LookupTable m_e_alphaO2;     // Ionization rate
      LookupTable m_e_alpha;       // Townsend alpha
      LookupTable m_e_eta;         // Townsend eta
      LookupTable m_e_temperature; // Electron temperature (K)


      Vector<int> m_wallBc; // Wall boundary conditions

      bool m_alpha_corr;          // Alpha correction
      bool m_isMobile_electrons;    // Mobile electrons or not
      bool m_isMobile_ions;         // Mobile ions or not
      bool m_isDiffusive_ions;      // Diffusive ions or not
      bool m_isDiffusive_electrons; // Diffusive electrons or not
  
      Real m_rng_seed;
      Real m_poiss_exp_swap;

      // Random number generator and some distributions
      std::mt19937_64*                      m_rng; 
      std::uniform_real_distribution<Real>* m_udist01;
      std::uniform_real_distribution<Real>* m_udist11;
      std::normal_distribution<Real>        m_gauss;

      int m_uniform_entries; // # of entries in uniform-made BOLSIG tables

      // Indices for various species
      int m_elec_idx;     // Electrons
      int m_N2plus_idx;   // Positive ions
      int m_O2plus_idx;   // Negative ions
      int m_N4plus_idx;   // Positive ions
      int m_O4plus_idx;   // Negative ions
      int m_O2plusN2_idx; // Negative ions
      int m_O2minus_idx;  // Negative ions
      int m_phot_idx;     // Photon idx

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

      Real m_photoi_factor;
  
      Real m_chemistry_dt;

      Real m_ion_mobility;  // Ion mobility
      Real m_ion_diffusion; // Ion diffusion coefficient

      int poisson_reaction(const Real a_propensity, const Real a_dt) const;

      void read_file_entries(LookupTable& a_table, const std::string a_string);
      void parse_transport_file();
      void parse_transport();
      void parse_gas_params();
      void parse_electron_mobility();
      void parse_electron_diffco();
      void parse_alpha();
      void parse_eta();
      void parse_temperature();
      void parse_see();
      void parseDomainBc();
      void parse_chemistry();
      void instantiate_species();
      void init_rng();


      Vector<Real> computeCdrFluxes(const Real         a_time,
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
				   Vector<Real>&          a_Photon_sources,
				   Vector<Real>&          a_particle_densities,
				   const Vector<RealVect> a_particle_gradients,
				   const Vector<Real>     a_Photon_densities,
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
#include "CD_NamespaceFooter.H"

#endif
