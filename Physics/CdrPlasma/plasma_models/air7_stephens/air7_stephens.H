/*!
  @file   air7_stephens.H
  @brief  7-species model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR7_STEPHENS_
#define _AIR7_STEPHENS_

#include <CD_CdrPlasmaPhysics.H>
#include "CD_LookupTable.H"

#include <random>

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace CdrPlasma {

    class air7_stephens : public CdrPlasmaPhysics {
    public:

      class electron;
      class N2plus;
      class O2plus;
      class N4plus;
      class O4plus;
      class O2plusN2;
      class O2minus;

      class phot_c4v0_X1v0;
      class phot_c4v0_X1v1;
      class phot_c4v1_X1v0;
      class phot_c4v1_X1v1;
      class phot_c4v1_X1v2;
      class phot_c4v1_X1v3;
      class phot_b1v1_X1v0;
      class phot_b1v1_X1v1;


      air7_stephens();
      ~air7_stephens();
  

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
      static std::string s_bolsig_b1_exc;
      static std::string s_bolsig_c4_exc;


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
      LookupTable m_b1_exc;        // Lookup table for excitation rates
      LookupTable m_c4_exc;        // Lookup table for excitation rates


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
      Real m_pq;     // Quenching pressure
      Real m_T;      // Temperature
      Real m_O2frac; // O2 fraction
      Real m_N2frac; // N2 fraction
      Real m_pO2;    // Partial O2 pressure

      Real m_photoi_factor;

      // Excitation efficiencies into various states
      Real m_c4v0_exc_rep;
      Real m_c4v1_exc_rep;
      Real m_b1v1_exc_rep;

      // Predissociation rates for various transitions
      Real m_c4v0_kp;
      Real m_c4v1_kp;
      Real m_b1v1_kp;

      // Radiation rates for various transitions
      Real m_c4v0_X1v0_kr;
      Real m_c4v0_X1v1_kr;
      Real m_c4v1_X1v0_kr;
      Real m_c4v1_X1v1_kr;
      Real m_c4v1_X1v2_kr;
      Real m_c4v1_X1v3_kr;
      Real m_b1v1_X1v0_kr;
      Real m_b1v1_X1v1_kr;

      // Total deactivation rates for excited states
      Real m_c4v0_k;
      Real m_c4v1_k;
      Real m_b1v1_k;

      // Total radiation rates for excited states
      Real m_c4v0_kr;
      Real m_c4v1_kr;
      Real m_b1v1_kr;

      // Photoionization efficiency for various transitions
      Real m_c4v0_X1v0_photoi_eff; // Transition c4(v=0) -> X1(v=0) + electron-ion pair
      Real m_c4v0_X1v1_photoi_eff; // Transition c4(v=0) -> X1(v=0) + electron-ion pair
      Real m_c4v1_X1v0_photoi_eff; // Transition c4(v=1) -> X1(v=0) + electron-ion pair
      Real m_c4v1_X1v1_photoi_eff; // Transition c4(v=1) -> X1(v=1) + electron-ion pair
      Real m_c4v1_X1v2_photoi_eff; // Transition c4(v=1) -> X1(v=2) + electron-ion pair
      Real m_c4v1_X1v3_photoi_eff; // Transition c4(v=1) -> X1(v=3) + electron-ion pair
      Real m_b1v1_X1v0_photoi_eff; // Transition b1(v=1) -> X1(v=0) + electron-ion pair
      Real m_b1v1_X1v1_photoi_eff; // Transition b1(v=1) -> X1(v=1) + electron-ion pair

      Real m_kq; // Quenching rate constant
  
      Real m_chemistry_dt;

      Real m_ion_mobility;  // Ion mobility
      Real m_ion_diffusion; // Ion diffusion coefficient

      int poisson_reaction(const Real a_propensity, const Real a_dt) const;
      int binomial_trials(const int a_trials, const Real a_p) const; // a_p is chance of success

      void read_file_entries(LookupTable& a_table, const std::string a_string);
      void parse_transport_file();
      void parse_transport();
      void parse_gas_params();
      void parse_electron_mobility();
      void parse_electron_diffco();
      void parse_alpha();
      void parse_eta();
      void parse_photoi();
      void parse_excitations();
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

    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
