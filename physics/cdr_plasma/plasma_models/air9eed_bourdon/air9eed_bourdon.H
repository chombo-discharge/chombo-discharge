/*!
  @file   air9eed_bourdon.H
  @brief  Declaration of a BOLSIG-based transport model for air
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR9EED_BOURDON_
#define _AIR9EED_BOURDON_

#include "cdr_plasma_physics.H"
#include "lookup_table.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {
    
    /*!
      @brief 9-species plasma model for air (79/21 N2/O2 ratio) at 300K

      @details This class uses functions from BOLSIG+ to compute ionization and transport coefficients for air at room temperature. BOLSIG+ has been pre-run for an electron energy range of \f$[10^{-3}\textrm{eV}, 20\textrm{eV}]\f$ which corresponds to electron temperatures in the range of 1K to 170000K. If you need transport data for other mixing ratios and gas temperatures, you should really, really re-run BOLSIG before running the kinetics in this program. 
    */
    class air9eed_bourdon : public cdr_plasma_physics {
    public:

      // Forward declaration of cdr species
      class eed;      // Electron energy distribution
      class electron;
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
      static void parse_gas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac);
  
      air9eed_bourdon();

      virtual ~air9eed_bourdon();

      Real compute_alpha(const RealVect a_E) const;

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

      Vector<Real> compute_cdr_source_terms(const Real              a_time,
					    const Real             a_kappa,
					    const Real             a_dx,
					    const RealVect         a_pos,
					    const RealVect         a_E,
					    const RealVect         a_gradE,
					    const Vector<Real>     a_cdr_densities,
					    const Vector<Real>     a_rte_densities,
					    const Vector<RealVect> a_grad_cdr) const;

      Vector<Real> compute_rte_source_terms(const Real         a_time,
					    const Real         a_kappa,
					    const Real         a_dx,
					    const RealVect     a_pos,
					    const RealVect     a_E,
					    const Vector<Real> a_cdr_densities) const;

      Vector<Real> compute_cdr_diffusion_coefficients(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const;
  
      Vector<RealVect> compute_cdr_velocities(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities) const;
  


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


      /*!
	@brief Set the initial surface charge
      */
      Real initial_sigma(const Real a_time, const RealVect a_pos) const; 

    protected:

      static std::string s_bolsig_energy_E;
      static std::string s_bolsig_mobility;
      static std::string s_bolsig_N2_alpha;
      static std::string s_bolsig_O2_alpha;
      static std::string s_bolsig_townsend;

      std::string m_transport_file;

      bool m_mobile_ions;
      bool m_diffusive_ions;
      bool m_diffusive_electrons;

      int m_uniform_tables;

      int m_eed_idx;
      int m_electron_idx;
      int m_N2plus_idx;
      int m_N4plus_idx;
      int m_O2plus_idx;
      int m_O4plus_idx;
      int m_O2plusN2_idx;
      int m_O2minus_idx;
      int m_Ominus_idx;
  
      int m_photon1_idx;
      int m_photon2_idx;
      int m_photon3_idx;
  
      Real m_Tg;       // Gas and ion temperature
      Real m_p;        // Gas pressure
      Real m_pq;       // Quenching pressure
      Real m_N;        // Neutral number density
      Real m_O2frac;   // O2 fraction
      Real m_N2frac;   // N2 fraction

      Real m_ion_mobility; // Ion mobility
      Real m_ion_diffusion; // Ion mobility

      Real m_cathode_work;    // Work function for cathodes
      Real m_dielectric_work; // Work function for dielectrics

      Real m_townsend2_electrode;            // 2nd Townsend on electrodes
      Real m_townsend2_dielectric;           // 2nd Townsend on dielectrics
      Real m_electrode_quantum_efficiency;   // Quantum efficiency on electrodes
      Real m_dielectric_quantum_efficiency;  // Quantum efficiency on dielectrics
      Real m_photoionization_efficiency;     // Photoionization efficiency
      Real m_excitation_efficiency;          // Excitation efficiency

      // Tables built from BOLSIG+
      lookup_table m_init_eed;
      lookup_table m_e_mobility;
      lookup_table m_e_N2_alpha;
      lookup_table m_e_O2_alpha;
      lookup_table m_e_O2_to_e_2O_c1;
      lookup_table m_e_O2_to_e_2O_c2;
      lookup_table m_e_O2_to_Ominus_O;
      lookup_table m_e_N2_to_e_N2;
      lookup_table m_e_O2_to_e_O2;
      lookup_table m_alpha_townsend;

      void parse_transport_file();
      void parse_photoi();
      void parse_see();
      void parse_transport();
      void instantiate_species();

      // Compute electron energy in eV
      inline Real electron_energy(const Real a_energy, const Real a_density);
  
      // Mobility computations
      inline Real compute_eed_mobility(const Real a_energy) const;
      inline Real compute_e_mobility(const Real a_energy) const;
      inline Real compute_N2plus_mobility(const Real a_EbyN) const;
      inline Real compute_N4plus_mobility(const Real a_EbyN) const;
      inline Real compute_O2plus_mobility(const Real a_EbyN) const;
      inline Real compute_O4plus_mobility(const Real a_EbyN) const;
      inline Real compute_O2plusN2_mobility(const Real a_EbyN) const;
      inline Real compute_O2minus_mobility(const Real a_EbyN) const;
      inline Real compute_Ominus_mobility(const Real a_EbyN) const;

      // Diffusion coefficient computatations
      inline Real compute_eed_diffco(const Real a_energy) const;
      inline Real compute_e_diffco(const Real a_energy) const;
      inline Real compute_N2plus_diffco() const;
      inline Real compute_N4plus_diffco() const;
      inline Real compute_O2plus_diffco() const;
      inline Real compute_O4plus_diffco() const;
      inline Real compute_O2plusN2_diffco() const;
      inline Real compute_O2minus_diffco() const;
      inline Real compute_Ominus_diffco() const;

      // Reactions
      inline Real compute_electron_N2_alpha(const Real a_energy) const;
      inline Real compute_electron_O2_alpha(const Real a_energy) const;
      inline Real compute_N2plus_N2_M_to_N4plus_M() const;
      inline Real compute_N4plus_O2_to_O2_2N2() const;
      inline Real compute_N2plus_O2_to_O2plus_N2(const Real a_Tg) const;
      inline Real compute_O2plus_2N2_to_O2plusN2_N2(const Real a_Tg) const;
      inline Real compute_O2plusN2_N2_to_O2plus_2N2(const Real a_Tg) const;
      inline Real compute_O2plusN2_O2_to_O4plus_N2() const;
      inline Real compute_O2plus_O2_M_to_O4plus_M(const Real a_Tg) const;
      inline Real compute_e_O4plus_to_2O2(const Real a_Te) const;
      inline Real compute_e_O2plus_to_O2(const Real a_Te) const;
      inline Real compute_e_2O2_to_O2minus_O2(const Real a_Te) const;
      inline Real compute_O2minus_O4plus_to_3O2() const;
      inline Real compute_O2minus_O4plus_M_to_3O2_M(const Real a_Tg) const;
      inline Real compute_O2minus_O2plus_M_to_2O2_M(const Real a_Tg) const;
      inline Real compute_e_O2_to_e_2O_c1(const Real a_energy) const;
      inline Real compute_e_O2_to_e_2O_c2(const Real a_energy) const;
      inline Real compute_e_O2_to_Ominus_O(const Real a_energy) const;
      inline Real compute_Oplus_O2_to_O_O2(const Real a_Tg) const;
      inline Real compute_e_N2_to_e_N2(const Real a_energy) const;
      inline Real compute_e_O2_to_e_O2(const Real a_energy) const;

      // Energy losses
      inline Real compute_e_N2_ionization_loss() const;
      inline Real compute_e_O2_ionization_loss() const;
      inline Real compute_e_O2_dissociation_loss_c1() const;
      inline Real compute_e_O2_dissociation_loss_c2() const;
      inline Real compute_e_O2_dissociative_attachment_loss() const;
      inline Real compute_e_O2_scattering_loss() const;
      inline Real compute_e_N2_scattering_loss() const;

      // How to compute the electorn energy
      inline Real compute_electron_energy(const Real a_energy_density, const Real a_electron_density) const;

      Real init_eed(const RealVect a_pos, const Real a_time, const RealVect a_E);

      void read_file_entries(lookup_table& a_table, const std::string a_string);
      void read_electron_mobility();
      void read_init_eed();
      void read_e_N2_alpha();
      void read_e_O2_alpha();
      void read_townsend();
  
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

      Vector<Real> compute_anode_fluxes(const Real         a_time,
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

      Vector<Real> compute_cathode_fluxes(const Real         a_time,
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
