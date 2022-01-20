/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir9EedBourdon.H
  @brief  Declaration of a BOLSIG-based transport model for air with electron energy transport
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir9EedBourdon_H
#define CD_CdrPlasmaAir9EedBourdon_H

#include <CD_CdrPlasmaPhysics.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
    
    /*!
      @brief 9-species plasma model for air (79/21 N2/O2 ratio) at 300K

      @details This class uses functions from BOLSIG+ to compute ionization and transport coefficients for air at room temperature. BOLSIG+ has been pre-run for an Electron energy range of \f$[10^{-3}\textrm{eV}, 20\textrm{eV}]\f$ which corresponds to Electron temperatures in the range of 1K to 170000K. If you need transport data for other mixing ratios and gas temperatures, you should really, really re-run BOLSIG before running the kinetics in this program. 
    */
    class CdrPlasmaAir9EedBourdon : public CdrPlasmaPhysics {
    public:

      // Forward declaration of cdr species
      class eed;      // Electron energy distribution
      class Electron;
      class N2plus;   
      class N4plus;   
      class O2plus;
      class O4plus;
      class O2plusN2;
      class O2minus;
      class Ominus;
      class O;

      // Forward declaration of Photon excitations
      class PhotonOne;
      class PhotonTwo;
      class PhotonThree;

      /*!
	@brief Get gas temperature, pressure, density, and molecular fractions
      */
      static void parseGas_parameters(Real& a_Tg, Real& a_p, Real& a_N, Real& a_O2frac, Real& a_N2frac);
  
      CdrPlasmaAir9EedBourdon();

      virtual ~CdrPlasmaAir9EedBourdon();

      Real computeAlpha(const Real a_E, const RealVect a_pos) const;

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

      Vector<Real> computeCdrSourceTerms(const Real              a_time,
					 const Real             a_kappa,
					 const Real             a_dx,
					 const RealVect         a_pos,
					 const RealVect         a_E,
					 const RealVect         a_gradE,
					 const Vector<Real>     a_cdr_densities,
					 const Vector<Real>     a_rte_densities,
					 const Vector<RealVect> a_grad_cdr) const;

      Vector<Real> computeRtSourceTerms(const Real         a_time,
					const Real         a_kappa,
					const Real         a_dx,
					const RealVect     a_pos,
					const RealVect     a_E,
					const Vector<Real> a_cdr_densities) const;

      Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
						   const RealVect     a_pos,
						   const RealVect     a_E,
						   const Vector<Real> a_cdr_densities) const;
  
      Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities) const;
  


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


      /*!
	@brief Set the initial surface charge
      */
      Real initialSigma(const Real a_time, const RealVect a_pos) const; 

    protected:

      static std::string s_bolsig_energy_E;
      static std::string s_bolsig_mobility;
      static std::string s_bolsig_N2_alpha;
      static std::string s_bolsig_O2_alpha;
      static std::string s_bolsig_townsend;

      std::string m_transport_file;

      bool m_isMobile_ions;
      bool m_isDiffusive_ions;
      bool m_isDiffusive_Electrons;

      int m_uniform_tables;

      int m_eed_idx;
      int m_Electron_idx;
      int m_N2plus_idx;
      int m_N4plus_idx;
      int m_O2plus_idx;
      int m_O4plus_idx;
      int m_O2plusN2_idx;
      int m_O2minus_idx;
      int m_Ominus_idx;
  
      int m_Photon1_idx;
      int m_Photon2_idx;
      int m_Photon3_idx;
  
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
      LookupTable<2> m_initEed;
      LookupTable<2> m_e_mobility;
      LookupTable<2> m_e_N2_alpha;
      LookupTable<2> m_e_O2_alpha;
      LookupTable<2> m_e_O2_to_e_2O_c1;
      LookupTable<2> m_e_O2_to_e_2O_c2;
      LookupTable<2> m_e_O2_to_Ominus_O;
      LookupTable<2> m_e_N2_to_e_N2;
      LookupTable<2> m_e_O2_to_e_O2;
      LookupTable<2> m_alpha_townsend;

      void parseTransportFile();
      void parsePhotoi();
      void parseSEE();
      void parseTransport();
      void initSpecies();

      // Compute Electron energy in eV
      inline Real Electron_energy(const Real a_energy, const Real a_density);
  
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
      inline Real compute_Electron_N2_alpha(const Real a_energy) const;
      inline Real compute_Electron_O2_alpha(const Real a_energy) const;
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
      inline Real compute_Electron_energy(const Real a_energy_density, const Real a_Electron_density) const;

      Real initEed(const RealVect a_pos, const Real a_time, const RealVect a_E);

      void read_Electron_mobility();
      void read_initEed();
      void read_e_N2_alpha();
      void read_e_O2_alpha();
      void read_townsend();
  
      /*!
	@brief Common function for dielectric and electrode fluxes, only difference is the secondary emission coefficients
      */
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

      Vector<Real> computeAnodeFluxes(const Real         a_time,
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

      Vector<Real> computeCathodeFluxes(const Real         a_time,
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

#include <CD_NamespaceFooter.H>

#endif
