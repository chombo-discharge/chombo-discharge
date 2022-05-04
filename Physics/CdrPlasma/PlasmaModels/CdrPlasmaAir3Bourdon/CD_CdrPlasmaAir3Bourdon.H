/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3Bourdon.H
  @brief  Standard three-species model for air discharges. 
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir3Bourdon_H
#define CD_CdrPlasmaAir3Bourdon_H

// Std includes
#include <random>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaAir3Bourdon : public CdrPlasmaPhysics
    {
    public:
      class Electron;
      class MPlus;
      class MMinus;

      class PhotonOne;
      class PhotonTwo;
      class PhotonThree;

      CdrPlasmaAir3Bourdon();
      ~CdrPlasmaAir3Bourdon();

      void
      advanceReactionNetwork(Vector<Real>&          a_particle_sources,
                             Vector<Real>&          a_Photon_sources,
                             const Vector<Real>     a_particle_densities,
                             const Vector<RealVect> a_particle_gradients,
                             const Vector<Real>     a_Photon_densities,
                             const RealVect         a_E,
                             const RealVect         a_pos,
                             const Real             a_dx,
                             const Real             a_dt,
                             const Real             a_time,
                             const Real             a_kappa) const override;

      Vector<Real>
      computeCdrDiffusionCoefficients(const Real         a_time,
                                      const RealVect     a_pos,
                                      const RealVect     a_E,
                                      const Vector<Real> a_cdr_densities) const override;

      Vector<RealVect>
      computeCdrDriftVelocities(const Real         a_time,
                                const RealVect     a_pos,
                                const RealVect     a_E,
                                const Vector<Real> a_cdr_densities) const override;

      Vector<Real>
      computeCdrDomainFluxes(const Real           a_time,
                             const RealVect       a_pos,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const RealVect       a_E,
                             const Vector<Real>   a_cdr_densities,
                             const Vector<Real>   a_cdr_velocities,
                             const Vector<Real>   a_cdr_gradients,
                             const Vector<Real>   a_rte_fluxes,
                             const Vector<Real>   a_extrap_cdr_fluxes) const override;

      Vector<Real>
      computeCdrElectrodeFluxes(const Real         a_time,
                                const RealVect     a_pos,
                                const RealVect     a_normal,
                                const RealVect     a_E,
                                const Vector<Real> a_cdr_densities,
                                const Vector<Real> a_cdr_velocities,
                                const Vector<Real> a_cdr_gradients,
                                const Vector<Real> a_rte_fluxes,
                                const Vector<Real> a_extrap_cdr_fluxes) const override;

      Vector<Real>
      computeCdrDielectricFluxes(const Real         a_time,
                                 const RealVect     a_pos,
                                 const RealVect     a_normal,
                                 const RealVect     a_E,
                                 const Vector<Real> a_cdr_densities,
                                 const Vector<Real> a_cdr_velocities,
                                 const Vector<Real> a_cdr_gradients,
                                 const Vector<Real> a_rte_fluxes,
                                 const Vector<Real> a_extrap_cdr_fluxes) const override;

      Real
      initialSigma(const Real a_time, const RealVect a_pos) const override;

      Real
      computeAlpha(const Real a_E, const RealVect a_pos) const override;

    protected:
      enum
      {
        ALPHA  = 1,
        ETA    = 1,
        DIFFCO = 1,
        MU     = 1
      };

      static std::string s_bolsig_mobility;
      static std::string s_bolsig_diffco;
      static std::string s_bolsig_alpha;
      static std::string s_bolsig_eta;
      static std::string s_bolsig_energy;

      std::string m_transport_file;

      // Tables built from BOLSIG+
      LookupTable<2> m_e_mobility;
      LookupTable<2> m_e_diffco;
      LookupTable<2> m_e_alpha;
      LookupTable<2> m_e_eta;
      LookupTable<2> m_e_energy;

      Vector<int> m_wallBc; // Wall boundary conditions

      bool m_alpha_corr;            // Alpha correction
      bool m_isMobile_Electrons;    // Mobile Electrons or not
      bool m_isMobile_ions;         // Mobile ions or not
      bool m_isDiffusive_ions;      // Diffusive ions or not
      bool m_isDiffusive_Electrons; // Diffusive Electrons or not

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
      Real m_photoi_eff;   // Photoionizationn efficiency

      void
      parseDomainBc();
      void
      initSpecies();

      Vector<Real>
      computeCdrFluxes(const Real         a_time,
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
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
