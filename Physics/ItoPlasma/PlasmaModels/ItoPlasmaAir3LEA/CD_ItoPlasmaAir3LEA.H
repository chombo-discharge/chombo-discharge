/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaAir3LEA.H
  @brief  Declaration of a three-species air model with energy transport
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaAir3LEA_H
#define CD_ItoPlasmaAir3LEA_H

// Our includes
#include <CD_ItoPlasmaPhysicsLEA.H>
#include <CD_LookupTable1D.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ItoPlasmaAir3LEA : public ItoPlasmaPhysicsLEA
    {
    public:
      class Electron;
      class Positive;
      class Negative;
      class PhotonZ;

      ItoPlasmaAir3LEA();
      ~ItoPlasmaAir3LEA();

      Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const override;

      Real
      computeAlpha(const RealVect a_E) const override;

      Real
      computeEta(const RealVect a_E) const;

      void
      updateReactionRatesLEA(const RealVect     a_E,
                             const Vector<Real> a_mean_energies,
                             const Real         a_dx,
                             const Real         a_kappa) const override;

    protected:
      // Stuff for initial conditions
      Real     m_blob_radius;
      Real     m_particle_weight;
      int      m_num_particles;
      RealVect m_blob_center;

      // Indices for species
      int m_electronIdx;
      int m_positiveIdx;
      int m_negativeIdx;
      int m_PhotonZ_idx;

      // Pressure, neutral number density etc
      Real m_N;
      Real m_p;
      Real m_pq;
      Real m_T;
      Real m_O2frac;
      Real m_N2frac;
      Real m_pO2;
      Real m_photoi_factor;

      // Zheleznyak photoionization rates
      Real
      excitationRates(const Real a_E) const;
      Real
      sergeyFactor(const Real a_O2frac) const;
      Real
      photoionizationRate(const Real a_E) const;

      // Read file into table
      void
      readTables();
    };

    class ItoPlasmaAir3LEA::Electron : public ItoSpecies
    {
    public:
      Electron() = default;
      Electron(const LookupTable1D& a_mobility, const LookupTable1D& a_diffusion);
      ~Electron();

      Real
      mobility(const Real a_energy) const override;
      Real
      diffusion(const Real a_energy) const override;

    protected:
      LookupTable1D m_mobility;
      LookupTable1D m_diffusion;
    };

    class ItoPlasmaAir3LEA::Positive : public ItoSpecies
    {
    public:
      Positive();
      ~Positive();
    };

    class ItoPlasmaAir3LEA::Negative : public ItoSpecies
    {
    public:
      Negative();
      ~Negative();
    };

    class ItoPlasmaAir3LEA::PhotonZ : public RtSpecies
    {
    public:
      PhotonZ();
      ~PhotonZ();

      Real
      getKappa(const RealVect a_pos) const override;

    protected:
      Real m_pO2;
      Real m_f1;
      Real m_f2;
      Real m_K1;
      Real m_K2;

      Real                                  m_seed;
      std::mt19937_64*                      m_rng;
      std::uniform_real_distribution<Real>* m_udist01;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
