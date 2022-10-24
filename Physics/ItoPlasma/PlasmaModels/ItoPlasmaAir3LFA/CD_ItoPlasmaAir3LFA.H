/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaAir3LFA.H
  @brief  Ito plasma three-species model for air
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaAir3LFA_H
#define CD_ItoPlasmaAir3LFA_H

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_LookupTable1D.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ItoPlasmaAir3LFA : public ItoPlasmaPhysics
    {
    public:
      class Electron;
      class Positive;
      class Negative;
      class PhotonZ;

      ItoPlasmaAir3LFA();
      ~ItoPlasmaAir3LFA();

      Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_cdr_densities) const override;

      Real
      computeAlpha(const RealVect a_E) const override;

      Real
      computeEta(const RealVect a_E) const;

      void
      updateReactionRates(const RealVect a_E, const Real a_dx, const Real a_kappa) const noexcept override;

      Vector<Real>
      computeItoMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept override;

      Vector<Real>
      computeItoDiffusion(const Real         a_time,
                          const RealVect     a_pos,
                          const RealVect     a_E,
                          const Vector<Real> a_cdr_densities) const noexcept override;

    protected:
      // for dt
      Real m_deltaX;

      // Stuff for initial conditions
      Real     m_blob_radius;
      Real     m_particle_weight;
      int      m_num_particles;
      RealVect m_blob_center;

      int m_numPlasmaSpecies;
      int m_numRtSpecies;

      // Indices for species
      int m_ElectronIdx;
      int m_PositiveIdx;
      int m_NegativeIdx;
      int m_PhotonZ_idx;

      int m_table_entries;
      int m_init_weights;

      // Mobile ions or not
      bool m_isMobile_ions;

      // Ion mobility
      Real m_ion_mu;
      Real m_ion_D;

      // Pressure, neutral number density etc
      Real m_N;
      Real m_p;
      Real m_pq;
      Real m_T;
      Real m_O2frac;
      Real m_N2frac;
      Real m_pO2;
      Real m_photoi_factor;

      // Lookup tables if you need them.
      std::map<std::string, LookupTable1D<2>> m_tables;

      // Functions for drawing initial particles
      void
      drawInitialParticles();

      // Zheleznyak photoionization rates
      Real
      excitationRates(const Real a_E) const;
      Real
      sergeyFactor(const Real a_O2frac) const;
      Real
      photoionizationRate(const Real a_E) const;

      // Electron velocity
      RealVect
      computeElectronDriftVelocity(const RealVect a_E) const;

      /*!
	@brief Add a table
      */
      void
      addTable(const std::string a_table_name, const std::string a_file);

      // Read file into table
      void
      readTables();
    };

    class ItoPlasmaAir3LFA::Electron : public ItoSpecies
    {
    public:
      Electron();
      ~Electron();
    };

    class ItoPlasmaAir3LFA::Positive : public ItoSpecies
    {
    public:
      Positive();
      ~Positive();
    };

    class ItoPlasmaAir3LFA::Negative : public ItoSpecies
    {
    public:
      Negative();
      ~Negative();
    };

    class ItoPlasmaAir3LFA::PhotonZ : public RtSpecies
    {
    public:
      PhotonZ();
      ~PhotonZ();

      Real
      getAbsorptionCoefficient(const RealVect a_pos) const override;

    protected:
      Real m_pO2;
      Real m_f1;
      Real m_f2;
      Real m_K1;
      Real m_K2;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
