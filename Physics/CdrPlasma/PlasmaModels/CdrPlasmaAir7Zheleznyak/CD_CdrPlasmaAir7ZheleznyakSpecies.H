/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir7ZheleznyakSpecies.H
  @brief  Declaration of charged species and Photons for the air zheleznyak model
  @author Robert Marskar
*/

#ifndef CdrPlasmaAir7ZheleznyakSpecies_H
#define CdrPlasmaAir7ZheleznyakSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_CdrPlasmaAir7Zheleznyak.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaAir7Zheleznyak::Electron : public CdrSpecies
    {
    public:
      Electron();
      ~Electron() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real     m_uniform_density;
      Real     m_seed_density;
      Real     m_seed_rad;
      RealVect m_seed_pos;
    };

    class CdrPlasmaAir7Zheleznyak::N2plus : public CdrSpecies
    {
    public:
      N2plus();
      ~N2plus() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real     m_frac;
      Real     m_uniform_density;
      Real     m_seed_density;
      Real     m_seed_rad;
      RealVect m_seed_pos;
    };

    class CdrPlasmaAir7Zheleznyak::O2plus : public CdrSpecies
    {
    public:
      O2plus();
      ~O2plus() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real     m_frac;
      Real     m_uniform_density;
      Real     m_seed_density;
      Real     m_seed_rad;
      RealVect m_seed_pos;
    };

    class CdrPlasmaAir7Zheleznyak::N4plus : public CdrSpecies
    {
    public:
      N4plus();
      ~N4plus() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const
      {
        return 0.0;
      }

    protected:
    };

    class CdrPlasmaAir7Zheleznyak::O4plus : public CdrSpecies
    {
    public:
      O4plus();
      ~O4plus() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const
      {
        return 0.0;
      }

    protected:
    };

    class CdrPlasmaAir7Zheleznyak::O2plusN2 : public CdrSpecies
    {
    public:
      O2plusN2();
      ~O2plusN2() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const
      {
        return 0.0;
      }

    protected:
    };

    class CdrPlasmaAir7Zheleznyak::O2minus : public CdrSpecies
    {
    public:
      O2minus();
      ~O2minus() {}
      Real
      initialData(const RealVect a_pos, const Real a_time) const
      {
        return 0.0;
      }
    };

    class CdrPlasmaAir7Zheleznyak::uv_Photon : public RtSpecies
    {
    public:
      uv_Photon();
      ~uv_Photon();
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const;
      Real
      getRandomKappa() const;
      Real
      getPO2() const
      {
        return m_pO2;
      }

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
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
