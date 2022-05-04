/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir7StephensSpecies.H
  @brief  Declaration of charged species and Photons for the air stephens model
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir7StephensSpecies_H
#define CD_CdrPlasmaAir7StephensSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_CdrPlasmaAir7Stephens.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaAir7Stephens::Electron : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::N2plus : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::O2plus : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::N4plus : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::O4plus : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::O2plusN2 : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::O2minus : public CdrSpecies
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

    class CdrPlasmaAir7Stephens::phot_c4v0_X1v0 : public RtSpecies
    {
    public:
      phot_c4v0_X1v0();
      ~phot_c4v0_X1v0() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_c4v0_X1v1 : public RtSpecies
    {
    public:
      phot_c4v0_X1v1();
      ~phot_c4v0_X1v1() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_c4v1_X1v0 : public RtSpecies
    {
    public:
      phot_c4v1_X1v0();
      ~phot_c4v1_X1v0(){};
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_c4v1_X1v1 : public RtSpecies
    {
    public:
      phot_c4v1_X1v1();
      ~phot_c4v1_X1v1() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_c4v1_X1v2 : public RtSpecies
    {
    public:
      phot_c4v1_X1v2();
      ~phot_c4v1_X1v2() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_c4v1_X1v3 : public RtSpecies
    {
    public:
      phot_c4v1_X1v3();
      ~phot_c4v1_X1v3() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_b1v1_X1v0 : public RtSpecies
    {
    public:
      phot_b1v1_X1v0();
      ~phot_b1v1_X1v0() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };

    class CdrPlasmaAir7Stephens::phot_b1v1_X1v1 : public RtSpecies
    {
    public:
      phot_b1v1_X1v1();
      ~phot_b1v1_X1v1() {}
      Real
      getAbsorptionCoefficient(const RealVect a_pos) const
      {
        return m_kappa;
      }
      Real
      getRandomKappa() const
      {
        return m_kappa;
      }

    protected:
      Real m_kappa;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
