/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3BourdonSpecies.H
  @brief  Declaration of species for CdrPlasmaAir3Bourdon model
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir3BourdonSpecies_H
#define CD_CdrPlasmaAir3BourdonSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_CdrPlasmaAir3Bourdon.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaAir3Bourdon::Electron : public CdrSpecies {
    public:
      Electron();
      ~Electron() {}
      Real initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class CdrPlasmaAir3Bourdon::MPlus : public CdrSpecies {
    public:
      MPlus();
      ~MPlus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class CdrPlasmaAir3Bourdon::MMinus : public CdrSpecies {
    public:
      MMinus();
      ~MMinus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class CdrPlasmaAir3Bourdon::PhotonOne : public RtSpecies {
    public:
      PhotonOne();
      ~PhotonOne(){}
      Real getAbsorptionCoefficient(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
      Real get_lambda() const { return m_lambda;}
      Real get_A()      const { return m_A;}
      Real get_pO2()    const { return m_pO2; }

    protected:
      Real m_lambda;
      Real m_A;
      Real m_pO2;
    };

    class CdrPlasmaAir3Bourdon::PhotonTwo : public RtSpecies {
    public:
      PhotonTwo();
      ~PhotonTwo(){}
      Real getAbsorptionCoefficient(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
      Real get_lambda() const { return m_lambda;}
      Real get_A()      const { return m_A;}
      Real get_pO2()    const { return m_pO2; }

    protected:
      Real m_lambda;
      Real m_A;
      Real m_pO2;
    };

    class CdrPlasmaAir3Bourdon::PhotonThree : public RtSpecies {
    public:
      PhotonThree();
      ~PhotonThree(){}
      Real getAbsorptionCoefficient(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
      Real get_lambda() const { return m_lambda;}
      Real get_A()      const { return m_A;}
      Real get_pO2()    const { return m_pO2; }

    protected:
      Real m_lambda;
      Real m_A;
      Real m_pO2;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
