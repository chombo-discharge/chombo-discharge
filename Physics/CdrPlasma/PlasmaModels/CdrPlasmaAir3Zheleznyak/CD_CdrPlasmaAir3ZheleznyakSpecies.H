/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3ZheleznyakSpecies.H
  @brief  Declaration of charged species and photons for the air zheleznyak model
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaAir3ZheleznyakSpecies_H
#define CD_CdrPlasmaAir3ZheleznyakSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_CdrPlasmaAir3Zheleznyak.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaAir3Zheleznyak::Electron : public CdrSpecies {
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

    class CdrPlasmaAir3Zheleznyak::MPlus : public CdrSpecies {
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

    class CdrPlasmaAir3Zheleznyak::MMinus : public CdrSpecies {
    public:
      MMinus();
      ~MMinus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class CdrPlasmaAir3Zheleznyak::uv_Photon : public RtSpecies {
    public:
      uv_Photon();
      ~uv_Photon();
      Real getKappa(const RealVect a_pos) const;
      Real getRandomKappa() const;
      Real get_pO2() const { return m_pO2;}
  
    protected:
  
      Real m_pO2;
      Real m_f1;
      Real m_f2;
      Real m_K1;
      Real m_K2;

      Real m_seed;
      std::mt19937_64* m_rng; 
      std::uniform_real_distribution<Real>* m_udist01;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
