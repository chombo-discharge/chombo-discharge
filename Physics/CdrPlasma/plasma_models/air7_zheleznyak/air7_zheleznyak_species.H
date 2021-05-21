/*!
  @file   air7_zheleznyak_species.H
  @brief  Declaration of charged species and Photons for the air zheleznyak model
  @author Robert Marskar
*/

#ifndef _AIR7_ZHELEZNYAK_SPECIES_
#define _AIR7_ZHELEZNYAK_SPECIES_

#include "air7_zheleznyak.H"
#include <random>

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace CdrPlasma {
  
    class air7_zheleznyak::electron : public CdrSpecies {
    public:
      electron();
      ~electron() {}
      Real initialData(const RealVect a_pos, const Real a_time) const;

    protected:
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::N2plus : public CdrSpecies {
    public:
      N2plus();
      ~N2plus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_frac;
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::O2plus : public CdrSpecies {
    public:
      O2plus();
      ~O2plus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_frac;
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::N4plus : public CdrSpecies {
    public:
      N4plus();
      ~N4plus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const{return 0.0;}
    protected:
    };

    class air7_zheleznyak::O4plus : public CdrSpecies {
    public:
      O4plus();
      ~O4plus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const{return 0.0;}
    protected:
    };

    class air7_zheleznyak::O2plusN2 : public CdrSpecies {
    public:
      O2plusN2();
      ~O2plusN2(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    protected:
    };

    class air7_zheleznyak::O2minus : public CdrSpecies {
    public:
      O2minus();
      ~O2minus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class air7_zheleznyak::uv_Photon : public RtSpecies {
    public:
      uv_Photon();
      ~uv_Photon();
      Real getKappa(const RealVect a_pos) const;
      Real get_random_kappa() const;
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
#include "CD_NamespaceFooter.H"
#endif
