/*!
  @file   air3_zheleznyak_species.H
  @brief  Declaration of charged species and photons for the air zheleznyak model
  @author Robert Marskar
*/

#ifndef _AIR3_ZHELEZNYAK_SPECIES_
#define _AIR3_ZHELEZNYAK_SPECIES_

#include "air3_zheleznyak.H"
#include <random>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {

    class air3_zheleznyak::electron : public cdr_species {
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

    class air3_zheleznyak::M_plus : public cdr_species {
    public:
      M_plus();
      ~M_plus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air3_zheleznyak::M_minus : public cdr_species {
    public:
      M_minus();
      ~M_minus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class air3_zheleznyak::uv_photon : public rte_species {
    public:
      uv_photon();
      ~uv_photon();
      Real get_kappa(const RealVect a_pos) const;
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
