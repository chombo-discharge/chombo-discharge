/*!
  @file   air7_zheleznyak_species.H
  @brief  Declaration of charged species and photons for the air zheleznyak model
  @author Robert Marskar
*/

#ifndef _AIR7_ZHELEZNYAK_SPECIES_
#define _AIR7_ZHELEZNYAK_SPECIES_

#include "air7_zheleznyak.H"
#include <random>

namespace physics {
  namespace cdr_plasma {
  
    class air7_zheleznyak::electron : public cdr_species {
    public:
      electron();
      ~electron() {}
      Real initial_data(const RealVect a_pos, const Real a_time) const;

    protected:
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::N2plus : public cdr_species {
    public:
      N2plus();
      ~N2plus(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_frac;
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::O2plus : public cdr_species {
    public:
      O2plus();
      ~O2plus(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const;
    protected:
      Real m_frac;
      Real m_uniform_density;
      Real m_seed_density;
      Real m_seed_rad;
      RealVect m_seed_pos;
    };

    class air7_zheleznyak::N4plus : public cdr_species {
    public:
      N4plus();
      ~N4plus(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const{return 0.0;}
    protected:
    };

    class air7_zheleznyak::O4plus : public cdr_species {
    public:
      O4plus();
      ~O4plus(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const{return 0.0;}
    protected:
    };

    class air7_zheleznyak::O2plusN2 : public cdr_species {
    public:
      O2plusN2();
      ~O2plusN2(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
    protected:
    };

    class air7_zheleznyak::O2minus : public cdr_species {
    public:
      O2minus();
      ~O2minus(){}
      Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class air7_zheleznyak::uv_photon : public rte_species {
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
#endif
