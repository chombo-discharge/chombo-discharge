#ifndef _AIR3_BOURDON_SPECIES_
#define _AIR3_BOURDON_SPECIES_

#include "air3_bourdon.H"
#include <random>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {

    class air3_bourdon::electron : public CdrSpecies {
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

    class air3_bourdon::M_plus : public CdrSpecies {
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

    class air3_bourdon::M_minus : public CdrSpecies {
    public:
      M_minus();
      ~M_minus(){}
      Real initialData(const RealVect a_pos, const Real a_time) const {return 0.0;}
    };

    class air3_bourdon::Photon_one : public RtSpecies {
    public:
      Photon_one();
      ~Photon_one(){}
      Real getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
      Real get_lambda() const { return m_lambda;}
      Real get_A()      const { return m_A;}
      Real get_pO2()    const { return m_pO2; }

    protected:
      Real m_lambda;
      Real m_A;
      Real m_pO2;
    };

    class air3_bourdon::Photon_two : public RtSpecies {
    public:
      Photon_two();
      ~Photon_two(){}
      Real getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
      Real get_lambda() const { return m_lambda;}
      Real get_A()      const { return m_A;}
      Real get_pO2()    const { return m_pO2; }

    protected:
      Real m_lambda;
      Real m_A;
      Real m_pO2;
    };

    class air3_bourdon::Photon_three : public RtSpecies {
    public:
      Photon_three();
      ~Photon_three(){}
      Real getKappa(const RealVect a_pos) const { return m_lambda*m_pO2/sqrt(3.0);}
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
#include "CD_NamespaceFooter.H"
#endif
