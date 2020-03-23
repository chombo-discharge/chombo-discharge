/*!
  @file   air7_stephens_species.H
  @brief  Declaration of charged species and photons for the air stephens model
  @author Robert Marskar
*/

#ifndef _AIR7_STEPHENS_SPECIES_
#define _AIR7_STEPHENS_SPECIES_

#include "air7_stephens.H"
#include <random>

class air7_stephens::electron : public cdr_species {
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

class air7_stephens::N2plus : public cdr_species {
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

class air7_stephens::O2plus : public cdr_species {
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

class air7_stephens::N4plus : public cdr_species {
public:
  N4plus();
  ~N4plus(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const{return 0.0;}
protected:
};

class air7_stephens::O4plus : public cdr_species {
public:
  O4plus();
  ~O4plus(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const{return 0.0;}
protected:
};

class air7_stephens::O2plusN2 : public cdr_species {
public:
  O2plusN2();
  ~O2plusN2(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
protected:
};

class air7_stephens::O2minus : public cdr_species {
public:
  O2minus();
  ~O2minus(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air7_stephens::phot_c4v0_X1v0 : public rte_species {
public:
  phot_c4v0_X1v0();
  ~phot_c4v0_X1v0(){}
  Real get_kappa(const RealVect a_pos) const {return m_kappa;}
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_c4v0_X1v1 : public rte_species {
public:
  phot_c4v0_X1v1();
  ~phot_c4v0_X1v1(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_c4v1_X1v0 : public rte_species {
public:
  phot_c4v1_X1v0();
  ~phot_c4v1_X1v0(){};
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_c4v1_X1v1 : public rte_species {
public:
  phot_c4v1_X1v1();
  ~phot_c4v1_X1v1(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_c4v1_X1v2 : public rte_species {
public:
  phot_c4v1_X1v2();
  ~phot_c4v1_X1v2(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_c4v1_X1v3 : public rte_species {
public:
  phot_c4v1_X1v3();
  ~phot_c4v1_X1v3(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_b1v1_X1v0 : public rte_species {
public:
  phot_b1v1_X1v0();
  ~phot_b1v1_X1v0(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

class air7_stephens::phot_b1v1_X1v1 : public rte_species {
public:
  phot_b1v1_X1v1();
  ~phot_b1v1_X1v1();
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
  Real get_random_kappa() const {return m_kappa;}
protected:
  Real m_kappa;
};

#endif
