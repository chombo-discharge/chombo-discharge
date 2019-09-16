#include "air6_mc8.H"

#ifndef _AIR6_MC8_SPECIES_
#define _AIR6_MC8_SPECIES_

class air6_mc8::electron : public species {
public:
  electron();
  ~electron() {}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::M_plus : public species {
public:
  M_plus();
  ~M_plus(){}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::M_minus : public species {
public:
  M_minus();
  ~M_minus(){}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::N2_c4v0 : public species {
public:
  N2_c4v0();
  ~N2_c4v0(){}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::N2_c4v1 : public species {
public:
  N2_c4v1();
  ~N2_c4v1(){}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::N2_b1v1 : public species {
public:
  N2_b1v1();
  ~N2_b1v1(){}
  virtual Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air6_mc8::phot_c4v0_X1v0 : public photon_group {
public:
  phot_c4v0_X1v0();
  ~phot_c4v0_X1v0(){}
  Real get_kappa(const RealVect a_pos) const {return m_kappa;}
protected:
  Real m_kappa;
};

class air6_mc8::phot_c4v0_X1v1 : public photon_group {
public:
  phot_c4v0_X1v1();
  ~phot_c4v0_X1v1(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_c4v1_X1v0 : public photon_group {
public:
  phot_c4v1_X1v0();
  ~phot_c4v1_X1v0(){};
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_c4v1_X1v1 : public photon_group {
public:
  phot_c4v1_X1v1();
  ~phot_c4v1_X1v1(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_c4v1_X1v2 : public photon_group {
public:
  phot_c4v1_X1v2();
  ~phot_c4v1_X1v2(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_c4v1_X1v3 : public photon_group {
public:
  phot_c4v1_X1v3();
  ~phot_c4v1_X1v3(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_b1v1_X1v0 : public photon_group {
public:
  phot_b1v1_X1v0();
  ~phot_b1v1_X1v0(){}
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

class air6_mc8::phot_b1v1_X1v1 : public photon_group {
public:
  phot_b1v1_X1v1();
  ~phot_b1v1_X1v1();
  Real get_kappa(const RealVect a_pos) const { return m_kappa; }
protected:
  Real m_kappa;
};

#endif
