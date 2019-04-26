/*!
  @file   air2.H
  @brief  Declaration of species for the BOLSIG-based 8-species air model
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air2.H"

#ifndef _AIR2_SPECIES_
#define _AIR2_SPECIES_

class air2::electron : public species {
public:
  electron();
  ~electron();
  Real initial_data(const RealVect a_pos, const Real a_time) const;
protected:
  Real     m_initial_ionization;
  Real     m_seed_density;
  Real     m_seed_radius;
  RealVect m_seed_position;
};

class air2::positive_ion : public species {
public:
  positive_ion();
  ~positive_ion();
  Real initial_data(const RealVect a_pos, const Real a_time) const;
protected:

  Real     m_initial_ionization;
  Real     m_seed_density;
  Real     m_seed_radius;
  RealVect m_seed_position;
};

class air2::photon_one : public photon_group {
public:
  photon_one();
  ~photon_one();
  Real get_kappa(const RealVect a_pos) const;
  Real get_lambda() const;
  Real get_A() const;
  Real get_pO2() const;
protected:
  Real m_lambda;
  Real m_A;
  Real m_pO2;
};

class air2::photon_two : public photon_group {
public:
  photon_two();
  ~photon_two();
  Real get_kappa(const RealVect a_pos) const;
  Real get_lambda() const;
  Real get_A() const;
  Real get_pO2() const;
protected:
  Real m_lambda;
  Real m_A;
  Real m_pO2;
};

class air2::photon_three : public photon_group {
public:
  photon_three();
  ~photon_three();
  Real get_kappa(const RealVect a_pos) const;
  Real get_lambda() const;
  Real get_A() const;
  Real get_pO2() const;
protected:
  Real m_lambda;
  Real m_A;
  Real m_pO2;
};



#endif
