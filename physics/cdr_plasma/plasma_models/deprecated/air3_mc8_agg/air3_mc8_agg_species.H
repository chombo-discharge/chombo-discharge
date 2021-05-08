#ifndef _AIR3_MC8_AGG_SPECIES_
#define _AIR3_MC8_AGG_SPECIES_

#include "air3_mc8_agg.H"

#include <chrono>
#include <random>

class air3_mc8_agg::electron : public species {
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

class air3_mc8_agg::M_plus : public species {
public:
  M_plus();
  ~M_plus(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const;
protected:
  Real m_uniform_density;
  Real m_seed_density;
  Real m_seed_rad;
  RealVect m_seed_pos;
};

class air3_mc8_agg::M_minus : public species {
public:
  M_minus();
  ~M_minus(){}
  Real initial_data(const RealVect a_pos, const Real a_time) const {return 0.0;}
};

class air3_mc8_agg::agg_photon : public photon_group {
public:
  agg_photon();
  ~agg_photon(){}
  Real get_kappa(const RealVect a_pos) const {return m_kappa;}
  Real get_random_kappa() const;
protected:
  Real m_kappa;

  // RNG stuff
  int m_rng_seed;
  std::mt19937_64*                 m_rng;
  std::uniform_real_distribution<Real>* m_udist;

  // Excitation efficiencies for various transitions
  Real m_c4v0_exc_eff;
  Real m_c4v1_exc_eff;
  Real m_b1v1_exc_eff;

  // Photoionization efficiency for various transitions
  Real m_c4v0_X1v0_photoi_eff; // Transition c4(v=0) -> X1(v=0) + electron-ion pair
  Real m_c4v0_X1v1_photoi_eff; // Transition c4(v=0) -> X1(v=0) + electron-ion pair
  Real m_c4v1_X1v0_photoi_eff; // Transition c4(v=1) -> X1(v=0) + electron-ion pair
  Real m_c4v1_X1v1_photoi_eff; // Transition c4(v=1) -> X1(v=1) + electron-ion pair
  Real m_c4v1_X1v2_photoi_eff; // Transition c4(v=1) -> X1(v=2) + electron-ion pair
  Real m_c4v1_X1v3_photoi_eff; // Transition c4(v=1) -> X1(v=3) + electron-ion pair
  Real m_b1v1_X1v0_photoi_eff; // Transition b1(v=1) -> X1(v=0) + electron-ion pair
  Real m_b1v1_X1v1_photoi_eff; // Transition b1(v=1) -> X1(v=1) + electron-ion pair

  Real m_prob_c4v0_X1v0; // Relative probability for the photon being a c4v0->X1v0 photon
  Real m_prob_c4v0_X1v1; // Relative probability for the photon being a c4v0->X1v1 photon
  Real m_prob_c4v1_X1v0; // Relative probability for the photon being a c4v1->X1v0 photon
  Real m_prob_c4v1_X1v1; // Relative probability for the photon being a c4v1->X1v1 photon
  Real m_prob_c4v1_X1v2; // Relative probability for the photon being a c4v1->X1v2 photon
  Real m_prob_c4v1_X1v3; // Relative probability for the photon being a c4v1->X1v3 photon
  Real m_prob_b1v1_X1v0; // Relative probability for the photon being a b1v1->X1v0 photon
  Real m_prob_b1v1_X1v1; // Relative probability for the photon being a b1v1->X1v1 photon

  // Probabilities and kappas
  Vector<Real> m_probs;
  Vector<Real> m_kappas;

  int draw_photon_type() const;
#include "CD_NamespaceFooter.H"

#endif
