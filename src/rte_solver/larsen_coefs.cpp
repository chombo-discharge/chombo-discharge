/*!
  @file   larsen_coefs.cpp
  @brief  Implementation of larsen_coefs.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "larsen_coefs.H"

larsen_coefs::larsen_coefs(const RefCountedPtr<rte_species>& a_rte_species,
			   const Real                         a_r1,
			   const Real                         a_r2){
  m_rte_species = a_rte_species;
  m_r1 = a_r1;
  m_r2 = a_r2;
}

larsen_coefs::~larsen_coefs(){

}

Real larsen_coefs::aco(const RealVect a_pos) const {

  Real val = m_rte_species->get_kappa(a_pos)*m_rte_species->get_kappa(a_pos);
  val *= 3.0/2.0;
  val *= (1 + 3*m_r2)/(1 - 2*m_r1);

  return val;
}

Real larsen_coefs::bco(const RealVect a_pos) const {
  return -m_rte_species->get_kappa(a_pos);
}

Real larsen_coefs::rhs(const RealVect a_pos) const {
  return 0.;
}
