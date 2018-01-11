/*!
  @file physical_domain.cpp
  @brief Implementation of physical_domain.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "physical_domain.H"

physical_domain::physical_domain(){
  m_probLo = RealVect::Zero;
  m_probHi = RealVect::Zero;
}

physical_domain::physical_domain(const RealVect& a_probLo, const RealVect& a_probHi){
  this->define(a_probLo, a_probHi);
}

void physical_domain::define(const RealVect& a_probLo, const RealVect& a_probHi){
  CH_assert(a_probLo < a_probHi);
  
  m_probLo = a_probLo;
  m_probHi = a_probHi;
}

//
const RealVect& physical_domain::get_prob_lo() const{
  return m_probLo;
}

//
const RealVect& physical_domain::get_prob_hi() const{
  return m_probHi;
}
