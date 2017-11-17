/*!
  @file PhysicalDomain.cpp
  @brief Implementation of PhysicalDomain.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "PhysicalDomain.H"

PhysicalDomain::PhysicalDomain(){
  m_probLo = RealVect::Zero;
  m_probHi = RealVect::Zero;
}

PhysicalDomain::PhysicalDomain(const RealVect& a_probLo, const RealVect& a_probHi){
  this->define(a_probLo, a_probHi);
}

void PhysicalDomain::define(const RealVect& a_probLo, const RealVect& a_probHi){
  CH_assert(a_probLo < a_probHi);
  
  m_probLo = a_probLo;
  m_probHi = a_probHi;
}

//
const RealVect& PhysicalDomain::get_prob_lo() const{
  return m_probLo;
}

//
const RealVect& PhysicalDomain::get_prob_hi() const{
  return m_probHi;
}
