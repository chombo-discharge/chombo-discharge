/*!
  @file physical_domain.cpp
  @brief Implementation of physical_domain.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "physical_domain.H"

#include <ParmParse.H>

physical_domain::physical_domain(){
  m_probLo = -RealVect::Unit;
  m_probHi =  RealVect::Unit;

  ParmParse pp("physical_domain");
  
  if(pp.contains("lo_corner")){
    Vector<Real> domlo;
    domlo.resize(pp.countval("lo_corner"));
    CH_assert(domlo.size() >= SpaceDim);
    pp.getarr("lo_corner", domlo, 0, SpaceDim);
    m_probLo = RealVect(D_DECL(domlo[0], domlo[1], domlo[2]));
  }

  if(pp.contains("hi_corner")){
    Vector<Real> domhi;
    domhi.resize(pp.countval("hi_corner"));
    CH_assert(domhi.size() >= SpaceDim);
    pp.getarr("hi_corner", domhi, 0, SpaceDim);
    m_probHi = RealVect(D_DECL(domhi[0], domhi[1], domhi[2]));
  }
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
