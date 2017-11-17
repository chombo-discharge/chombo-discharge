/*!
  @file Dielectric.cpp
  @brief Implementation of Dielectric.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "Dielectric.H"

//
Dielectric::Dielectric(){
}
  
//
Dielectric::Dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  this->define(a_baseif, a_permittivity);
}

//
Dielectric::~Dielectric(){

}

//
void Dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  m_tuple = std::pair<RefCountedPtr<BaseIF>, Real>(a_baseif, a_permittivity);
}

//
const RefCountedPtr<BaseIF>& Dielectric::get_function() const {
  return m_tuple.first;
}
  
//
const Real& Dielectric::get_permittivity() const {
  return m_tuple.second;
}
