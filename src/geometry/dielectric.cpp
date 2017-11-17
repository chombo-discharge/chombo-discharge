/*!
  @file dielectric.cpp
  @brief Implementation of dielectric.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "dielectric.H"

//
dielectric::dielectric(){
}
  
//
dielectric::dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  this->define(a_baseif, a_permittivity);
}

//
dielectric::~dielectric(){

}

//
void dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  m_tuple = std::pair<RefCountedPtr<BaseIF>, Real>(a_baseif, a_permittivity);
}

//
const RefCountedPtr<BaseIF>& dielectric::get_function() const {
  return m_tuple.first;
}
  
//
const Real& dielectric::get_permittivity() const {
  return m_tuple.second;
}
