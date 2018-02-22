/*!
  @file dielectric.cpp
  @brief Implementation of dielectric.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "dielectric.H"

dielectric::dielectric(){
}
  
dielectric::dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  this->define(a_baseif, a_permittivity);
}

dielectric::dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos)){
  this->define(a_baseif, a_permittivity);
}

dielectric::~dielectric(){
}

void dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  m_baseif       = a_baseif;
  m_permittivity = a_permittivity;

  m_constant = true;
}

void dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos)){
  m_baseif               = a_baseif;
  m_variablepermittivity = a_permittivity;
  
  m_constant = false;
}

const RefCountedPtr<BaseIF>& dielectric::get_function() const {
  return m_baseif;
}
  
Real dielectric::get_permittivity(const RealVect a_pos) const {
  if(m_constant){
    return m_permittivity;
  }
  else{
    return m_variablepermittivity(a_pos);
  }
}
