/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_Dielectric.cpp
  @brief  Implementation of CD_Dielectric.H
  @author Robert marskar
  @date Nov. 2017
*/

// Our includes
#include <CD_Dielectric.H>
#include <CD_NamespaceHeader.H>
  
Dielectric::Dielectric(){
}
  
Dielectric::Dielectric(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  this->define(a_baseif, a_permittivity);
}

Dielectric::Dielectric(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos)){
  this->define(a_baseif, a_permittivity);
}

Dielectric::~Dielectric(){
}

void Dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real a_permittivity){
  m_baseif       = a_baseif;
  m_permittivity = a_permittivity;

  m_constant = true;
}

void Dielectric::define(RefCountedPtr<BaseIF> a_baseif, Real (*a_permittivity)(const RealVect a_pos)){
  m_baseif               = a_baseif;
  m_variablepermittivity = a_permittivity;
  
  m_constant = false;
}

const RefCountedPtr<BaseIF>& Dielectric::getImplicitFunction() const {
  return m_baseif;
}
  
Real Dielectric::getPermittivity(const RealVect a_pos) const {
  if(m_constant){
    return m_permittivity;
  }
  else{
    return m_variablepermittivity(a_pos);
  }
}

#include <CD_NamespaceFooter.H>
