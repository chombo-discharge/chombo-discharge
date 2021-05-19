/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_Electrode.cpp
  @brief  Implementation of CD_Electrode.H
  @author Robert marskar
*/


// Our includes
#include <CD_Electrode.H>
#include <CD_NamespaceHeader.H>

Electrode::Electrode(){
}
  
Electrode::Electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction){
  this->define(a_baseif, a_live, a_fraction);
}

Electrode::~Electrode(){

}

void Electrode::define(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction){
  m_tuple    = std::pair<RefCountedPtr<BaseIF>, bool>(a_baseif, a_live);

  if(a_live){
    m_fraction = a_fraction;
  }
  else{
    m_fraction = 0.0;
  }
}

const RefCountedPtr<BaseIF>& Electrode::getImplicitFunction() const {
  return m_tuple.first;
}
  
const bool& Electrode::isLive() const {
  return m_tuple.second;
}

const Real& Electrode::getFraction() const {
  return m_fraction;
}

#include <CD_NamespaceFooter.H>
