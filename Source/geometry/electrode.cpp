/*!
  @file electrode.cpp
  @brief Implementation of electrode.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "electrode.H"

#include "CD_NamespaceHeader.H"

electrode::electrode(){
}
  
electrode::electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction){
  this->define(a_baseif, a_live, a_fraction);
}

electrode::~electrode(){

}

void electrode::define(RefCountedPtr<BaseIF> a_baseif, bool a_live, Real a_fraction){
  m_tuple    = std::pair<RefCountedPtr<BaseIF>, bool>(a_baseif, a_live);

  if(a_live){
    m_fraction = a_fraction;
  }
  else{
    m_fraction = 0.0;
  }
}

const RefCountedPtr<BaseIF>& electrode::get_function() const {
  return m_tuple.first;
}
  
const bool& electrode::is_live() const {
  return m_tuple.second;
}

const Real& electrode::get_fraction() const {
  return m_fraction;
}
#include "CD_NamespaceFooter.H"
