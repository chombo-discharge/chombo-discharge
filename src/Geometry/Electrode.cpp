/*!
  @file Electrode.cpp
  @brief Implementation of Electrode.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "Electrode.H"

//
Electrode::Electrode(){
}
  
//
Electrode::Electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live){
  this->define(a_baseif, a_live);
}

//
Electrode::~Electrode(){

}

//
void Electrode::define(RefCountedPtr<BaseIF> a_baseif, bool a_live){
  m_tuple = std::pair<RefCountedPtr<BaseIF>, bool>(a_baseif, a_live);
}

//
const RefCountedPtr<BaseIF>& Electrode::get_function() const {
  return m_tuple.first;
}
  
//
const bool& Electrode::is_live() const {
  return m_tuple.second;
}
