/*!
  @file electrode.cpp
  @brief Implementation of electrode.H
  @author Robert marskar
  @date Nov. 2017
*/

#include "electrode.H"

//
electrode::electrode(){
}
  
//
electrode::electrode(RefCountedPtr<BaseIF> a_baseif, bool a_live){
  this->define(a_baseif, a_live);
}

//
electrode::~electrode(){

}

//
void electrode::define(RefCountedPtr<BaseIF> a_baseif, bool a_live){
  m_tuple = std::pair<RefCountedPtr<BaseIF>, bool>(a_baseif, a_live);
}

//
const RefCountedPtr<BaseIF>& electrode::get_function() const {
  return m_tuple.first;
}
  
//
const bool& electrode::is_live() const {
  return m_tuple.second;
}
