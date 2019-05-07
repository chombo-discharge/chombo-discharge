/*!
  @file   photon_group.cpp
  @brief  Implementation of photon_group.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "photon_group.H"

photon_group::photon_group(){
  this->define("default_photon", 1.0);

  m_scatter = 0.0;
  m_kappa = 0.0;
}

photon_group::photon_group(const std::string a_name, const Real a_kappa){
  this->define(a_name, a_kappa);
}

photon_group::photon_group(const std::string a_name, Real (*a_kappa) (const RealVect a_pos)){
  this->define(a_name, a_kappa);
}

void photon_group::define(const std::string a_name, const Real a_kappa){
  m_name  = a_name;
  m_kappa = a_kappa;

  m_constant = true;
}

void photon_group::define(const std::string a_name, Real (*a_kappa)(const RealVect a_pos)){
  m_name     = a_name;
  m_varkappa = a_kappa;

  m_constant = false;
}

std::string photon_group::get_name() const{
  return m_name;
}

Real photon_group::get_kappa(const RealVect a_pos) const{
  if(m_constant){
    return m_kappa;
  }
  else{
    return m_varkappa(a_pos);
  }
}

Real photon_group::get_scatter(const RealVect a_pos) const{
  if(m_constant){
    return m_scatter;
  }
  else{
    return m_varscatter(a_pos);
  }
}
