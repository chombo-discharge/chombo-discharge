/*!
  @file   rte_species.cpp
  @brief  Implementation of rte_species.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rte_species.H"

#include "CD_NamespaceHeader.H"
  
rte_species::rte_species(){
  this->define("default_photon", 1.0);

  m_scatter = 0.0;
  m_kappa = 0.0;
}

rte_species::rte_species(const std::string a_name, const Real a_kappa){
  this->define(a_name, a_kappa);
}

rte_species::rte_species(const std::string a_name, Real (*a_kappa) (const RealVect a_pos)){
  this->define(a_name, a_kappa);
}

void rte_species::define(const std::string a_name, const Real a_kappa){
  m_name  = a_name;
  m_kappa = a_kappa;

  m_constant = true;
}

void rte_species::define(const std::string a_name, Real (*a_kappa)(const RealVect a_pos)){
  m_name     = a_name;
  m_varkappa = a_kappa;

  m_constant = false;
}

std::string rte_species::getName() const{
  return m_name;
}

bool rte_species::constant_kappa() const {
  return m_constant;
}

Real rte_species::getKappa(const RealVect a_pos) const{
  if(m_constant){
    return m_kappa;
  }
  else{
    return m_varkappa(a_pos);
  }
}

Real rte_species::get_scatter(const RealVect a_pos) const{
  if(m_constant){
    return m_scatter;
  }
  else{
    return m_varscatter(a_pos);
  }
}
#include "CD_NamespaceFooter.H"
