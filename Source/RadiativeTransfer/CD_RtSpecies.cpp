/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSpecies.cpp
  @brief  Implementation of CD_RtSpecies.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>
  
RtSpecies::RtSpecies(){
  this->define("default_Photon", 1.0);

  m_scatter = 0.0;
  m_kappa   = 0.0;
}

RtSpecies::RtSpecies(const std::string a_name, const Real a_kappa){
  this->define(a_name, a_kappa);
}

RtSpecies::RtSpecies(const std::string a_name, Real (*a_kappa) (const RealVect a_pos)){
  this->define(a_name, a_kappa);
}

void RtSpecies::define(const std::string a_name, const Real a_kappa){
  m_name  = a_name;
  m_kappa = a_kappa;

  m_constant = true;
}

void RtSpecies::define(const std::string a_name, Real (*a_kappa)(const RealVect a_pos)){
  m_name     = a_name;
  m_varkappa = a_kappa;

  m_constant = false;
}

std::string RtSpecies::getName() const{
  return m_name;
}

bool RtSpecies::isKappaConstant() const {
  return m_constant;
}

Real RtSpecies::getKappa(const RealVect a_pos) const{
  if(m_constant){
    return m_kappa;
  }
  else{
    return m_varkappa(a_pos);
  }
}

Real RtSpecies::getScatteringCoefficient(const RealVect a_pos) const{
  if(m_constant){
    return m_scatter;
  }
  else{
    return m_varscatter(a_pos);
  }
}

#include <CD_NamespaceFooter.H>
