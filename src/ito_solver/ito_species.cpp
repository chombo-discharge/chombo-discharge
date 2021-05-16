/*!
  @file   ito_species.cpp
  @brief  Implementation of ito_species.H
  @author Robert Marskar
  @date   March 2020
*/

#include "ito_species.H"

#include "CD_NamespaceHeader.H"

ito_species::ito_species(){
  m_name = "ito_species.H";
}

ito_species::ito_species(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive){
  m_name      = a_name;
  m_charge    = a_charge;
  m_isMobile    = a_mobile;
  m_isDiffusive = a_diffusive;
}

ito_species::~ito_species(){

}

std::string ito_species::getName() const{
  return m_name;
}

int ito_species::get_charge() const{
  return m_charge;
}
    
bool ito_species::isDiffusive() const{
  return m_isDiffusive;
}


bool ito_species::isMobile() const {
  return m_isMobile;
}

List<ito_particle>& ito_species::get_initial_particles() {
  return m_initial_particles;
}

Real ito_species::mobility(const Real a_energy) const {
  return 0.0;
}

Real ito_species::diffusion(const Real a_energy) const {
  return 0.0;
}
#include "CD_NamespaceFooter.H"
