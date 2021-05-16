/*!
  @file cdr_species.cpp
  @brief Implementation of cdr_species.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "cdr_species.H"

#include "CD_NamespaceHeader.H"
  
cdr_species::cdr_species(){
  m_name         = "default_cdr_species";
  m_unit         = "default_unit";
  m_charge       = 0;
  m_isDiffusive    = true;
  m_isMobile       = true;
  m_force_output = false;

  m_init_with_function  = true;
  m_init_with_particles = true;
  
  m_deposition = DepositionType::NGP;

  m_initial_particles.clear();
}

cdr_species::cdr_species(const std::string a_name, const int a_charge, const bool a_mobile, const bool a_diffusive){
  m_name      = a_name;
  m_charge    = a_charge;
  m_isMobile    = a_mobile;
  m_isDiffusive = a_diffusive;

  m_init_with_function  = true;
  m_init_with_particles = true;
  
  m_deposition = DepositionType::NGP;
  m_initial_particles.clear();
}

cdr_species::~cdr_species(){

}

Real cdr_species::initialData(const RealVect a_pos, const Real a_time) const{
  return 0.;
}

std::string cdr_species::getName() const {
  return m_name;
}

std::string cdr_species::get_unit() const {
  return m_unit;
}

int cdr_species::get_charge() const {
  return m_charge;
}

bool cdr_species::isDiffusive() const {
  return m_isDiffusive;
}

bool cdr_species::isMobile() const {
  return m_isMobile;
}

bool cdr_species::force_output() const {
  return m_force_output;
}

bool cdr_species::init_with_particles() const {
  return m_init_with_particles;
}

bool cdr_species::init_with_function() const{
  return m_init_with_function;
}

DepositionType::Which cdr_species::get_deposition() {
  return m_deposition;
}

List<Particle>& cdr_species::get_initial_particles() {
  return m_initial_particles;
}
#include "CD_NamespaceFooter.H"
