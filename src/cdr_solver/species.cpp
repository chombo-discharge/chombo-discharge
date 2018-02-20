/*!
  @file species.cpp
  @brief Implementation of species.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "species.H"

species::species(){
  m_name      = "default_species";
  m_unit      = "default_unit";
  m_charge    = 0;
  m_diffusive = true;
}

species::species(const std::string a_name, const int a_charge, const bool a_diffusive){
  m_name      = a_name;
  m_charge    = a_charge;
  m_diffusive = a_diffusive;
}

species::~species(){

}

const Real species::initial_data(const RealVect a_pos, const Real a_time) const{
  return 0.;
}

const std::string species::get_name() const {
  return m_name;
}

const std::string species::get_unit() const {
  return m_unit;
}

const int species::get_charge() const {
  return m_charge;
}

const bool species::is_diffusive() const {
  return m_diffusive;
}
