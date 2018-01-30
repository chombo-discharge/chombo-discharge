/*!
  @file species.cpp
  @brief Implementation of species.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "species.H"

species::species(){
  m_name      = "default_species";
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
#if 0 
  return 0.;
#else // This is a test
  const RealVect orig = 0.0*RealVect::Unit;

  Real R = 0.1;
  Real r = (a_pos-orig).vectorLength();
  
  return exp(-r*r/(2*R*R));
#endif
}

const std::string species::get_name() const {
  return m_name;
}

const int species::get_charge() const {
  return m_charge;
}

const bool species::is_diffusive() const {
  return m_diffusive;
}
