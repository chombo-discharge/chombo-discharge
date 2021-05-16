/*!
  @file   air3_species.cpp
  @brief  Implementation of air3_species.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air3_species.H"
#include "units.H" 

#include <ParmParse.H>

air3::electron::electron(){
  m_name      = "electron density";
  m_unit      = "m-3";
  m_charge    = -1;
  m_diffusive = true;

  {// Get initial parameter
    ParmParse pp("air3");
    std::string str;
    m_initial_ionization = 1.E10;
    pp.query("initial_ionization", m_initial_ionization);
    if(pp.contains("electron_diffusion")){
      pp.get("electron_diffusion", str);
      if(str == "false"){
	m_diffusive = false;
      }
    }
  }
}

air3::electron::~electron(){

}

air3::positive_species::positive_species(){
  m_name   = "positive species density";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;

  {
    ParmParse pp("air3");
    m_initial_ionization = 1.E10;
    pp.query("initial_ionization", m_initial_ionization);
  }
}

air3::positive_species::~positive_species(){

}

air3::negative_species::negative_species(){
  m_name   = "negative species density";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = false;
}

air3::negative_species::~negative_species(){

}

air3::photon_one::photon_one(){
  m_name   = "photon_one";
  m_A      = 1.12E-4; // Default parameters
  m_lambda = 4.15E-2;

  { // Override from input script
    ParmParse pp("air3");
    pp.query("photon1_A_coeff",      m_A);
    pp.query("photon1_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air3::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air3::photon_one::~photon_one(){

}

air3::photon_two::photon_two(){
  m_name   = "photon_two";
  m_A      = 2.88E-3; // Default parameters
  m_lambda = 1.09E-1;

  { // Override from input script
    ParmParse pp("air3");
    pp.query("photon2_A_coeff",      m_A);
    pp.query("photon2_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air3::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air3::photon_two::~photon_two(){

}

air3::photon_three::photon_three(){
  m_name   = "photon_three";
  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Override from input script
    ParmParse pp("air3");
    pp.query("photon3_A_coeff",      m_A);
    pp.query("photon3_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air3::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air3::photon_three::~photon_three(){

}

Real air3::electron::initialData(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air3::positive_species::initialData(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air3::negative_species::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air3::photon_one::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air3::photon_one::get_lambda() const {
  return m_lambda;
}

Real air3::photon_one::get_A() const {
  return m_A;
}

Real air3::photon_one::get_pO2() const {
  return m_pO2;
}

Real air3::photon_two::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air3::photon_two::get_lambda() const {
  return m_lambda;
}

Real air3::photon_two::get_A() const {
  return m_A;
}

Real air3::photon_two::get_pO2() const {
  return m_pO2;
}

Real air3::photon_three::get_kappa(const RealVect a_pos) const{
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air3::photon_three::get_lambda() const {
  return m_lambda;
}

Real air3::photon_three::get_A() const {
  return m_A;
}

Real air3::photon_three::get_pO2() const {
  return m_pO2;
#include "CD_NamespaceFooter.H"

