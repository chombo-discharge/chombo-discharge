/*!
  @file   air9eed_species.cpp
  @brief  Implementation of air9eed_species.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Could definitely see ways to cut down on the typic 
*/

#include "air9eed_species.H"
#include "units.H" 

#include <ParmParse.H>

air9eed::eed::eed(){
  m_name      = "electron energy density";
  m_unit      = "eVm-3";
  m_charge    = 0;
  m_diffusive = false;
  m_mobile    = true;

  // Get gas parameters
  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
}

air9eed::eed::~eed(){

}

air9eed::electron::electron(){
  m_name   = "electron density";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = false;
  m_mobile = true;

  {// Get initial parameter
    ParmParse pp("air9eed");
    pp.get("initial_ionization", m_initial_ionization);
  }
}

air9eed::electron::~electron(){

}

air9eed::N2plus::N2plus() {
  m_name   = "N2plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  ParmParse pp("air9eed");
  pp.get("initial_ionization", m_initial_ionization);
  m_initial_ionization *= N2frac;

  std::string str;
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::N2plus::~N2plus(){

}

air9eed::N4plus::N4plus(){
  m_name   = "N4plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

    std::string str;
  ParmParse pp("air9eed");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::N4plus::~N4plus(){

}

air9eed::O2plus::O2plus(){
  m_name   = "O2plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  std::string str;
  ParmParse pp("air9eed");
  pp.get("initial_ionization", m_initial_ionization);
  m_initial_ionization *= O2frac;

  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::O2plus::~O2plus(){

}

air9eed::O4plus::O4plus(){
  m_name   = "O4plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air9eed");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::O4plus::~O4plus(){

}

air9eed::O2plusN2::O2plusN2() {
  m_name   = "O2plusN2";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air9eed");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::O2plusN2::~O2plusN2(){

}

air9eed::O2minus::O2minus(){
  m_name   = "O2minus";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = false;
  m_mobile = true;

    std::string str;
  ParmParse pp("air9eed");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::O2minus::~O2minus(){

}

air9eed::Ominus::Ominus(){
  m_name   = "Ominus";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = false;
  m_mobile = true;

    std::string str;
  ParmParse pp("air9eed");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
}

air9eed::Ominus::~Ominus(){

}

air9eed::photon_one::photon_one(){
  m_name   = "photon_one";
  m_A      = 1.12E-4; // Default parameters
  m_lambda = 4.15E-2;

  { // Override from input script
    ParmParse pp("air9eed");
    pp.query("photon1_A_coeff",      m_A);
    pp.query("photon1_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air9eed::photon_one::~photon_one(){

}

air9eed::photon_two::photon_two(){
  m_name   = "photon_two";
  m_A      = 2.88E-3; // Default parameters
  m_lambda = 1.09E-1;

  { // Override from input script
    ParmParse pp("air9eed");
    pp.query("photon2_A_coeff",      m_A);
    pp.query("photon2_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air9eed::photon_two::~photon_two(){

}

air9eed::photon_three::photon_three(){
  m_name   = "photon_three";
  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Override from input script
    ParmParse pp("air9eed");
    pp.query("photon3_A_coeff",      m_A);
    pp.query("photon3_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air9eed::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air9eed::photon_three::~photon_three(){

}

Real air9eed::eed::initial_data(const RealVect a_pos, const Real a_time) const{
  return 1.E10;
}

Real air9eed::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air9eed::N2plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air9eed::N4plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air9eed::O2plus::initial_data(const RealVect a_pos, const Real a_time) const{
  return m_initial_ionization;
}

Real air9eed::O4plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air9eed::O2plusN2::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air9eed::O2minus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air9eed::Ominus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air9eed::photon_one::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air9eed::photon_one::get_lambda() const {
  return m_lambda;
}

Real air9eed::photon_one::get_A() const {
  return m_A;
}

Real air9eed::photon_one::get_pO2() const {
  return m_pO2;
}

Real air9eed::photon_two::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air9eed::photon_two::get_lambda() const {
  return m_lambda;
}

Real air9eed::photon_two::get_A() const {
  return m_A;
}

Real air9eed::photon_two::get_pO2() const {
  return m_pO2;
}

Real air9eed::photon_three::get_kappa(const RealVect a_pos) const{
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air9eed::photon_three::get_lambda() const {
  return m_lambda;
}

Real air9eed::photon_three::get_A() const {
  return m_A;
}

Real air9eed::photon_three::get_pO2() const {
  return m_pO2;
}

