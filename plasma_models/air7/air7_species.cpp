/*!
  @file   air7_species.cpp
  @brief  Implementation of air7_species.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Could definitely see ways to cut down on the typic 
*/

#include "air7_species.H"
#include "units.H" 

#include <ParmParse.H>

air7::electron::electron(){
  m_name   = "electron density";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = true;
  m_mobile = true;

  {// Get initial parameter
    ParmParse pp("air7");
    pp.get("initial_ionization", m_initial_ionization);

    std::string str;
    if(pp.contains("mobile_electrons")){
      pp.get("mobile_electrons", str);
      if(str == "true"){
	m_mobile = true;
      }
      else if(str == "false"){
	m_mobile = false;
      }
    }
    if(pp.contains("diffusive_electrons")){
      pp.get("diffusive_electrons", str);
      if(str == "true"){
	m_diffusive = true;
      }
      else if(str == "false"){
	m_diffusive = false;
      }
    }

  }
}

air7::electron::~electron(){

}

air7::N2plus::N2plus() {
  m_name   = "N2plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  Real Tg, p, N, O2frac, N2frac;
  air7::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  ParmParse pp("air7");
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
  if(pp.contains("diffusive_ions")){
    pp.get("diffusive_ions", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
  }
}

air7::N2plus::~N2plus(){

}

air7::N4plus::N4plus(){
  m_name   = "N4plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air7");
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

air7::N4plus::~N4plus(){

}

air7::O2plus::O2plus(){
  m_name   = "O2plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  Real Tg, p, N, O2frac, N2frac;
  air7::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  std::string str;
  ParmParse pp("air7");
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
  if(pp.contains("diffusive_ions")){
    pp.get("diffusive_ions", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
  }
}

air7::O2plus::~O2plus(){

}

air7::O4plus::O4plus(){
  m_name   = "O4plus";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air7");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
  if(pp.contains("diffusive_ions")){
    pp.get("diffusive_ions", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
  }
}

air7::O4plus::~O4plus(){

}

air7::O2plusN2::O2plusN2() {
  m_name   = "O2plusN2";
  m_unit   = "m-3";
  m_charge = 1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air7");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
  if(pp.contains("diffusive_ions")){
    pp.get("diffusive_ions", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
  }
}

air7::O2plusN2::~O2plusN2(){

}

air7::O2minus::O2minus(){
  m_name   = "O2minus";
  m_unit   = "m-3";
  m_charge = -1;
  m_diffusive = false;
  m_mobile = true;

  std::string str;
  ParmParse pp("air7");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_mobile = true;
    }
    else if(str == "false"){
      m_mobile = false;
    }
  }
  if(pp.contains("diffusive_ions")){
    pp.get("diffusive_ions", str);
    if(str == "true"){
      m_diffusive = true;
    }
    else if(str == "false"){
      m_diffusive = false;
    }
  }
}

air7::O2minus::~O2minus(){

}

air7::photon_one::photon_one(){
  m_name   = "photon_one";
  m_A      = 1.12E-4; // Default parameters
  m_lambda = 4.15E-2;

  { // Override from input script
    ParmParse pp("air7");
    pp.query("photon1_A_coeff",      m_A);
    pp.query("photon1_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air7::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air7::photon_one::~photon_one(){

}

air7::photon_two::photon_two(){
  m_name   = "photon_two";
  m_A      = 2.88E-3; // Default parameters
  m_lambda = 1.09E-1;

  { // Override from input script
    ParmParse pp("air7");
    pp.query("photon2_A_coeff",      m_A);
    pp.query("photon2_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air7::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air7::photon_two::~photon_two(){

}

air7::photon_three::photon_three(){
  m_name   = "photon_three";
  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Override from input script
    ParmParse pp("air7");
    pp.query("photon3_A_coeff",      m_A);
    pp.query("photon3_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  air7::get_gas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

air7::photon_three::~photon_three(){

}

Real air7::electron::initial_data(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air7::N2plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real air7::N4plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air7::O2plus::initial_data(const RealVect a_pos, const Real a_time) const{
  return m_initial_ionization;
}

Real air7::O4plus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air7::O2plusN2::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air7::O2minus::initial_data(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real air7::photon_one::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air7::photon_one::get_lambda() const {
  return m_lambda;
}

Real air7::photon_one::get_A() const {
  return m_A;
}

Real air7::photon_one::get_pO2() const {
  return m_pO2;
}

Real air7::photon_two::get_kappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air7::photon_two::get_lambda() const {
  return m_lambda;
}

Real air7::photon_two::get_A() const {
  return m_A;
}

Real air7::photon_two::get_pO2() const {
  return m_pO2;
}

Real air7::photon_three::get_kappa(const RealVect a_pos) const{
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real air7::photon_three::get_lambda() const {
  return m_lambda;
}

Real air7::photon_three::get_A() const {
  return m_A;
}

Real air7::photon_three::get_pO2() const {
  return m_pO2;
}

