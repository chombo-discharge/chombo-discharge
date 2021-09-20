/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CdrPlasmaAir9EedBourdonSpecies.cpp
  @brief  Implementation of CdrPlasmaAir9EedBourdonSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaAir9EedBourdonSpecies.H>
#include <CD_Units.H> 
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaAir9EedBourdon::eed::eed(){
  m_name         = "eed";
  m_chargeNumber = 0;
  m_isDiffusive  = true;
  m_isMobile     = true;

  // Get gas parameters
  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
}

CdrPlasmaAir9EedBourdon::eed::~eed(){

}

CdrPlasmaAir9EedBourdon::Electron::Electron(){
  m_name   = "Electron";
  m_chargeNumber = -1;
  m_isDiffusive = true;
  m_isMobile = true;

  {// Get initial parameter
    ParmParse pp("CdrPlasmaAir9EedBourdon");
    pp.get("initial_ionization", m_initial_ionization);
  }
}

CdrPlasmaAir9EedBourdon::Electron::~Electron(){

}

CdrPlasmaAir9EedBourdon::N2plus::N2plus() {
  m_name   = "N2plus";
  m_chargeNumber = 1;
  m_isDiffusive = false;
  m_isMobile = true;

  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("initial_ionization", m_initial_ionization);
  m_initial_ionization *= N2frac;

  std::string str;
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::N2plus::~N2plus(){

}

CdrPlasmaAir9EedBourdon::N4plus::N4plus(){
  m_name   = "N4plus";
  m_chargeNumber = 1;
  m_isDiffusive = false;
  m_isMobile = true;

  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::N4plus::~N4plus(){

}

CdrPlasmaAir9EedBourdon::O2plus::O2plus(){
  m_name   = "O2plus";
  m_chargeNumber = 1;
  m_isDiffusive = false;
  m_isMobile = true;

  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  pp.get("initial_ionization", m_initial_ionization);
  m_initial_ionization *= O2frac;

  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::O2plus::~O2plus(){

}

CdrPlasmaAir9EedBourdon::O4plus::O4plus(){
  m_name   = "O4plus";
  m_chargeNumber = 1;
  m_isDiffusive = false;
  m_isMobile = true;

  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::O4plus::~O4plus(){

}

CdrPlasmaAir9EedBourdon::O2plusN2::O2plusN2() {
  m_name   = "O2plusN2";
  m_chargeNumber = 1;
  m_isDiffusive = false;
  m_isMobile = true;

  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::O2plusN2::~O2plusN2(){

}

CdrPlasmaAir9EedBourdon::O2minus::O2minus(){
  m_name   = "O2minus";
  m_chargeNumber = -1;
  m_isDiffusive = false;
  m_isMobile = true;

  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::O2minus::~O2minus(){

}

CdrPlasmaAir9EedBourdon::Ominus::Ominus(){
  m_name   = "Ominus";
  m_chargeNumber = -1;
  m_isDiffusive = false;
  m_isMobile = true;

  std::string str;
  ParmParse pp("CdrPlasmaAir9EedBourdon");
  if(pp.contains("mobile_ions")){
    pp.get("mobile_ions", str);
    if(str == "true"){
      m_isMobile = true;
    }
    else if(str == "false"){
      m_isMobile = false;
    }
  }
}

CdrPlasmaAir9EedBourdon::Ominus::~Ominus(){

}

CdrPlasmaAir9EedBourdon::PhotonOne::PhotonOne(){
  m_name   = "PhotonOne";
  m_A      = 1.12E-4; // Default parameters
  m_lambda = 4.15E-2;

  { // Override from input script
    ParmParse pp("CdrPlasmaAir9EedBourdon");
    pp.query("Photon1_A_coeff",      m_A);
    pp.query("Photon1_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

CdrPlasmaAir9EedBourdon::PhotonOne::~PhotonOne(){

}

CdrPlasmaAir9EedBourdon::PhotonTwo::PhotonTwo(){
  m_name   = "PhotonTwo";
  m_A      = 2.88E-3; // Default parameters
  m_lambda = 1.09E-1;

  { // Override from input script
    ParmParse pp("CdrPlasmaAir9EedBourdon");
    pp.query("Photon2_A_coeff",      m_A);
    pp.query("Photon2_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

CdrPlasmaAir9EedBourdon::PhotonTwo::~PhotonTwo(){

}

CdrPlasmaAir9EedBourdon::PhotonThree::PhotonThree(){
  m_name   = "PhotonThree";
  m_A      = 2.76E-1;
  m_lambda = 6.69E-1;

  { // Override from input script
    ParmParse pp("CdrPlasmaAir9EedBourdon");
    pp.query("Photon3_A_coeff",      m_A);
    pp.query("Photon3_lambda_coeff", m_lambda);
  }

  // Get gas stuff from input script
  Real Tg, p, N, O2frac, N2frac;
  CdrPlasmaAir9EedBourdon::parseGas_parameters(Tg, p, N, O2frac, N2frac);
  m_pO2 = p*O2frac;
}

CdrPlasmaAir9EedBourdon::PhotonThree::~PhotonThree(){

}

Real CdrPlasmaAir9EedBourdon::eed::initialData(const RealVect a_pos, const Real a_time) const{
  return 1.E10;
}

Real CdrPlasmaAir9EedBourdon::Electron::initialData(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real CdrPlasmaAir9EedBourdon::N2plus::initialData(const RealVect a_pos, const Real a_time) const {
  return m_initial_ionization;
}

Real CdrPlasmaAir9EedBourdon::N4plus::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real CdrPlasmaAir9EedBourdon::O2plus::initialData(const RealVect a_pos, const Real a_time) const{
  return m_initial_ionization;
}

Real CdrPlasmaAir9EedBourdon::O4plus::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real CdrPlasmaAir9EedBourdon::O2plusN2::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real CdrPlasmaAir9EedBourdon::O2minus::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real CdrPlasmaAir9EedBourdon::Ominus::initialData(const RealVect a_pos, const Real a_time) const {
  return 0.0;
}

Real CdrPlasmaAir9EedBourdon::PhotonOne::getKappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real CdrPlasmaAir9EedBourdon::PhotonOne::getLambda() const {
  return m_lambda;
}

Real CdrPlasmaAir9EedBourdon::PhotonOne::getA() const {
  return m_A;
}

Real CdrPlasmaAir9EedBourdon::PhotonOne::getPO2() const {
  return m_pO2;
}

Real CdrPlasmaAir9EedBourdon::PhotonTwo::getKappa(const RealVect a_pos) const {
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real CdrPlasmaAir9EedBourdon::PhotonTwo::getLambda() const {
  return m_lambda;
}

Real CdrPlasmaAir9EedBourdon::PhotonTwo::getA() const {
  return m_A;
}

Real CdrPlasmaAir9EedBourdon::PhotonTwo::getPO2() const {
  return m_pO2;
}

Real CdrPlasmaAir9EedBourdon::PhotonThree::getKappa(const RealVect a_pos) const{
  return m_lambda*m_pO2/(sqrt(3.0));
}

Real CdrPlasmaAir9EedBourdon::PhotonThree::getLambda() const {
  return m_lambda;
}

Real CdrPlasmaAir9EedBourdon::PhotonThree::getA() const {
  return m_A;
}

Real CdrPlasmaAir9EedBourdon::PhotonThree::getPO2() const {
  return m_pO2;
}

#include <CD_NamespaceFooter.H>
