/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir3BourdonSpecies.cpp
  @brief  Implementation of CD_CdrPlasmaAir3BourdonSpecies.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaAir3Bourdon.H>
#include <CD_CdrPlasmaAir3BourdonSpecies.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaAir3Bourdon::Electron::Electron(){
  m_name = "Electron";
  m_chargeNumber = -1;
  ParmParse pp("CdrPlasmaAir3Bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons",    m_isMobile);    
  pp.get("diffusive_electrons", m_isDiffusive);
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir3Bourdon::MPlus::MPlus(){
  m_name = "MPlus";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir3Bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions",     m_isMobile);    
  pp.get("diffusive_ions",  m_isDiffusive);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density",    m_seed_density);
  pp.get("seed_radius",     m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir3Bourdon::MMinus::MMinus(){
  m_name = "MMinus";
  m_chargeNumber = -1;
  ParmParse pp("CdrPlasmaAir3Bourdon");
  std::string str;
  
  pp.get("mobile_ions",    m_isMobile);    
  pp.get("diffusive_ions", m_isDiffusive);
}

Real CdrPlasmaAir3Bourdon::Electron::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real CdrPlasmaAir3Bourdon::MPlus::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

CdrPlasmaAir3Bourdon::PhotonOne::PhotonOne(){
  m_name     = "PhotonOne";

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("CdrPlasmaAir3Bourdon");
  pp.get("Photon1_A_coeff",      m_A);
  pp.get("Photon1_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

CdrPlasmaAir3Bourdon::PhotonTwo::PhotonTwo(){
  m_name     = "PhotonTwo";

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("CdrPlasmaAir3Bourdon");
  pp.get("Photon2_A_coeff",      m_A);
  pp.get("Photon2_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

CdrPlasmaAir3Bourdon::PhotonThree::PhotonThree(){
  m_name     = "PhotonThree";

  Real O2_frac = 0.2;
  Real pressure;
  ParmParse pp("CdrPlasmaAir3Bourdon");
  pp.get("Photon3_A_coeff",      m_A);
  pp.get("Photon3_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

#include <CD_NamespaceFooter.H>
