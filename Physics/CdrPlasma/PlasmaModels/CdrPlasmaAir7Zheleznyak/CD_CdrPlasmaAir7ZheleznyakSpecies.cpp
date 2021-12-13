/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaAir7ZheleznyakSpecies.cpp
  @brief  Implementation of CD_CdrPlasmaAir7ZheleznyakSpecies.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Our includes
#include <CD_CdrPlasmaAir7Zheleznyak.H>
#include <CD_CdrPlasmaAir7ZheleznyakSpecies.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaAir7Zheleznyak::Electron::Electron(){
  m_name = "Electron";
  m_chargeNumber = -1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive = (str == "true") ? true : false;

  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Zheleznyak::N2plus::N2plus(){
  m_name = "N2plus";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;

  pp.get("frac_N2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Zheleznyak::O2plus::O2plus(){
  m_name = "O2plus";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;

  pp.get("frac_O2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

CdrPlasmaAir7Zheleznyak::N4plus::N4plus(){
  m_name = "N4plus";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Zheleznyak::O4plus::O4plus(){
  m_name = "O4plus";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Zheleznyak::O2plusN2::O2plusN2(){
  m_name = "O2plusN2";
  m_chargeNumber = 1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
}

CdrPlasmaAir7Zheleznyak::O2minus::O2minus(){
  m_name = "O2minus";
  m_chargeNumber = -1;
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  std::string str;
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
}

Real CdrPlasmaAir7Zheleznyak::Electron::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real CdrPlasmaAir7Zheleznyak::N2plus::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
}

Real CdrPlasmaAir7Zheleznyak::O2plus::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
}

CdrPlasmaAir7Zheleznyak::uv_Photon::uv_Photon(){
  m_name   = "uv_Photon";

  Real pressure, O2_frac;
  
  ParmParse pp("CdrPlasmaAir7Zheleznyak");
  pp.get("photoi_f1",      m_f1);
  pp.get("photoi_f2",      m_f2);
  pp.get("photoi_K1",      m_K1);
  pp.get("photoi_K2",      m_K2);

  pp.get("frac_O2",  O2_frac);
  pp.get("pressure", pressure);
  pp.get("photoi_seed",     m_seed);

  // Convert units
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
  m_K1  = m_K1*m_pO2;
  m_K2  = m_K2*m_pO2;

  // Seed the RNG
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();
  m_rng = new std::mt19937_64(m_seed);
  m_udist01 = new std::uniform_real_distribution<Real>(0.0, 1.0);
}

CdrPlasmaAir7Zheleznyak::uv_Photon::~uv_Photon(){
  
}

Real CdrPlasmaAir7Zheleznyak::uv_Photon::getAbsorptionCoefficient(const RealVect a_pos) const {
  return getRandomKappa();
}

Real CdrPlasmaAir7Zheleznyak::uv_Photon::getRandomKappa() const {
  const Real f = m_f1 + (*m_udist01)(*m_rng)*(m_f2 - m_f1);
  return m_K1*pow(m_K2/m_K1, (f-m_f1)/(m_f2-m_f1));
}

#include <CD_NamespaceFooter.H>
