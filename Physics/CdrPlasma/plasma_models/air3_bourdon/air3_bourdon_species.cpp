#include "air3_bourdon.H"
#include "air3_bourdon_species.H"
#include <CD_Units.H>
#include "ParmParse.H"

#include <chrono>

#include "CD_NamespaceHeader.H"
using namespace Physics::CdrPlasma;

air3_bourdon::electron::electron(){
  m_name = "electron";
  m_chargeNumber = -1;
  ParmParse pp("air3_bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons",    m_isMobile);    
  pp.get("diffusive_electrons", m_isDiffusive);
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_bourdon::M_plus::M_plus(){
  m_name = "M_plus";
  m_chargeNumber = 1;
  ParmParse pp("air3_bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions",     m_isMobile);    
  pp.get("diffusive_ions",  m_isDiffusive);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density",    m_seed_density);
  pp.get("seed_radius",     m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_bourdon::M_minus::M_minus(){
  m_name = "M_minus";
  m_chargeNumber = -1;
  ParmParse pp("air3_bourdon");
  std::string str;
  
  pp.get("mobile_ions",    m_isMobile);    
  pp.get("diffusive_ions", m_isDiffusive);
}

Real air3_bourdon::electron::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real air3_bourdon::M_plus::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

air3_bourdon::Photon_one::Photon_one(){
  m_name     = "Photon_one";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("air3_bourdon");
  pp.get("Photon1_A_coeff",      m_A);
  pp.get("Photon1_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

air3_bourdon::Photon_two::Photon_two(){
  m_name     = "Photon_two";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("air3_bourdon");
  pp.get("Photon2_A_coeff",      m_A);
  pp.get("Photon2_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}

air3_bourdon::Photon_three::Photon_three(){
  m_name     = "Photon_three";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  ParmParse pp("air3_bourdon");
  pp.get("Photon3_A_coeff",      m_A);
  pp.get("Photon3_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*Units::atm2pascal;
}
#include "CD_NamespaceFooter.H"
