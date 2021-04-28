#include "air3_bourdon.H"
#include "air3_bourdon_species.H"
#include "units.H"
#include "ParmParse.H"

#include <chrono>

using namespace physics::cdr_plasma;

air3_bourdon::electron::electron(){
  m_name = "electron";
  m_charge = -1;
  ParmParse pp("air3_bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons",    m_mobile);    
  pp.get("diffusive_electrons", m_diffusive);
  pp.get("uniform_density",     m_uniform_density);
  pp.get("seed_density",        m_seed_density);
  pp.get("seed_radius",         m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_bourdon::M_plus::M_plus(){
  m_name = "M_plus";
  m_charge = 1;
  ParmParse pp("air3_bourdon");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions",     m_mobile);    
  pp.get("diffusive_ions",  m_diffusive);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density",    m_seed_density);
  pp.get("seed_radius",     m_seed_rad);
  
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_bourdon::M_minus::M_minus(){
  m_name = "M_minus";
  m_charge = -1;
  ParmParse pp("air3_bourdon");
  std::string str;
  
  pp.get("mobile_ions", m_mobile);    
  pp.get("diffusive_ions", m_diffusive);
}

Real air3_bourdon::electron::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real air3_bourdon::M_plus::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

air3_bourdon::photon_one::photon_one(){
  m_name     = "photon_one";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("air3_bourdon");
  pp.get("photon1_A_coeff",      m_A);
  pp.get("photon1_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

air3_bourdon::photon_two::photon_two(){
  m_name     = "photon_two";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  
  ParmParse pp("air3_bourdon");
  pp.get("photon2_A_coeff",      m_A);
  pp.get("photon2_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}

air3_bourdon::photon_three::photon_three(){
  m_name     = "photon_three";
  m_constant = true;

  Real O2_frac = 0.2;
  Real pressure;
  ParmParse pp("air3_bourdon");
  pp.get("photon3_A_coeff",      m_A);
  pp.get("photon3_lambda_coeff", m_lambda);
  pp.get("pressure",             pressure);
  
  m_pO2 = pressure*O2_frac*units::s_atm2pascal;
}
