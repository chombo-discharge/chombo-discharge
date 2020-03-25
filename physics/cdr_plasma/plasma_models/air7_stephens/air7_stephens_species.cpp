/*!
  @file air_stephens_species.H
*/

#include "air7_stephens.H"
#include "air7_stephens_species.H"
#include "units.H"

#include <chrono>

using namespace physics::cdr_plasma;

air7_stephens::electron::electron(){
  m_name = "electron";
  m_charge = -1;
  ParmParse pp("air7_stephens");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;

  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air7_stephens::N2plus::N2plus(){
  m_name = "N2plus";
  m_charge = 1;
  ParmParse pp("air7_stephens");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;

  pp.get("frac_N2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air7_stephens::O2plus::O2plus(){
  m_name = "O2plus";
  m_charge = 1;
  ParmParse pp("air7_stephens");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;

  pp.get("frac_O2", m_frac);
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air7_stephens::N4plus::N4plus(){
  m_name = "N4plus";
  m_charge = 1;
  ParmParse pp("air7_stephens");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

air7_stephens::O4plus::O4plus(){
  m_name = "O4plus";
  m_charge = 1;
  ParmParse pp("air7_stephens");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

air7_stephens::O2plusN2::O2plusN2(){
  m_name = "O2plusN2";
  m_charge = 1;
  ParmParse pp("air7_stephens");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

air7_stephens::O2minus::O2minus(){
  m_name = "O2minus";
  m_charge = -1;
  ParmParse pp("air7_stephens");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

Real air7_stephens::electron::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real air7_stephens::N2plus::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
}

Real air7_stephens::O2plus::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_frac*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
}

air7_stephens::phot_c4v0_X1v0::phot_c4v0_X1v0(){
  m_name      = "phot_c4v0_X1v0";

  Real p;
  ParmParse pp("air7_stephens");


  pp.get("pressure", p);
  pp.get("c4v0_X1v0_beer", m_kappa);
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_c4v0_X1v1::phot_c4v0_X1v1(){
  m_name      = "phot_c4v0_X1v1";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("c4v0_X1v1_beer", m_kappa);
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_c4v1_X1v0::phot_c4v1_X1v0(){
  m_name      = "phot_c4v1_X1v0";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_c4v1_X1v1::phot_c4v1_X1v1(){
  m_name      = "phot_c4v1_X1v1";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_c4v1_X1v2::phot_c4v1_X1v2(){
  m_name      = "phot_c4v1_X1v2";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v2_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_c4v1_X1v3::phot_c4v1_X1v3(){
  m_name      = "phot_c4v1_X1v3";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v3_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_b1v1_X1v0::phot_b1v1_X1v0(){
  m_name      = "phot_b1v1_X1v0";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air7_stephens::phot_b1v1_X1v1::phot_b1v1_X1v1(){
  m_name      = "phot_b1v1_X1v1";

  Real p;
  ParmParse pp("air7_stephens");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}
