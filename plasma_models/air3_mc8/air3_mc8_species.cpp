#include "air3_mc8.H"
#include "air3_mc8_species.H"

air3_mc8::electron::electron(){
  m_name = "electron";
  m_charge = -1;
  ParmParse pp("air3_mc8");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;

  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_mc8::M_plus::M_plus(){
  m_name = "M_plus";
  m_charge = 1;
  ParmParse pp("air3_mc8");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;

  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air3_mc8::M_minus::M_minus(){
  m_name = "M_minus";
  m_charge = -1;
  ParmParse pp("air3_mc8");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

Real air3_mc8::electron::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real air3_mc8::M_plus::initial_data(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

air3_mc8::phot_c4v0_X1v0::phot_c4v0_X1v0(){
  m_name      = "phot_c4v0_X1v0";

  Real p;
  ParmParse pp("air3_mc8");


  pp.get("pressure", p);
  pp.get("c4v0_X1v0_beer", m_kappa);
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_c4v0_X1v1::phot_c4v0_X1v1(){
  m_name      = "phot_c4v0_X1v1";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v0_X1v1_beer", m_kappa);
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_c4v1_X1v0::phot_c4v1_X1v0(){
  m_name      = "phot_c4v1_X1v0";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_c4v1_X1v1::phot_c4v1_X1v1(){
  m_name      = "phot_c4v1_X1v1";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_c4v1_X1v2::phot_c4v1_X1v2(){
  m_name      = "phot_c4v1_X1v2";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v2_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_c4v1_X1v3::phot_c4v1_X1v3(){
  m_name      = "phot_c4v1_X1v3";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v3_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_b1v1_X1v0::phot_b1v1_X1v0(){
  m_name      = "phot_b1v1_X1v0";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air3_mc8::phot_b1v1_X1v1::phot_b1v1_X1v1(){
  m_name      = "phot_b1v1_X1v1";

  Real p;
  ParmParse pp("air3_mc8");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

