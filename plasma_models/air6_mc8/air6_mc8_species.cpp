#include "air6_mc8.H"
#include "air6_mc8_species.H"

air6_mc8::electron::electron(){
  m_name = "electron";
  
  ParmParse pp("air6_mc8");
  std::string str;
  
  pp.get("mobile_electrons", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_diffusive = (str == "true") ? true : false;
}

air6_mc8::M_plus::M_plus(){
  m_name = "M_plus";
  
  ParmParse pp("air6_mc8");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

air6_mc8::M_minus::M_minus(){
  m_name = "M_minus";
  
  ParmParse pp("air6_mc8");
  std::string str;
  
  pp.get("mobile_ions", str);    m_mobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_diffusive = (str == "true") ? true : false;
}

air6_mc8::N2_c4v0::N2_c4v0(){
  m_name      = "N2_c4v0";
  m_mobile    = false;
  m_diffusive = false;
}

air6_mc8::N2_c4v1::N2_c4v1(){
  m_name      = "N2_c4v1";
  m_mobile    = false;
  m_diffusive = false;
}

air6_mc8::N2_b1v1::N2_b1v1(){
  m_name      = "N2_b1v1";
  m_mobile    = false;
  m_diffusive = false;
}

air6_mc8::phot_c4v0_X1v0::phot_c4v0_X1v0(){
  m_name      = "phot_c4v0_X1v0";

  Real p;
  ParmParse pp("air6_mc8");


  pp.get("pressure", p);
  pp.get("c4v0_X1v0_beer", m_kappa);

  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_c4v0_X1v1::phot_c4v0_X1v1(){
  m_name      = "phot_c4v0_X1v1";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v0_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_c4v1_X1v0::phot_c4v1_X1v0(){
  m_name      = "phot_c4v1_X1v0";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_c4v1_X1v1::phot_c4v1_X1v1(){
  m_name      = "phot_c4v1_X1v1";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_c4v1_X1v2::phot_c4v1_X1v2(){
  m_name      = "phot_c4v1_X1v2";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v2_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_c4v1_X1v3::phot_c4v1_X1v3(){
  m_name      = "phot_c4v1_X1v3";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("c4v1_X1v3_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_b1v1_X1v0::phot_b1v1_X1v0(){
  m_name      = "phot_b1v1_X1v0";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v0_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}

air6_mc8::phot_b1v1_X1v1::phot_b1v1_X1v1(){
  m_name      = "phot_b1v1_X1v1";

  Real p;
  ParmParse pp("air6_mc8");
  
  pp.get("pressure", p);
  pp.get("b1v1_X1v1_beer", m_kappa);
  
  m_kappa = 1./(m_kappa*p);
}
