/*!
  @file   robinconductivitydomainbcfactory.cpp
  @brief  Implementation of robinconductivitydomainbcfactory.H
  @author Robert Marskar
  @date   Jan. 2018
*/


#include "robinconductivitydomainbcfactory.H"

robinconductivitydomainbcfactory::robinconductivitydomainbcfactory(){
  this->set_coefs(1., 1., 0.);
}


robinconductivitydomainbcfactory::~robinconductivitydomainbcfactory(){

}

void robinconductivitydomainbcfactory::set_coefs(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivitydomainbcfactory::set_coefs(const RefCountedPtr<robin_coef> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void robinconductivitydomainbcfactory::set_coefs(const EBAMRFluxData& a_aco,
						 const EBAMRFluxData& a_bco,
						 const EBAMRFluxData& a_rhs){

  MayDay::Abort("robinconductivitydomainbcfactory::set_coefs - data based coefficients not supported (yet)");

  m_acodata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}


BaseDomainBC* robinconductivitydomainbcfactory::create(const ProblemDomain& a_domain,
						       const EBISLayout&    a_ebisl,
						       const RealVect&      a_dx){

  robinconductivitydomainbc* fresh = new robinconductivitydomainbc();

  if(m_const_coeff){
    fresh->set_coefs(m_aco, m_bco, m_rhs);
  }
  else if(m_func_coeff){

    fresh->set_coefs(m_robinco);
  }
  else if(m_data_coeff){
    for (int lvl = 0; lvl < m_acodata.size(); lvl++){
      if(a_domain == (m_acodata[lvl]->disjointBoxLayout()).physDomain()){
	fresh->set_coefs(m_acodata[lvl], m_bcodata[lvl], m_rhsdata[lvl]);
      }
    }
  }
  else{
    MayDay::Abort("robinconductivityebbcfactory::create - must set coefficients first");
  }


  return static_cast<BaseDomainBC*> (fresh);
}
