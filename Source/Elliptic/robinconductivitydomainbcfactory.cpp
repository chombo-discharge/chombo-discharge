/*!
  @file   robinconductivitydomainbcfactory.cpp
  @brief  Implementation of robinconductivitydomainbcfactory.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "robinconductivitydomainbcfactory.H"

#include "CD_NamespaceHeader.H"

robinconductivitydomainbcfactory::robinconductivitydomainbcfactory(){
  this->setCoefficientss(1., 1., 0.);
}

robinconductivitydomainbcfactory::~robinconductivitydomainbcfactory(){

}

void robinconductivitydomainbcfactory::setCoefficientss(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aCoefficient = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivitydomainbcfactory::setCoefficientss(const RefCountedPtr<robin_coef> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void robinconductivitydomainbcfactory::setCoefficientss(const EBAMRFluxData& a_aco,
						 const EBAMRFluxData& a_bco,
						 const EBAMRFluxData& a_rhs){

  MayDay::Abort("robinconductivitydomainbcfactory::setCoefficientss - data based coefficients not supported (yet)");

  m_aCoefficientdata = a_aco;
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
    fresh->setCoefficientss(m_aCoefficient, m_bco, m_rhs);
  }
  else if(m_func_coeff){

    fresh->setCoefficientss(m_robinco);
  }
  else if(m_data_coeff){
    for (int lvl = 0; lvl < m_aCoefficientdata.size(); lvl++){
      if(a_domain == (m_aCoefficientdata[lvl]->disjointBoxLayout()).physDomain()){
	fresh->setCoefficientss(m_aCoefficientdata[lvl], m_bcodata[lvl], m_rhsdata[lvl]);
      }
    }
  }
  else{
    MayDay::Abort("robinconductivityebbcfactory::create - must set coefficients first");
  }

  return static_cast<BaseDomainBC*> (fresh);
}
#include "CD_NamespaceFooter.H"
