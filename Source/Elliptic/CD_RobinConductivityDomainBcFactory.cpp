/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RobinConductivityDomainBcFactory.cpp
  @brief  Implementation of CD_RobinConductivityDomainBcFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RobinConductivityDomainBcFactory.H>
#include <CD_NamespaceHeader.H>

RobinConductivityDomainBcFactory::RobinConductivityDomainBcFactory(){
  this->setCoefficientss(1., 1., 0.);
}

RobinConductivityDomainBcFactory::~RobinConductivityDomainBcFactory(){

}

void RobinConductivityDomainBcFactory::setCoefficientss(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aCoefficient = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void RobinConductivityDomainBcFactory::setCoefficientss(const RefCountedPtr<RobinCoefficients> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void RobinConductivityDomainBcFactory::setCoefficientss(const EBAMRFluxData& a_aco,
						 const EBAMRFluxData& a_bco,
						 const EBAMRFluxData& a_rhs){

  MayDay::Abort("RobinConductivityDomainBcFactory::setCoefficientss - data based coefficients not supported (yet)");

  m_aCoefficientdata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}

BaseDomainBC* RobinConductivityDomainBcFactory::create(const ProblemDomain& a_domain,
						       const EBISLayout&    a_ebisl,
						       const RealVect&      a_dx){

  RobinConductivityDomainBc* fresh = new RobinConductivityDomainBc();

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

#include <CD_NamespaceFooter.H>
