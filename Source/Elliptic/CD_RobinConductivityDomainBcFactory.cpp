/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
  this->setCoefficients(1., 1., 0.);
}

RobinConductivityDomainBcFactory::~RobinConductivityDomainBcFactory(){

}

void RobinConductivityDomainBcFactory::setCoefficients(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void RobinConductivityDomainBcFactory::setCoefficients(const RefCountedPtr<RobinCoefficients> a_robinco){
  m_robinCoefficients = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void RobinConductivityDomainBcFactory::setCoefficients(const EBAMRFluxData& a_aco,
						 const EBAMRFluxData& a_bco,
						 const EBAMRFluxData& a_rhs){

  MayDay::Abort("RobinConductivityDomainBcFactory::setCoefficients - data based coefficients not supported (yet)");

  m_acodata = a_aco;
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
    fresh->setCoefficients(m_aco, m_bco, m_rhs);
  }
  else if(m_func_coeff){

    fresh->setCoefficients(m_robinCoefficients);
  }
  else if(m_data_coeff){
    for (int lvl = 0; lvl < m_acodata.size(); lvl++){
      if(a_domain == (m_acodata[lvl]->disjointBoxLayout()).physDomain()){
	fresh->setCoefficients(m_acodata[lvl], m_bcodata[lvl], m_rhsdata[lvl]);
      }
    }
  }
  else{
    MayDay::Abort("RobinConductivityEbBcFactory::create - must set coefficients first");
  }

  return static_cast<BaseDomainBC*> (fresh);
}

#include <CD_NamespaceFooter.H>
