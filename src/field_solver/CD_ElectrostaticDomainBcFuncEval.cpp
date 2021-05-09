/*!
  @file   CD_ElectrostaticDomainBcFuncEval.cpp
  @brief  Implementation of CD_ElectrostaticDomainBcFuncEval.H
  @author Robert Marskar
  @date   May 2021
*/

#include "CD_ElectrostaticDomainBcFuncEval.H"
#include "CD_NamespaceHeader.H"

ElectrostaticDomainBcFuncEval::ElectrostaticDomainBcFuncEval(){
  m_isFunctionSet = false;
}

ElectrostaticDomainBcFuncEval::~ElectrostaticDomainBcFuncEval(){
}

void ElectrostaticDomainBcFuncEval::setTime(const Real a_time){
  m_time = a_time;
}

void ElectrostaticDomainBcFuncEval::setFunction(ElectrostaticDomainBc::BcFunction a_bcFunc){
  m_bcFunc        = a_bcFunc;
  m_isFunctionSet = true;
}

Real ElectrostaticDomainBcFuncEval::value(const RealVect& a_point, const int& a_comp) const {
  if(!m_isFunctionSet) MayDay::Abort("ElectrostaticDomainBcFuncEval -- function has not been set!");
  
  return m_bcFunc(a_point, m_time);
}

Real ElectrostaticDomainBcFuncEval::derivative(const RealVect& a_point, const int& a_comp, const int& a_derivDir) const {
  MayDay::Abort("ElectrostaticDomainBcFuncEval::derivative -- calling this is an error. How did you get here?");
}

#include "CD_NamespaceFooter.H"
