/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RobinConductivityEbBcFactory.cpp
  @brief  Implementation of CD_RobinConductivityEbBcFactory.H
  @author Robert Marskar
*/

#include <CD_RobinConductivityEbBcFactory.H>
#include <CD_NamespaceHeader.H>

RobinConductivityEbBcFactory::RobinConductivityEbBcFactory(const RealVect a_origin){
  this->setStencilType(IrregStencil::StencilType::TaylorExtrapolation);
  this->setCoefficients(1., 1., 0.);

  m_origin = a_origin;
}

RobinConductivityEbBcFactory::~RobinConductivityEbBcFactory(){
}

void RobinConductivityEbBcFactory::setStencilType(const IrregStencil::StencilType a_type){
  if(a_type == IrregStencil::StencilType::TaylorExtrapolation || a_type == IrregStencil::StencilType::LeastSquares){
    m_type = a_type;
  }
  else{
    MayDay::Abort("RobinConductivityEbBcFactory::set_type - only taylor and least squares supported for now");
  }
}

void RobinConductivityEbBcFactory::setCoefficients(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void RobinConductivityEbBcFactory::setCoefficients(const RefCountedPtr<RobinCoefficients> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void RobinConductivityEbBcFactory::setCoefficients(const EBAMRIVData& a_aco,
						   const EBAMRIVData& a_bco,
						   const EBAMRIVData& a_rhs){

  MayDay::Abort("RobinConductivityEbBcFactory::setCoefficients - data based coefficients not supported (yet)");

  m_acodata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}

RobinConductivityEbBc* RobinConductivityEbBcFactory::create(const ProblemDomain& a_domain,
							    const EBISLayout&    a_ebisl,
							    const RealVect&      a_dx,
							    const IntVect*       a_ghost_phi,
							    const IntVect*       a_ghost_rhs){

  RobinConductivityEbBc* fresh = new RobinConductivityEbBc(a_dx, m_origin);

  if(m_const_coeff){
    fresh->setCoefficients(m_aco, m_bco, m_rhs);
  }
  else if(m_func_coeff){

    fresh->setCoefficients(m_robinco);
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

  fresh->setStencilType(m_type);

  return fresh;
}

#include <CD_NamespaceFooter.H>
