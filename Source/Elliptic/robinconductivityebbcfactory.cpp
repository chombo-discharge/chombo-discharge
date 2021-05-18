/*!
  @file robinconductivityebbc.cpp
  @brief Implementation of robinconductivityebbc.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "robinconductivityebbcfactory.H"

#include "CD_NamespaceHeader.H"

robinconductivityebbcfactory::robinconductivityebbcfactory(const RealVect a_origin){
  this->set_type(IrregStencil::StencilType::TaylorExtrapolation);
  this->setCoefficientss(1., 1., 0.);

  m_origin = a_origin;
}

robinconductivityebbcfactory::~robinconductivityebbcfactory(){
}

void robinconductivityebbcfactory::set_type(const IrregStencil::StencilType a_type){
  if(a_type == IrregStencil::StencilType::TaylorExtrapolation || a_type == IrregStencil::StencilType::LeastSquares){
    m_type = a_type;
  }
  else{
    MayDay::Abort("robinconductivityebbcfactory::set_type - only taylor and least squares supported for now");
  }
}

void robinconductivityebbcfactory::setCoefficientss(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aCoefficient = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivityebbcfactory::setCoefficientss(const RefCountedPtr<RobinCoefficients> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void robinconductivityebbcfactory::setCoefficientss(const EBAMRIVData& a_aco,
					     const EBAMRIVData& a_bco,
					     const EBAMRIVData& a_rhs){

  MayDay::Abort("robinconductivityebbcfactory::setCoefficientss - data based coefficients not supported (yet)");

  m_aCoefficientdata = a_aco;
  m_bcodata = a_bco;
  m_rhsdata = a_rhs;

  m_const_coeff = false;
  m_func_coeff  = false;
  m_data_coeff  = true;
}

robinconductivityebbc* robinconductivityebbcfactory::create(const ProblemDomain& a_domain,
							    const EBISLayout&    a_ebisl,
							    const RealVect&      a_dx,
							    const IntVect*       a_ghost_phi,
							    const IntVect*       a_ghost_rhs){

  robinconductivityebbc* fresh = new robinconductivityebbc(a_dx, m_origin);

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

  fresh->set_type(m_type);

  return fresh;
}
#include "CD_NamespaceFooter.H"
