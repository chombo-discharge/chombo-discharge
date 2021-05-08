/*!
  @file robinconductivityebbc.cpp
  @brief Implementation of robinconductivityebbc.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "robinconductivityebbcfactory.H"

#include "CD_NamespaceHeader.H"

robinconductivityebbcfactory::robinconductivityebbcfactory(const RealVect a_origin){
  this->set_type(stencil_type::taylor);
  this->set_coefs(1., 1., 0.);

  m_origin = a_origin;
}

robinconductivityebbcfactory::~robinconductivityebbcfactory(){
}

void robinconductivityebbcfactory::set_type(const stencil_type a_type){
  if(a_type == stencil_type::taylor || a_type == stencil_type::lsq){
    m_type = a_type;
  }
  else{
    MayDay::Abort("robinconductivityebbcfactory::set_type - only taylor and least squares supported for now");
  }
}

void robinconductivityebbcfactory::set_coefs(const Real a_aco, const Real a_bco, const Real a_rhs){
  m_aco = a_aco;
  m_bco = a_bco;
  m_rhs = a_rhs;

  m_const_coeff = true;
  m_func_coeff  = false;
  m_data_coeff  = false;
}

void robinconductivityebbcfactory::set_coefs(const RefCountedPtr<robin_coef> a_robinco){
  m_robinco = a_robinco;

  m_const_coeff = false;
  m_func_coeff  = true;
  m_data_coeff  = false;
}

void robinconductivityebbcfactory::set_coefs(const EBAMRIVData& a_aco,
					     const EBAMRIVData& a_bco,
					     const EBAMRIVData& a_rhs){

  MayDay::Abort("robinconductivityebbcfactory::set_coefs - data based coefficients not supported (yet)");

  m_acodata = a_aco;
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

  fresh->set_type(m_type);

  return fresh;
}
#include "CD_NamespaceFooter.H"
