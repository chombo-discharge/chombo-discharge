/*!
  @file   air3_zheleznyak_tagger.cpp
  @brief  Implementation air3_zheleznyak_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air3_zheleznyak_tagger.H"
#include "air3_zheleznyak.H"

#include <ParmParse.H>

air3_zheleznyak_tagger::air3_zheleznyak_tagger(){
  m_name = "air3_zheleznyak_tagger";
  m_num_tracers = 2;
}

air3_zheleznyak_tagger::~air3_zheleznyak_tagger(){

}

void air3_zheleznyak_tagger::parse_options(){
  parseVerbosity();
  parse_boxes();
  parse_buffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_alpha",      m_refi_alpha);
  pp.get("coarsen_alpha",     m_coar_alpha);
  pp.get("max_coarsen_lvl",   m_max_coarsen_level);
}


Vector<Real> air3_zheleznyak_tagger::tracer(const RealVect         a_pos,
					    const Real             a_time,
					    const Real             a_dx,
					    const RealVect         a_E,
					    const Real             a_min_E,
					    const Real             a_max_E,
					    const RealVect         a_grad_E,
					    const Real             a_min_grad_E,
					    const Real             a_max_grad_E){

  const air3_zheleznyak* plaskin = static_cast<air3_zheleznyak*> (&(*m_plaskin));
  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = plaskin->compute_alpha_eff(a_E);


  return tracers;
}
bool air3_zheleznyak_tagger::coarsen_cell(const RealVect         a_pos,
					  const Real             a_time,
					  const Real             a_dx,
					  const int              a_lvl,
					  const Vector<Real>     a_tracer,
					  const Vector<RealVect> a_grad_tracer){

  bool coarsen = false;

  if(a_lvl >= m_max_coarsen_level){
    coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[1]*a_dx < m_coar_alpha;
  }
  else{
    coarsen = false;
  }


  return coarsen;
}


bool air3_zheleznyak_tagger::refine_cell(const RealVect         a_pos,
					 const Real             a_time,
					 const Real             a_dx,
					 const int              a_lvl,
					 const Vector<Real>     a_tracer,
					 const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;

  return refine1 || refine2;
#include "CD_NamespaceFooter.H"
