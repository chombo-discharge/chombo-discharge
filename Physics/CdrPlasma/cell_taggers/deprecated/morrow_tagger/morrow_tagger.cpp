/*!
  @file   morrow_tagger.cpp
  @brief  Implementation morrow_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "morrow_tagger.H"
#include "morrow_bourdon.H"

#include <ParmParse.H>

morrow_tagger::morrow_tagger(){
  m_name = "morrow_tagger";
  m_num_tracers = 2;
}




morrow_tagger::~morrow_tagger(){

}

void morrow_tagger::parseOptions(){
  parseVerbosity();
  parseBoxes();
  parseBuffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_alpha",      m_refi_alpha);
  pp.get("coarsen_alpha",     m_coar_alpha);
  pp.get("max_coarsen_lvl",   m_max_coarsen_level);
}


Vector<Real> morrow_tagger::tracer(const RealVect         a_pos,
				   const Real             a_time,
				   const Real             a_dx,
				   const RealVect         a_E,
				   const Real             a_min_E,
				   const Real             a_max_E,
				   const RealVect         a_grad_E,
				   const Real             a_min_grad_E,
				   const Real             a_max_grad_E){
  
  const morrow_bourdon* plaskin = static_cast<morrow_bourdon*> (&(*m_plaskin));
  const Real alpha = plaskin->computeAlpha(a_E);
  const Real eta   = plaskin->compute_eta(a_E);

  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = (alpha-eta);


  return tracers;
}
bool morrow_tagger::coarsen_cell(const RealVect         a_pos,
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


bool morrow_tagger::refine_cell(const RealVect         a_pos,
				const Real             a_time,
				const Real             a_dx,
				const int              a_lvl,
				const Vector<Real>     a_tracer,
				const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;
    
  return refine1 || refine2;
#include "CD_NamespaceFooter.H"
