/*!
  @file   ml_tagger.cpp
  @brief  Implementation ml_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "ml_tagger.H"
#include "morrow_zheleznyak.H"

#include <ParmParse.H>

ml_tagger::ml_tagger(){
  m_name = "ml_tagger";
  m_num_tracers = 2;
}




ml_tagger::~ml_tagger(){

}

void ml_tagger::parseOptions(){
  parseVerbosity();
  parseBoxes();
  parseBuffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_alpha",      m_refi_alpha);
  pp.get("coarsen_alpha",     m_coar_alpha);
}


Vector<Real> ml_tagger::tracer(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const RealVect         a_E,
			       const Real             a_min_E,
			       const Real             a_max_E,
			       const RealVect         a_grad_E,
			       const Real             a_min_grad_E,
			       const Real             a_max_grad_E){
  
  const morrow_zheleznyak* plaskin = static_cast<morrow_zheleznyak*> (&(*m_plaskin));
  const Real alpha = plaskin->computeAlpha(a_E);
  const Real eta   = plaskin->compute_eta(a_E);
  const Real ve    = plaskin->compute_ve(a_E).vectorLength();

  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = (alpha-eta);


  return tracers;
}
bool ml_tagger::coarsenCell(const RealVect         a_pos,
			     const Real             a_time,
			     const Real             a_dx,
			     const int              a_lvl,
			     const Vector<Real>     a_tracer,
			     const Vector<RealVect> a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[1]*a_dx < m_coar_alpha;
  
  return coarsen;
}


bool ml_tagger::refineCell(const RealVect         a_pos,
			    const Real             a_time,
			    const Real             a_dx,
			    const int              a_lvl,
			    const Vector<Real>     a_tracer,
			    const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;
    
  return refine1 || refine2;
#include "CD_NamespaceFooter.H"
