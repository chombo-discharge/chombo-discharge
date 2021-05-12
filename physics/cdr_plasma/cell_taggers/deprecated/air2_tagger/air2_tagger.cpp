/*!
  @file   air2_tagger.cpp
  @brief  Implementation air2_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air2_tagger.H"
#include "air2.H"

#include <ParmParse.H>

air2_tagger::air2_tagger(){
  m_name = "air2_tagger";
  m_num_tracers = 2;
}




air2_tagger::~air2_tagger(){

}

void air2_tagger::parseOptions(){
  parseVerbosity();
  parse_boxes();
  parse_buffer();

  ParmParse pp(m_name.c_str());
  Vector<Real> v;
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_alpha",      m_refi_alpha);
  pp.get("coarsen_alpha",     m_coar_alpha);
  pp.get("refine_radius",     m_refine_radius);

  ParmParse p("air2");
  p.getarr("seed_position", v, 0, SpaceDim); m_seed_pos = RealVect(D_DECL(v[0], v[1], v[2]));
  p.get("seed_radius", m_seed_radius);
}


Vector<Real> air2_tagger::tracer(const RealVect         a_pos,
				 const Real             a_time,
				 const Real             a_dx,
				 const RealVect         a_E,
				 const Real             a_min_E,
				 const Real             a_max_E,
				 const RealVect         a_grad_E,
				 const Real             a_min_grad_E,
				 const Real             a_max_grad_E){
  
  const air2* plaskin = static_cast<air2*> (&(*m_plaskin));
  const Real alpha = plaskin->get_alpha(a_E.vectorLength());
  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = alpha;

  return tracers;
}
bool air2_tagger::coarsen_cell(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const int              a_lvl,
			       const Vector<Real>     a_tracer,
			       const Vector<RealVect> a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[1]*a_dx < m_coar_alpha;
  
  return coarsen;
}


bool air2_tagger::refine_cell(const RealVect         a_pos,
			      const Real             a_time,
			      const Real             a_dx,
			      const int              a_lvl,
			      const Vector<Real>     a_tracer,
			      const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;


  const Real dist = (a_pos-m_seed_pos).vectorLength()/m_seed_radius;
  const bool refine3 = dist <= m_refine_radius;

    
  return refine1 || refine2 || refine3;
#include "CD_NamespaceFooter.H"
