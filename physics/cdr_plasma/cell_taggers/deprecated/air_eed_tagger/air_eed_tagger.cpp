/*!
  @file   air_eed_tagger.cpp
  @brief  Implementation air_eed_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "air_eed_tagger.H"
#include "air_eed.H"

#include <ParmParse.H>

air_eed_tagger::air_eed_tagger(){
  m_name = "air_eed_tagger";
  m_num_tracers = 2;
}

air_eed_tagger::~air_eed_tagger(){

}

void air_eed_tagger::parse_options(){
  parse_verbosity();
  parse_boxes();
  parse_buffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_alpha",      m_refi_alpha);
  pp.get("coarsen_alpha",     m_coar_alpha);
  pp.get("max_coarsen_lvl",   m_max_coarsen_level);
}


Vector<Real> air_eed_tagger::tracer(const RealVect         a_pos,
				    const Real             a_time,
				    const Real             a_dx,
				    const RealVect         a_E,
				    const Real             a_min_E,
				    const Real             a_max_E,
				    const RealVect         a_grad_E,
				    const Real             a_min_grad_E,
				    const Real             a_max_grad_E){

  const air_eed* plaskin = static_cast<air_eed*> (&(*m_plaskin));
  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = plaskin->compute_alpha_eff(a_E.vectorLength());


  return tracers;
}
bool air_eed_tagger::coarsen_cell(const RealVect         a_pos,
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


bool air_eed_tagger::refine_cell(const RealVect         a_pos,
				 const Real             a_time,
				 const Real             a_dx,
				 const int              a_lvl,
				 const Vector<Real>     a_tracer,
				 const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;
    
  return refine1 || refine2;
}
