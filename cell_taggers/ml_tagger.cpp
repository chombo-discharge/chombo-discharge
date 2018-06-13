/*!
  @file   ml_tagger.cpp
  @brief  Implementation ml_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "ml_tagger.H"
#include "morrow_lowke.H"

#include <ParmParse.H>

ml_tagger::ml_tagger(){
  m_num_tracers = 2;

  m_coar_curv  = 0.1;
  m_refi_curv  = 0.1;
  m_refi_alpha = 1.2;
  m_refi_fudge = 1.0;

  // Get input parameters if we have them
  {
    ParmParse pp("ml_tagger");

    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("refine_alpha",      m_refi_alpha);
    pp.query("refine_fudge",      m_refi_fudge);
  }
}


ml_tagger::~ml_tagger(){

}


Vector<Real> ml_tagger::tracer(const RealVect&         a_pos,
			       const Real&             a_time,
			       const Real&             a_dx,
			       const RealVect&         a_E,
			       const Real&             a_min_E,
			       const Real&             a_max_E,
			       const RealVect&         a_grad_E,
			       const Real&             a_min_grad_E,
			       const Real&             a_max_grad_E){
  
  const morrow_lowke* plaskin = static_cast<morrow_lowke*> (&(*m_plaskin));
  const Real alpha = plaskin->compute_alpha(a_E);
  const Real eta   = plaskin->compute_eta(a_E);
  const Real ve    = plaskin->compute_velocities(a_E)[0].vectorLength();

  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = (eta > 0) ?  (alpha-eta) : 0.0;


  return tracers;
}


bool ml_tagger::coarsen_cell(const RealVect&         a_pos,
			     const Real&             a_time,
			     const Real&             a_dx,
			     const int&              a_lvl,
			     const Vector<Real>&     a_tracer,
			     const Vector<RealVect>& a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv;// && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool ml_tagger::refine_cell(const RealVect&         a_pos,
			    const Real&             a_time,
			    const Real&             a_dx,
			    const int&              a_lvl,
			    const Vector<Real>&     a_tracer,
			    const Vector<RealVect>& a_grad_tracer){
  //  const bool refine1 = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv || a_tracer[0] > m_refi_mag;
  const int max_depth = m_amr->get_max_amr_depth();
  const Real factor   = (max_depth > 0) ? pow(1 - 1./max_depth, a_lvl) : 1;
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2  = a_tracer[1]*a_dx > m_refi_alpha*factor;
    
  if(a_pos[1] < 2E-3){
    return refine1 || refine2;
  }
  else{
    return false;
  }
}
