/*!
  @file   field_tagger.cpp
  @brief  Implementation field_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "field_tagger.H"

#include <ParmParse.H>

field_tagger::field_tagger(){
  m_num_tracers = 1;

  m_coar_curv = 0.1;
  m_coar_mag  = 0.1;
  m_refi_curv = 0.2;
  m_refi_mag  = 0.75;

  // Get input parameters if we have them
  {
    ParmParse pp("field_tagger");

    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("coarsen_magnitude", m_coar_mag);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("refine_magnitude",  m_refi_mag);
  }
}


field_tagger::~field_tagger(){

}


Vector<Real> field_tagger::tracer(const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dx,
				  const RealVect&         a_E,
				  const Real&             a_min_E,
				  const Real&             a_max_E,
				  const RealVect&         a_grad_E,
				  const Real&             a_min_grad_E,
				  const Real&             a_max_grad_E){

  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;

  return tracers;
}


bool field_tagger::coarsen_cell(const RealVect&         a_pos,
				const Real&             a_time,
				const Real&             a_dx,
				const int&              a_lvl,
				const Vector<Real>&     a_tracer,
				const Vector<RealVect>& a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool field_tagger::refine_cell(const RealVect&         a_pos,
			       const Real&             a_time,
			       const Real&             a_dx,
			       const int&              a_lvl,
			       const Vector<Real>&     a_tracer,
			       const Vector<RealVect>& a_grad_tracer){
  const bool refine = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] >= m_refi_curv || a_tracer[0] > m_refi_mag;

  return refine;
}
