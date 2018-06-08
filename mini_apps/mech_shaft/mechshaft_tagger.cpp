/*!
  @file   mechshaft_tagger.cpp
  @brief  Implementation mechshaft_tagger.cpp
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "mechshaft_tagger.H"

#include <ParmParse.H>
#include <PolyGeom.H>

mechshaft_tagger::mechshaft_tagger(){
  m_num_tracers = 1;

  m_coar_curv  = 0.1;
  m_coar_mag   = 0.1;
  m_refi_curv  = 0.2;
  m_refi_mag   = 0.75;
  m_lo         = -RealVect::Unit;
  m_hi         =  RealVect::Unit;

  // Get input parameters if we have them
  {
    ParmParse pp("mechshaft_tagger");

    pp.query("coarsen_curvature", m_coar_curv);
    pp.query("coarsen_magnitude", m_coar_mag);
    pp.query("refine_curvature",  m_refi_curv);
    pp.query("refine_magnitude",  m_refi_mag);
  }

  { // Geometrical inputs
    Vector<Real> vec1, vec2;
    ParmParse pp("mechshaft_tagger");
    pp.queryarr("tagbox_lo", vec1, 0, SpaceDim);
    pp.queryarr("tagbox_hi", vec2, 0, SpaceDim);
    
    const RealVect m_lo = RealVect(D_DECL(vec1[0], vec1[1], vec1[2]));
    const RealVect m_hi = RealVect(D_DECL(vec2[0], vec2[1], vec2[2]));
  }
}


mechshaft_tagger::~mechshaft_tagger(){

}


Vector<Real> mechshaft_tagger::tracer(const RealVect&         a_pos,
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


bool mechshaft_tagger::coarsen_cell(const RealVect&         a_pos,
				    const Real&             a_time,
				    const Real&             a_dx,
				    const int&              a_lvl,
				    const Vector<Real>&     a_tracer,
				    const Vector<RealVect>& a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool mechshaft_tagger::refine_cell(const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dx,
				   const int&              a_lvl,
				   const Vector<Real>&     a_tracer,
				   const Vector<RealVect>& a_grad_tracer){
  
  const bool refine_tag      = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv || a_tracer[0] > m_refi_mag;
  const bool refine_geom     = a_pos > m_lo && a_pos < m_hi;

  return (refine_tag && refine_geom);
}
