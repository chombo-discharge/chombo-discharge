/*!
  @file   field_tagger.cpp
  @brief  Implementation field_tagger.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "field_tagger.H"

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace physics::cdr_plasma;

field_tagger::field_tagger(){
  m_name = "field_tagger";
  m_num_tracers = 1;
}

field_tagger::field_tagger(const RefCountedPtr<cdr_plasma_physics>&     a_plaskin,
			   const RefCountedPtr<cdr_plasma_stepper>&     a_timeStepper,
			   const RefCountedPtr<AmrMesh>&               a_amr,
			   const RefCountedPtr<computational_geometry>& a_computationalGeometry) : field_tagger() {
  this->define(a_plaskin, a_timeStepper, a_amr, a_computationalGeometry);
}

field_tagger::~field_tagger(){

}

void field_tagger::parseOptions(){
  parseVerbosity();
  parseBoxes();
  parseBuffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("coarsen_magnitude", m_coar_mag);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_magnitude",  m_refi_mag);
}


Vector<Real> field_tagger::tracer(const RealVect          a_pos,
				  const Real              a_time,
				  const Real              a_dx,
				  const RealVect          a_E,
				  const Real              a_min_E,
				  const Real              a_max_E,
				  const RealVect          a_grad_E,
				  const Real              a_min_grad_E,
				  const Real              a_max_grad_E){

  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;

  return tracers;
}


bool field_tagger::coarsen_cell(const RealVect         a_pos,
				const Real             a_time,
				const Real             a_dx,
				const int              a_lvl,
				const Vector<Real>     a_tracer,
				const Vector<RealVect> a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool field_tagger::refine_cell(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const int              a_lvl,
			       const Vector<Real>     a_tracer,
			       const Vector<RealVect> a_grad_tracer){
  const bool refine = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] >= m_refi_curv || a_tracer[0] > m_refi_mag;

  return refine;
}
#include "CD_NamespaceFooter.H"
