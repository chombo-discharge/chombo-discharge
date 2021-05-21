/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrPlasmaStreamerTagger.cpp
  @brief  Implementation CD_CdrPlasmaStreamerTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaStreamerTagger.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(){
  m_name = "CdrPlasmaStreamerTagger";

  m_num_tracers = 2;
}

CdrPlasmaStreamerTagger::~CdrPlasmaStreamerTagger(){

}

CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(const RefCountedPtr<CdrPlasmaPhysics>&     a_physics,
						 const RefCountedPtr<cdr_plasma_stepper>&     a_timeStepper,
						 const RefCountedPtr<AmrMesh>&               a_amr,
						 const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry) : CdrPlasmaStreamerTagger() {
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);
}

void CdrPlasmaStreamerTagger::parseOptions(){
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


Vector<Real> CdrPlasmaStreamerTagger::tracer(const RealVect         a_pos,
					     const Real             a_time,
					     const Real             a_dx,
					     const RealVect         a_E,
					     const Real             a_min_E,
					     const Real             a_max_E,
					     const RealVect         a_grad_E,
					     const Real             a_min_grad_E,
					     const Real             a_max_grad_E){


  
  Vector<Real> tracers(m_num_tracers);
  tracers[0] = a_E.vectorLength()/a_max_E;
  tracers[1] = m_physics->computeAlpha(a_E);


  return tracers;
}
bool CdrPlasmaStreamerTagger::coarsenCell(const RealVect         a_pos,
					  const Real             a_time,
					  const Real             a_dx,
					  const int              a_lvl,
					  const Vector<Real>     a_tracer,
					  const Vector<RealVect> a_grad_tracer){

  bool coarsen = false;

  if(a_lvl >= m_max_coarsen_level){
    coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv || a_tracer[1]*a_dx < m_coar_alpha;
  }
  else{
    coarsen = false;
  }


  return coarsen;
}


bool CdrPlasmaStreamerTagger::refineCell(const RealVect         a_pos,
					 const Real             a_time,
					 const Real             a_dx,
					 const int              a_lvl,
					 const Vector<Real>     a_tracer,
					 const Vector<RealVect> a_grad_tracer){
  const bool refine1  = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1]*a_dx > m_refi_alpha;

  return refine1 || refine2;
}

#include <CD_NamespaceFooter.H>
