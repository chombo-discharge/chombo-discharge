/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaCurvatureTagger.cpp
  @brief  Implementation CD_CdrPlasmaCurvatureTagger.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaCurvatureTagger.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaCurvatureTagger::CdrPlasmaCurvatureTagger(){
  m_name = "CdrPlasmaCurvatureTagger";
  m_num_tracers = 1;
}

CdrPlasmaCurvatureTagger::CdrPlasmaCurvatureTagger(const RefCountedPtr<CdrPlasmaPhysics>&     a_plaskin,
						   const RefCountedPtr<CdrPlasmaStepper>&     a_timeStepper,
						   const RefCountedPtr<AmrMesh>&               a_amr,
						   const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry) : CdrPlasmaCurvatureTagger() {
  this->define(a_plaskin, a_timeStepper, a_amr, a_computationalGeometry);
}

CdrPlasmaCurvatureTagger::~CdrPlasmaCurvatureTagger(){

}

void CdrPlasmaCurvatureTagger::parseOptions(){
  parseVerbosity();
  parseBoxes();
  parseBuffer();

  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("coarsen_magnitude", m_coar_mag);
  pp.get("refine_curvature",  m_refi_curv);
  pp.get("refine_magnitude",  m_refi_mag);
}


Vector<Real> CdrPlasmaCurvatureTagger::tracer(const RealVect          a_pos,
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


bool CdrPlasmaCurvatureTagger::coarsenCell(const RealVect         a_pos,
					   const Real             a_time,
					   const Real             a_dx,
					   const int              a_lvl,
					   const Vector<Real>     a_tracer,
					   const Vector<RealVect> a_grad_tracer){

  const bool coarsen = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] < m_coar_curv && a_tracer[0] < m_coar_mag;
  
  return coarsen;
}


bool CdrPlasmaCurvatureTagger::refineCell(const RealVect         a_pos,
					  const Real             a_time,
					  const Real             a_dx,
					  const int              a_lvl,
					  const Vector<Real>     a_tracer,
					  const Vector<RealVect> a_grad_tracer){
  const bool refine = a_grad_tracer[0].vectorLength()*a_dx/a_tracer[0] >= m_refi_curv || a_tracer[0] > m_refi_mag;

  return refine;
}

#include <CD_NamespaceFooter.H>
