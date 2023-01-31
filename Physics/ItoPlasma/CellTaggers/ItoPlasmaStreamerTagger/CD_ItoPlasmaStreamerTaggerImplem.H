/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaStreamerTaggerImplem.H
  @brief  Implementation CD_ItoPlasmaStreamerTagger.H
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaStreamerTaggerImplem_H
#define CD_ItoPlasmaStreamerTaggerImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaStreamerTagger.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

template <typename S>
ItoPlasmaStreamerTagger<S>::ItoPlasmaStreamerTagger()
{
  this->m_name        = "ItoPlasmaStreamerTagger";
  this->m_num_tracers = 2;
}

template <typename S>
ItoPlasmaStreamerTagger<S>::~ItoPlasmaStreamerTagger()
{}

template <typename S>
ItoPlasmaStreamerTagger<S>::ItoPlasmaStreamerTagger(const RefCountedPtr<ItoPlasmaPhysics>&      a_physics,
                                                    const RefCountedPtr<S>&                     a_timeStepper,
                                                    const RefCountedPtr<AmrMesh>&               a_amr,
                                                    const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
  : ItoPlasmaStreamerTagger()
{
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);
}

template <typename S>
void
ItoPlasmaStreamerTagger<S>::parseOptions()
{
  this->parseVerbosity();
  this->parseTagBoxes();
  this->parseRefinementBoxes();
  this->parseBuffer();

  ParmParse pp(this->m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature", m_refi_curv);
  pp.get("refine_alpha", m_refi_alpha);
  pp.get("coarsen_alpha", m_coar_alpha);
  pp.get("max_coarsen_lvl", m_max_coarsen_level);
}

template <typename S>
void
ItoPlasmaStreamerTagger<S>::parseRuntimeOptions()
{
  this->parseVerbosity();
  this->parseTagBoxes();
  this->parseRefinementBoxes();
  this->parseBuffer();

  ParmParse pp(this->m_name.c_str());
  pp.get("coarsen_curvature", m_coar_curv);
  pp.get("refine_curvature", m_refi_curv);
  pp.get("refine_alpha", m_refi_alpha);
  pp.get("coarsen_alpha", m_coar_alpha);
  pp.get("max_coarsen_lvl", m_max_coarsen_level);
}

template <typename S>
Vector<Real>
ItoPlasmaStreamerTagger<S>::tracer(const RealVect a_pos,
                                   const Real     a_time,
                                   const Real     a_dx,
                                   const RealVect a_E,
                                   const Real     a_min_E,
                                   const Real     a_max_E,
                                   const RealVect a_grad_E,
                                   const Real     a_min_grad_E,
                                   const Real     a_max_grad_E)
{

  Vector<Real> tracers(this->m_num_tracers);
  tracers[0] = a_E.vectorLength() / a_max_E;
  tracers[1] = this->m_physics->computeAlpha(a_E);

  return tracers;
}

template <typename S>
bool
ItoPlasmaStreamerTagger<S>::coarsenCell(const RealVect         a_pos,
                                        const Real             a_time,
                                        const Real             a_dx,
                                        const int              a_lvl,
                                        const Vector<Real>     a_tracer,
                                        const Vector<RealVect> a_grad_tracer)
{

  bool coarsen = false;

  if (a_lvl >= m_max_coarsen_level) {
    coarsen = a_grad_tracer[0].vectorLength() * a_dx / a_tracer[0] < m_coar_curv && a_tracer[1] * a_dx < m_coar_alpha;
  }
  else {
    coarsen = false;
  }

  return coarsen;
}

template <typename S>
bool
ItoPlasmaStreamerTagger<S>::refineCell(const RealVect         a_pos,
                                       const Real             a_time,
                                       const Real             a_dx,
                                       const int              a_lvl,
                                       const Vector<Real>     a_tracer,
                                       const Vector<RealVect> a_grad_tracer)
{
  const bool refine1 = a_grad_tracer[0].vectorLength() * a_dx / a_tracer[0] > m_refi_curv;
  const bool refine2 = a_tracer[1] * a_dx > m_refi_alpha;

  return refine1 || refine2;
}

#include <CD_NamespaceFooter.H>

#endif
