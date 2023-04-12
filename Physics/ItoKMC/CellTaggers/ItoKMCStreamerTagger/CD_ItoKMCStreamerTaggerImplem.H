/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCStreamerTaggerImplem.H
  @brief  Implementation CD_ItoKMCStreamerTagger.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCStreamerTaggerImplem_H
#define CD_ItoKMCStreamerTaggerImplem_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCStreamerTagger.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

template <typename S>
ItoKMCStreamerTagger<S>::ItoKMCStreamerTagger() noexcept
{
  CH_TIME("ItoKMCStreamerTagger::ItoKMCStreamerTagger");

  this->m_name         = "ItoKMCStreamerTagger";
  this->m_numTagFields = 2;
}

template <typename S>
ItoKMCStreamerTagger<S>::~ItoKMCStreamerTagger() noexcept
{}

template <typename S>
ItoKMCStreamerTagger<S>::ItoKMCStreamerTagger(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                                              const RefCountedPtr<S>&             a_timeStepper,
                                              const RefCountedPtr<AmrMesh>&       a_amr) noexcept
  : ItoKMCStreamerTagger<S>()
{
  CH_TIME("ItoKMCStreamerTagger::ItoKMCStreamerTagger");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::ItoKMCStreamerTagger" << endl;
  }

  this->define(a_physics, a_timeStepper, a_amr);
}

template <typename S>
void
ItoKMCStreamerTagger<S>::parseOptions() noexcept
{
  CH_TIME("ItoKMCStreamerTagger::parseOptions");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseTagBoxes();
  this->parseRefinementBoxes();
  this->parseBuffer();

  ParmParse pp(this->m_name.c_str());
  pp.get("coarsen_curvature", this->m_coarsenCurvature);
  pp.get("refine_curvature", this->m_refineCurvature);
  pp.get("refine_alpha", this->m_refineAlpha);
  pp.get("coarsen_alpha", this->m_coarsenAlpha);
  pp.get("max_coarsen_lvl", this->m_maxCoarsenLevel);
}

template <typename S>
void
ItoKMCStreamerTagger<S>::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCStreamerTagger::parseOptions");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseTagBoxes();
  this->parseRefinementBoxes();
  this->parseBuffer();

  ParmParse pp(this->m_name.c_str());
  pp.get("coarsen_curvature", this->m_coarsenCurvature);
  pp.get("refine_curvature", this->m_refineCurvature);
  pp.get("refine_alpha", this->m_refineAlpha);
  pp.get("coarsen_alpha", this->m_coarsenAlpha);
  pp.get("max_coarsen_lvl", this->m_maxCoarsenLevel);
}

template <typename S>
Vector<Real>
ItoKMCStreamerTagger<S>::computeTagFields(const RealVect a_pos,
                                          const Real     a_time,
                                          const Real     a_dx,
                                          const RealVect a_E,
                                          const Real     a_minE,
                                          const Real     a_maxE,
                                          const RealVect a_gradE,
                                          const Real     a_minGradE,
                                          const Real     a_maxGradE) const noexcept
{
  CH_TIME("ItoKMCStreamerTagger::computeTagFields");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::computeTagFields" << endl;
  }

  Vector<Real> tagFields(this->m_numTagFields);

  tagFields[0] = a_E.vectorLength() / a_maxE;
  tagFields[1] = this->m_physics->computeAlpha(a_E) * a_dx;

  return tagFields;
}

template <typename S>
bool
ItoKMCStreamerTagger<S>::coarsenCell(const RealVect         a_pos,
                                     const Real             a_time,
                                     const Real             a_dx,
                                     const int              a_lvl,
                                     const Vector<Real>     a_tagFields,
                                     const Vector<RealVect> a_gradTagFields) const noexcept
{
  CH_TIME("ItoKMCStreamerTagger::coarsenCell");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::coarsenCell" << endl;
  }

  bool coarsen = false;

  if (a_lvl >= m_maxCoarsenLevel) {
    const Real     E       = a_tagFields[0];
    const Real     alphaDx = a_tagFields[1];
    const RealVect gradE   = a_gradTagFields[0];

    const bool coarsenCurv  = gradE.vectorLength() * a_dx / E < m_coarsenCurvature;
    const bool coarsenAlpha = alphaDx < m_coarsenAlpha;

    coarsen = coarsenCurv && coarsenAlpha;
  }
  else {
    coarsen = false;
  }

  return coarsen;
}

template <typename S>
bool
ItoKMCStreamerTagger<S>::refineCell(const RealVect         a_pos,
                                    const Real             a_time,
                                    const Real             a_dx,
                                    const int              a_lvl,
                                    const Vector<Real>     a_tagFields,
                                    const Vector<RealVect> a_gradTagFields) const noexcept
{
  CH_TIME("ItoKMCStreamerTagger::refineCell");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCStreamerTagger::refineCell" << endl;
  }

  const Real     E       = a_tagFields[0];
  const Real     alphaDx = a_tagFields[1];
  const RealVect gradE   = a_gradTagFields[0];

  const bool refineCurv  = gradE.vectorLength() * a_dx / E > m_refineCurvature;
  const bool refineAlpha = alphaDx > m_refineAlpha;

  return refineCurv || refineAlpha;
}

#include <CD_NamespaceFooter.H>

#endif
