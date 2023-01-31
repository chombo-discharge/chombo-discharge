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
ItoPlasmaStreamerTagger<S>::ItoPlasmaStreamerTagger() noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::ItoPlasmaStreamerTagger");

  this->m_name         = "ItoPlasmaStreamerTagger";
  this->m_numTagFields = 2;
}

template <typename S>
ItoPlasmaStreamerTagger<S>::~ItoPlasmaStreamerTagger() noexcept
{}

template <typename S>
ItoPlasmaStreamerTagger<S>::ItoPlasmaStreamerTagger(const RefCountedPtr<ItoPlasmaPhysics>& a_physics,
                                                    const RefCountedPtr<S>&                a_timeStepper,
                                                    const RefCountedPtr<AmrMesh>&          a_amr) noexcept
  : ItoPlasmaStreamerTagger<S>()
{
  CH_TIME("ItoPlasmaStreamerTagger::ItoPlasmaStreamerTagger");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::ItoPlasmaStreamerTagger" << endl;
  }

  this->define(a_physics, a_timeStepper, a_amr);
}

template <typename S>
void
ItoPlasmaStreamerTagger<S>::parseOptions() noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::parseOptions");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::parseOptions" << endl;
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
ItoPlasmaStreamerTagger<S>::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::parseOptions");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::parseOptions" << endl;
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
ItoPlasmaStreamerTagger<S>::computeTagFields(const RealVect a_pos,
                                             const Real     a_time,
                                             const Real     a_dx,
                                             const RealVect a_E,
                                             const Real     a_minE,
                                             const Real     a_maxE,
                                             const RealVect a_gradE,
                                             const Real     a_minGradE,
                                             const Real     a_maxGradE) const noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::computeTagFields");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::computeTagFields" << endl;
  }

  Vector<Real> tagFields(this->m_numTagFields);

  tagFields[0] = a_E.vectorLength() / a_maxE;
  tagFields[1] = this->m_physics->computeAlpha(a_E) * a_dx;

  return tagFields;
}

template <typename S>
bool
ItoPlasmaStreamerTagger<S>::coarsenCell(const RealVect         a_pos,
                                        const Real             a_time,
                                        const Real             a_dx,
                                        const int              a_lvl,
                                        const Vector<Real>     a_tagFields,
                                        const Vector<RealVect> a_gradTagFields) const noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::coarsenCell");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::coarsenCell" << endl;
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
ItoPlasmaStreamerTagger<S>::refineCell(const RealVect         a_pos,
                                       const Real             a_time,
                                       const Real             a_dx,
                                       const int              a_lvl,
                                       const Vector<Real>     a_tagFields,
                                       const Vector<RealVect> a_gradTagFields) const noexcept
{
  CH_TIME("ItoPlasmaStreamerTagger::refineCell");
  if (this->m_verbosity > 5) {
    pout() << "ItoPlasmaStreamerTagger::refineCell" << endl;
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
