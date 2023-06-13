/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

CdrPlasmaStreamerTagger::~CdrPlasmaStreamerTagger()
{
  CH_TIME("CdrPlasmaStreamerTagger::~CdrPlasmaStreamerTagger()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaStreamerTagger::~CdrPlasmaStreamerTagger()" << endl;
  }
}

CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                                                 const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                                                 const RefCountedPtr<AmrMesh>&               a_amr,
                                                 const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(RefCountedPtr<...>x4)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaStreamerTagger::CdrPlasmaStreamerTagger(RefCountedPtr<...>x4)" << endl;
  }

  // Call parent define function for setting the input parameters.
  this->define(a_physics, a_timeStepper, a_amr, a_computationalGeometry);

  m_name       = "CdrPlasmaStreamerTagger";
  m_numTracers = 2; // 2 tracer fields -- the electric field and the Townsend ionization coefficient.
}

void
CdrPlasmaStreamerTagger::parseOptions()
{
  CH_TIME("CdrPlasmaStreamerTagger::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaStreamerTagger::parseOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseTagBoxes();
  this->parseRefinementBoxes();
  this->parseBuffer();

  // Parse class options.
  ParmParse pp(m_name.c_str());
  pp.get("coarsen_curvature", m_coarCurv);
  pp.get("refine_curvature", m_refiCurv);
  pp.get("refine_alpha", m_refiAlpha);
  pp.get("coarsen_alpha", m_coarAlpha);
  pp.get("max_coarsen_lvl", m_maxCoarsenLevel);
}

void
CdrPlasmaStreamerTagger::parseRuntimeOptions()
{
  CH_TIME("CdrPlasmaStreamerTagger::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaStreamerTagger::parseRuntimeOptions()" << endl;
  }
  this->parseOptions();
}

Vector<Real>
CdrPlasmaStreamerTagger::tracer(const RealVect a_pos,
                                const Real     a_time,
                                const Real     a_dx,
                                const RealVect a_electricField,
                                const Real     a_minElectricField,
                                const Real     a_maxElectricField,
                                const RealVect a_gradElectricField,
                                const Real     a_minGradElectricField,
                                const Real     a_maxGradElectricField) const
{
  Vector<Real> tracers(m_numTracers, 0.0);

  // We define the two tracer fields as the electric strength |E| and as the alpha-coefficient.

  const Real E = a_electricField.vectorLength();

  // Set tracer fields
  tracers[0] = E / a_maxElectricField;
  tracers[1] = m_physics->computeAlpha(E, a_pos);
  tracers[1] -= m_physics->computeEta(E, a_pos);

  return tracers;
}

bool
CdrPlasmaStreamerTagger::coarsenCell(const RealVect         a_pos,
                                     const Real             a_time,
                                     const Real             a_dx,
                                     const int              a_lvl,
                                     const Vector<Real>     a_tracers,
                                     const Vector<RealVect> a_gradTracers) const
{
  bool coarsen = false;

  // TLDR: Coarsen if both criteria are met.

  if (a_lvl >= m_maxCoarsenLevel) {
    coarsen = a_gradTracers[0].vectorLength() * a_dx / a_tracers[0] < m_coarCurv && a_tracers[1] * a_dx < m_coarAlpha;
  }
  else {
    coarsen = false;
  }

  return coarsen;
}

bool
CdrPlasmaStreamerTagger::refineCell(const RealVect         a_pos,
                                    const Real             a_time,
                                    const Real             a_dx,
                                    const int              a_lvl,
                                    const Vector<Real>     a_tracers,
                                    const Vector<RealVect> a_gradTracers) const
{
  // TLDR: Refine if either criterion are met.

  const bool refine1 = a_gradTracers[0].vectorLength() * a_dx / a_tracers[0] > m_refiCurv;
  const bool refine2 = a_tracers[1] * a_dx > m_refiAlpha;
  const bool refine3 = a_lvl < this->getManualRefinementLevel(a_pos);

  return refine1 || refine2 || refine3;
}

#include <CD_NamespaceFooter.H>
