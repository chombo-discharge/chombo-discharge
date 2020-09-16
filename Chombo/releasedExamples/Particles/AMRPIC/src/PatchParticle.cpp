#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PatchParticle.H"

#include <iostream>
#include <cmath>

#include "NamespaceHeader.H"

// Flag everything as not defined or set
PatchParticle::PatchParticle()
{
  m_isDefined        = false;
  m_isParametersSet  = false;
  m_isPusherCoeffSet = false;
  m_isMeshInterpSet  = false;

  m_meshInterp       = NULL;
}

PatchParticle::~PatchParticle()
{

  if (m_meshInterp != NULL)
    {
      delete m_meshInterp;
      m_meshInterp = NULL;
    }
}

// define with parameters
void PatchParticle::define(const ProblemDomain&     a_domain,
                           const Real&              a_dx,
                           const InterpType&        a_poissonInterpType,
                           const Cosmology&         a_cosmology)
{
  CH_assert(a_dx > 0.0);

  // Store the domain and grid spacing
  m_domain     = a_domain;
  m_dx         = a_dx;
  m_origin     = RealVect(D_DECL(0.0, 0.0, 0.0));
  m_cosmology  = a_cosmology;

  setParameters(a_poissonInterpType);
  setMeshInterp(a_poissonInterpType);

  m_isDefined = true;
}

/// get particle-mesh object
MeshInterp* PatchParticle::getMeshInterp() const
{
  CH_assert(m_isMeshInterpSet);
  return m_meshInterp;
}

// Set MeshInterp object
void PatchParticle::setMeshInterp(const InterpType& a_poissonInterpType)
{
  CH_assert(m_isParametersSet);
  CH_assert(a_poissonInterpType == m_poissonInterpType);

  if (m_meshInterp != NULL)
    {
      delete m_meshInterp;
    }

  m_meshInterp = static_cast<MeshInterp* > (new MeshInterp(m_domain.domainBox(),
                                                           RealVect(D_DECL(m_dx, m_dx, m_dx)),
                                                           m_origin));
  
  m_isMeshInterpSet = true;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new PatchGodunov object
PatchParticle* PatchParticle::new_patchParticle() const
{
  CH_assert(m_isParametersSet);

  // Make the new object
  PatchParticle* retval =
    static_cast<PatchParticle*>(new PatchParticle());

  // Define the new object
  retval->define(m_domain, m_dx, m_poissonInterpType, m_cosmology);

  // Return the new object
  return retval;
}

void PatchParticle::setParameters(const InterpType& a_poissonInterpType)
{
  m_poissonInterpType = a_poissonInterpType;
  m_isParametersSet = true;
}

// return interpolation scheme used by the patcher
const InterpType& PatchParticle::interpOrder() const
{
  CH_assert(m_isParametersSet);
  return m_poissonInterpType;
}

// Return true if everything is defined and setup
bool PatchParticle::isDefined() const
{
  bool rv = m_isDefined
    && m_isMeshInterpSet
    && m_isParametersSet;

  return rv;
}

///
void PatchParticle::resetPusherCoeff()
{
  m_isPusherCoeffSet = false;
}

// set coefficients for particle pusher;
void PatchParticle::setPusherCoeff(const Real& a_time, const Real& a_dt)
{
  CH_assert(isDefined());

  KDK_PusherCoeff(a_time, a_dt);
  
  m_isPusherCoeffSet = true;
}

// use hybrid (neither peculiar nor comoving v) KDK scheme
void PatchParticle::KDK_PusherCoeff(const Real& a_time, const Real& a_dt)
{
  CH_assert(isDefined());

  const Real a       = m_cosmology.a(a_time);
  const Real ap12    = m_cosmology.a(a_time + 0.5*a_dt);
  const Real ap1     = m_cosmology.a(a_time + 1.0*a_dt);
  const Real adotp12 = m_cosmology.adot(a_time + 0.5*a_dt);

  m_coef_exp_pred = a / ap12;
  m_coef_acc_pred = (1.0 - (a / ap12)) / adotp12;

  m_coef_vel = a_dt / ap12;
  m_coef_exp = ap12 / ap1;
  m_coef_acc = (1.0 - (ap12 / ap1)) / adotp12;

  m_isPusherCoeffSet = true;
}

#include "NamespaceFooter.H"
