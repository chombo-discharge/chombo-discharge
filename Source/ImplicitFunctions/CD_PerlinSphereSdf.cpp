/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSphereSdf.H
  @brief  Implementation of CD_PerlinSphereSdf.H
  @author Robert Marskar
*/

// Our includes
#include <CD_PerlinSphereSdf.H>
#include <CD_NamespaceHeader.H>

PerlinSphereSdf::PerlinSphereSdf(const Real&     a_rad,
                                 const RealVect& a_center,
                                 const bool&     a_inside,
                                 const Real&     a_noiseAmp,
                                 const RealVect& a_noiseFreq,
                                 const Real&     a_persistence,
                                 const int&      a_octaves,
                                 const bool&     a_reseed)
{

  //
  m_rad      = a_rad - a_noiseAmp;
  m_center   = a_center;
  m_inside   = a_inside;
  m_perlinIF = RefCountedPtr<BaseIF>(new PerlinSdf(a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed));
}

PerlinSphereSdf::PerlinSphereSdf(const PerlinSphereSdf& a_inputIF)
{
  m_rad      = a_inputIF.m_rad;
  m_center   = a_inputIF.m_center;
  m_inside   = a_inputIF.m_inside;
  m_perlinIF = a_inputIF.m_perlinIF;
}

PerlinSphereSdf::~PerlinSphereSdf() {}

Real
PerlinSphereSdf::value(const RealVect& a_pos) const
{

  RealVect v;
  Real     retval;

  const RealVect pos = a_pos - m_center;

  // Get noise on the circle/sphere
#if CH_SPACEDIM == 2
  const Real theta = atan2(pos[0], pos[1]);
  const Real x     = m_rad * sin(theta);
  const Real y     = m_rad * cos(theta);

  v = RealVect(x, y);
#elif CH_SPACEDIM == 3
  const Real xy    = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
  const Real theta = atan2(xy, pos[2]);
  const Real phi   = atan2(pos[1], pos[0]);
  const Real x     = m_rad * sin(theta) * sin(phi);
  const Real y     = m_rad * sin(theta) * cos(phi);
  const Real z     = m_rad * cos(theta);

  v = RealVect(x, y, z);
#endif

  // Random radius
  const Real R = m_rad + m_perlinIF->value(v);

  // Value function
  const Real dist2 = pos.vectorLength() * pos.vectorLength() - R * R;
  if (dist2 > 0.) {
    retval = sqrt(dist2);
  }
  else {
    retval = -sqrt(-dist2);
  }

  // Switch inside to outside
  if (!m_inside) {
    retval = -retval;
  }

  return retval;
}

BaseIF*
PerlinSphereSdf::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new PerlinSphereSdf(*this));
}

#include <CD_NamespaceFooter.H>
