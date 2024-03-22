/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinPlaneSdf.cpp
  @brief  CD_Implementation of CD_PerlinPlaneSdf.H
  @author Robert Marskar
*/

// Chombo includes
#include <PlaneIF.H>
#include <PolyGeom.H>

// Our includes
#include <CD_PerlinSdf.H>
#include <CD_PerlinPlaneSdf.H>
#include <CD_NamespaceHeader.H>

PerlinPlaneSdf::PerlinPlaneSdf(const RealVect a_normal,
                               const RealVect a_point,
                               const bool     a_inside,
                               const Real     a_noiseAmp,
                               const RealVect a_noiseFreq,
                               const Real     a_persistence,
                               const int      a_octaves,
                               const bool     a_reseed)
{
  // This is the maximum noise the Perlin will spit out.
  Real amp = 0.0;
  for (int i = 0; i < a_octaves; i++) {
    amp += a_noiseAmp * pow(a_persistence, i);
  }

  // Adjust point so noise amplitude does not affect the average position of the plane
  m_point  = a_point - 0.5 * a_normal * amp; // Perlin noise is on [0,1] but we want the scaling on [-0.5, 0.5]
  m_normal = a_normal;

  m_plane  = RefCountedPtr<BaseIF>(new PlaneIF(m_normal, m_point, a_inside));
  m_perlin = RefCountedPtr<BaseIF>(new PerlinSdf(a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed));
}

PerlinPlaneSdf::PerlinPlaneSdf(const PerlinPlaneSdf& a_inputIF)
{
  m_normal = a_inputIF.m_normal;
  m_point  = a_inputIF.m_point;
  m_plane  = a_inputIF.m_plane;
  m_perlin = a_inputIF.m_perlin;
}

PerlinPlaneSdf::~PerlinPlaneSdf()
{}

Real
PerlinPlaneSdf::value(const RealVect& a_pos) const
{
  // TLDR: To elevate the noise we displace the value along the normal (by an amount given by the Perlin noise function).

  const RealVect x0 = m_point;
  const RealVect x1 = a_pos;
  const RealVect xp = x1 - PolyGeom::dot((x1 - x0), m_normal) * m_normal;

  return m_plane->value(a_pos) + m_perlin->value(xp);
}

BaseIF*
PerlinPlaneSdf::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new PerlinPlaneSdf(*this));
}

#include <CD_NamespaceFooter.H>
