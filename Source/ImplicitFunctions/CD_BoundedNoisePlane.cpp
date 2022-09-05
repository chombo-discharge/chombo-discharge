/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BoundedNoisePlane.cpp
  @brief  Declaration of a signed distance function for a noisy plane
  @author Robert Marskar
*/

// Chombo includes
#include <PolyGeom.H>

// Our includes
#include <CD_BoundedNoisePlane.H>
#include <CD_NamespaceHeader.H>

BoundedNoisePlane::BoundedNoisePlane(const int      a_normal,
                                     const RealVect a_point,
                                     const RealVect a_clampLo,
                                     const RealVect a_clampHi,
                                     const Real     a_clampK,
                                     const bool     a_inside,
                                     const Real     a_noiseAmp,
                                     const RealVect a_noiseFreq,
                                     const Real     a_persistence,
                                     const int      a_octaves,
                                     const bool     a_reseed)
{

  // Maximum amplitude that the Perlin noise function can spit out.
  Real amp = 0.0;
  for (int i = 0; i < a_octaves; i++) {
    amp += a_noiseAmp * pow(a_persistence, i);
  }

  // Perlin noise is on [0,1] but we want the scaling on [-0.5, 0.5]
  m_normal              = a_normal;
  const RealVect normal = BASISREALV(a_normal);
  m_point               = a_point - 0.5 * normal * amp;

  m_plane  = RefCountedPtr<BaseIF>(new PlaneIF(normal, m_point, a_inside));
  m_perlin = RefCountedPtr<BaseIF>(new PerlinSdf(a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed));

  m_clampLo = a_clampLo;
  m_clampHi = a_clampHi;
  m_clampK  = a_clampK;
}

BoundedNoisePlane::BoundedNoisePlane(const BoundedNoisePlane& a_inputIF) {}

BoundedNoisePlane::~BoundedNoisePlane() {}

Real
BoundedNoisePlane::value(const RealVect& a_pos) const
{
  // TLDR: To elevate the noise we displace the value along the normal (by an amount given by the Perlin noise function),
  //       clamped with a boxcar function.

  const RealVect n  = BASISREALV(m_normal);
  const RealVect x0 = m_point;
  const RealVect x1 = a_pos;
  const RealVect xp = x1 - PolyGeom::dot((x1 - x0), n) * n;

  auto h = [k = m_clampK](const Real x) { return 1.0 / (1.0 + exp(-2 * k * x)); };

  Real boxCar = 1.0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != m_normal) {
      const Real& x = a_pos[dir];

      boxCar *= h(x - std::min(m_clampLo[dir], m_clampHi[dir])) - h(x - std::max(m_clampLo[dir], m_clampHi[dir]));
    }
  }

  return m_plane->value(a_pos) + m_perlin->value(xp) * boxCar;
}

BaseIF*
BoundedNoisePlane::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new BoundedNoisePlane(*this));
}

#include <CD_NamespaceFooter.H>
