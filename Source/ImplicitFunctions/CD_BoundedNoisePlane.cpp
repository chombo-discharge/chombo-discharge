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

BoundedNoisePlane::BoundedNoisePlane(const std::string a_orientation,
                                     const RealVect    a_point,
                                     const RealVect    a_clampLo,
                                     const RealVect    a_clampHi,
                                     const Real        a_clampK,
                                     const Real        a_noiseAmp,
                                     const RealVect    a_noiseFreq,
                                     const Real        a_persistence,
                                     const int         a_octaves,
                                     const bool        a_reseed)
{

  // Maximum amplitude that the Perlin noise function can spit out.
  m_maxAmp = 0.0;
  for (int i = 0; i < a_octaves; i++) {
    m_maxAmp += a_noiseAmp * pow(a_persistence, i);
  }

  // Set the normal vector
  if (a_orientation == "x+") {
    m_normal = std::make_pair(0, 1);
  }
  else if (a_orientation == "x-") {
    m_normal = std::make_pair(0, -1);
  }
  else if (a_orientation == "y+") {
    m_normal = std::make_pair(1, 1);
  }
  else if (a_orientation == "y-") {
    m_normal = std::make_pair(1, -1);
  }
#if CH_SPACEDIM == 3
  else if (a_orientation == "z+") {
    m_normal = std::make_pair(2, 1);
  }
  else if (a_orientation == "z-") {
    m_normal = std::make_pair(2, -1);
  }
#endif
  else {
    const std::string err = "BoundedNoisePlane::BoundedNoisePlane a_orientation = " + a_orientation + " not supported";
    MayDay::Error(err.c_str());
  }

  const RealVect n = 1.0 * m_normal.second * BASISREALV(m_normal.first);

  m_point = a_point;

  m_plane  = RefCountedPtr<BaseIF>(new PlaneIF(n, m_point, false));
  m_perlin = RefCountedPtr<BaseIF>(new PerlinSdf(a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed));

  m_clampLo = a_clampLo;
  m_clampHi = a_clampHi;
  m_clampK  = a_clampK;
}

BoundedNoisePlane::BoundedNoisePlane(const BoundedNoisePlane& a_inputIF)
{

  m_normal  = a_inputIF.m_normal;
  m_maxAmp  = a_inputIF.m_maxAmp;
  m_point   = a_inputIF.m_point;
  m_plane   = a_inputIF.m_plane;
  m_perlin  = a_inputIF.m_perlin;
  m_clampLo = a_inputIF.m_clampLo;
  m_clampHi = a_inputIF.m_clampHi;
  m_clampK  = a_inputIF.m_clampK;
}

BoundedNoisePlane::~BoundedNoisePlane() {}

Real
BoundedNoisePlane::value(const RealVect& a_pos) const
{
  // TLDR: To elevate the noise we displace the value along the normal (by an amount given by the Perlin noise function),
  //       clamped with a boxcar function.

  const RealVect n  = m_normal.second * BASISREALV(m_normal.first);
  const RealVect x0 = m_point;
  const RealVect x1 = a_pos;
  const RealVect xp = x1 - PolyGeom::dot((x1 - x0), n) * n;

  auto h = [k = m_clampK](const Real x) { return 1.0 / (1.0 + exp(-2 * k * x)); };

  Real boxCar = 1.0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (dir != m_normal.first) {
      const Real& x = xp[dir];

      boxCar *= h(x - std::min(m_clampLo[dir], m_clampHi[dir])) - h(x - std::max(m_clampLo[dir], m_clampHi[dir]));
    }
  }

  return (-m_plane->value(a_pos) + (m_perlin->value(xp) - 0.5 * m_maxAmp) * boxCar);
}

BaseIF*
BoundedNoisePlane::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new BoundedNoisePlane(*this));
}

#include <CD_NamespaceFooter.H>
