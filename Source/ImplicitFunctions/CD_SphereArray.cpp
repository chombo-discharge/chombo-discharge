/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_SphereArray.cpp
  @brief  Implementation of CD_SphereArray.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <BaseIF.H>

// Our includes
#include <CD_SphereArray.H>
#include <CD_NamespaceHeader.H>

/// @cond DOXYGEN_SKIP
using AABB   = EBGeometry::BoundingVolumes::AABBT<Real>;
using Vec3   = EBGeometry::Vec3T<Real>;
using Sphere = EBGeometry::SphereSDF<Real>;
/// @endcond

constexpr size_t SphereArray::K;

SphereArray::SphereArray(const Real     a_radius,
                         const RealVect a_loCenter,
                         const RealVect a_sphereGap,
                         const IntVect  a_numSpheres,
                         const bool     a_useFast,
                         const bool     a_flipInside,
                         const Real     a_zCoord)
{
  CH_TIME("SphereArray::SphereArray(full)");

  // Make a sphere array
  std::vector<std::shared_ptr<Sphere>> spheres;
  std::vector<AABB>                    boundingVolumes;

  for (size_t i = 0; i < a_numSpheres[0]; i++) {
    for (size_t j = 0; j < a_numSpheres[1]; j++) {
#if CH_SPACEDIM == 3
      for (size_t k = 0; k < a_numSpheres[2]; k++) {
#endif

        const Real x = a_loCenter[0] + i * (a_sphereGap[0] + 2 * a_radius);
        const Real y = a_loCenter[1] + j * (a_sphereGap[1] + 2 * a_radius);
#if CH_SPACEDIM == 2
        const Real z = a_zCoord;
#else
      const Real z = a_loCenter[2] + k * (a_sphereGap[2] + 2 * a_radius);
#endif

        const Vec3 center(x, y, z);

        spheres.emplace_back(std::make_shared<Sphere>(center, a_radius));
        boundingVolumes.emplace_back(center - a_radius * Vec3::one(), center + a_radius * Vec3::one());
#if CH_SPACEDIM == 3
      }
#endif
    }
  }

  // Make the slow and fast unions.
  m_slowUnion  = EBGeometry::Union<Real>(spheres);
  m_fastUnion  = EBGeometry::FastUnion<Real, Sphere, AABB, K>(spheres, boundingVolumes);
  m_useFast    = a_useFast;
  m_flipInside = a_flipInside;
}

SphereArray::SphereArray(const SphereArray& a_input)
  : m_slowUnion(a_input.m_slowUnion),
    m_fastUnion(a_input.m_fastUnion),
    m_flipInside(a_input.m_flipInside),
    m_useFast(a_input.m_useFast)
{
  CH_TIME("SphereArray::SphereArray(other)");
}

SphereArray::~SphereArray()
{
  CH_TIME("SphereArray::~SphereArray");
}

Real
SphereArray::value(const RealVect& a_point) const
{
  CH_TIME("SphereArray::value");

#if CH_SPACEDIM == 2
  const EBGeometry::Vec3T<Real> x(a_point[0], a_point[1], 0.0);
#else
  const EBGeometry::Vec3T<Real> x(a_point[0], a_point[1], a_point[2]);
#endif

  Real dist = 0.0;
  if (m_useFast) {
    dist = m_fastUnion->value(x);
  }
  else {
    dist = m_slowUnion->value(x);
  }

  // Chombo and EBGeometry use opposite sign convetions
  if (!m_flipInside) {
    dist = -dist;
  }

  return dist;
}

BaseIF*
SphereArray::newImplicitFunction() const
{
  CH_TIME("SphereArray::newImplicitFunction");

  return (BaseIF*)(new SphereArray(*this));
}

#include <CD_NamespaceFooter.H>
