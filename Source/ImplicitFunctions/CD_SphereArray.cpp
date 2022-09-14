/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
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

// EBGeometry aliases
using AABB   = EBGeometry::BoundingVolumes::AABBT<Real>;
using Vec3   = EBGeometry::Vec3T<Real>;
using SDF    = EBGeometry::SignedDistanceFunction<Real>;
using Sphere = EBGeometry::SphereSDF<Real>;

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
  std::vector<std::shared_ptr<SDF>> spheres;

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

        spheres.emplace_back(std::make_shared<Sphere>(center, a_radius, false));
#if CH_SPACEDIM == 3
      }
#endif
    }
  }

  // Make a fast union. To do this we must have the SDF objects (our vector of
  // spheres) as well as a way for enclosing these objects. We need to define
  // ourselves a lambda that creates an appropriate bounding volumes for each
  // SDF.
  EBGeometry::BVH::BVConstructorT<SDF, AABB> aabbConstructor = [](const std::shared_ptr<const SDF>& a_prim) {
    const Sphere& sph = static_cast<const Sphere&>(*a_prim);

    const Vec3& c = sph.getCenter();
    const Real& r = sph.getRadius();

    const Vec3 lo = c - r * Vec3::one();
    const Vec3 hi = c + r * Vec3::one();

    return AABB(lo, hi);
  };

  // Make the slow and fast unions.
  const auto& constSpheres = m_slowUnion =
    std::make_shared<EBGeometry::Union<Real, Sphere>>((const std::vector<std::shared_ptr<Sphere>>&)spheres, false);
  m_fastUnion = std::make_shared<
    EBGeometry::UnionBVH<Real, Sphere, AABB, SphereArray::K>>((const std::vector<std::shared_ptr<Sphere>>&)spheres,
                                                              false,
                                                              aabbConstructor);
  m_useFast    = a_useFast;
  m_flipInside = a_flipInside;
}

SphereArray::SphereArray(const SphereArray& a_input)
{
  CH_TIME("SphereArray::SphereArray(other)");

  m_slowUnion  = a_input.m_slowUnion;
  m_fastUnion  = a_input.m_fastUnion;
  m_useFast    = a_input.m_useFast;
  m_flipInside = a_input.m_flipInside;
}

SphereArray::~SphereArray() { CH_TIME("SphereArray::~SphereArray"); }

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
