/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SphereArray.H
  @brief  Declaration of a sphere array implicit function
  @author Robert Marskar
*/

#ifndef CD_SphereArray_H
#define CD_SphereArray_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <EBGeometry.hpp>
#include <CD_NamespaceHeader.H>

using AABB   = EBGeometry::BoundingVolumes::AABBT<Real>;
using Vec3   = EBGeometry::Vec3T<Real>;
using SDF    = EBGeometry::SignedDistanceFunction<Real>;
using Sphere = EBGeometry::SphereSDF<Real>;

/*!
  @brief A Cartesian array of spheres. 
*/
class SphereArray : public BaseIF
{
public:
  /*!
    @brief Disallowed (for now)
  */
  SphereArray() = delete;

  /*!
    @brief SphereArray full constructor. Constructs an array of spheres
    @param[in] a_radius     Sphere radius
    @param[in] a_lowerLeft  Lower-left sphere center     
    @param[in] a_sphereGap  Distance between spheres
    @param[in] a_numSpheres Number of spheres in each coordinate direction
    @param[in] a_useFast    Use fast union or not
    @param[in] a_flipInside Make outside to inside
    @param[in] a_zCoord     If 2D, this sets the z-coordinate
  */
  SphereArray(const Real     a_radius,
              const RealVect a_loCenter,
              const RealVect a_sphereGap,
              const IntVect  a_numSpheres,
              const bool     a_useFast,
              const bool     a_flipInside,
              const Real     a_zCoord = 0.0);

  /*!
    @brief Copy constructor. Sends the pointers to the object. 
  */
  SphereArray(const SphereArray& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~SphereArray();

  /*!
    @brief Value function
    @param[in] a_point Coordinate
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /*!
    @brief Factory method -- calls the copy constructor
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief Tree degree for BVH accelerator
  */
  static constexpr size_t K = 4;

  /*!
    @brief Switch for using fast or slow union
  */
  bool m_useFast;

  /*!
    @brief Make inside to outside
  */
  bool m_flipInside;

  /*!
    @brief Standard union ala Chombo
  */
  std::shared_ptr<EBGeometry::Union<Real, Sphere>> m_slowUnion;

  /*!
    @brief Union that uses bounding volume hierarchy acceleration
  */
  std::shared_ptr<EBGeometry::UnionBVH<Real, Sphere, AABB, K>> m_fastUnion;
};

#include <CD_NamespaceFooter.H>

#endif
