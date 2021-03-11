/*!
  @file   dcel_BoundingSphere.H
  @brief  Declaration of a BoundingSphere class. Sphere can be computed with various approximate or exact functions. 
  @author Robert Marskar
  @date   Ap. 2018
*/

#ifndef _DCEL_BOUNDINGVOLUMES_H_
#define _DCEL_BOUNDINGVOLUMES_H_

#include "dcel_vec.H"

#include <vector>

namespace dcel {

  enum class PointLocation{
    Inside,
    Outside,
    Boundary
  };

  enum class BoundingVolumeAlgorithm {
    Ritter,
  };

  template <class T>
  class BoundingVolume {
  public:

    using Vec3 = Vec3T<T>;
    
    BoundingVolume();
    ~BoundingVolume();

    inline
    auto getPointLocation(const Vec3& a_p) const noexcept -> PointLocation;

    inline
    auto isPointInside(const Vec3& a_x0) const -> bool;

    inline
    auto isPointOutside(const Vec3& a_x0) const -> bool;
    
    virtual auto getDistanceToPoint(const Vec3& a_x0) const -> T = 0;
  };

  template <class T>
  class BoundingSphere : public BoundingVolume<T> {
  public:

    using Vec3 = Vec3T<T>;

    BoundingSphere();
    BoundingSphere(const BoundingSphere& a_other);
    BoundingSphere(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_alg);
    ~BoundingSphere();

    inline
    auto define(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept -> void;

    inline
    auto intersects(const BoundingSphere& a_other) const noexcept-> bool;

    inline
    auto getDistanceToPoint(const Vec3& a_x0) const noexcept -> T override;

    inline
    auto getRadius() noexcept -> T&;

    inline
    auto getRadius() const noexcept -> const T&;

    inline
    auto getCenter() noexcept -> Vec3&;

    inline
    auto getCenter() const noexcept -> const Vec3&;

  protected:

    T m_radius;
    Vec3 m_center;

    inline
    auto buildRitter(const std::vector<Vec3>& a_points) noexcept -> void;
  };

  template <class T>
  class AABB : public BoundingVolume<T> {
  public:

    using Vec3 = Vec3T<T>;

    AABB();
    AABB(const AABB& a_other);
    AABB(const std::vector<AABB>& a_others);
    AABB(const std::vector<Vec3>& a_points);
    ~AABB();

    inline
    auto define(const std::vector<Vec3>& a_points) noexcept -> void;

    inline
    auto intersects(const AABB& a_other) const noexcept -> bool;

    inline
    auto getDistanceToPoint(const Vec3& a_x0) const noexcept -> T override;

    inline
    auto getLowCorner() noexcept -> Vec3&;

    inline
    auto getLowCorner() const noexcept -> const Vec3&;

    inline
    auto getHighCorner() noexcept -> Vec3&;

    inline
    auto getHighCorner() const noexcept -> const Vec3&;
    
  protected:

    Vec3 m_loCorner;
    Vec3 m_hiCorner;
  };
}

#include "dcel_BoundingVolumesI.H"

#endif
