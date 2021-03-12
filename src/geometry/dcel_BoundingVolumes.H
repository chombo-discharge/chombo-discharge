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
    decltype(auto) getPointLocation(const Vec3& a_p) const noexcept;

    inline
    decltype(auto) isPointInside(const Vec3& a_x0) const;

    inline
    decltype(auto) isPointOutside(const Vec3& a_x0) const;
    
    virtual T getDistanceToPoint(const Vec3& a_x0) const = 0;
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
    decltype(auto) define(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept;

    inline
    decltype(auto) intersects(const BoundingSphere& a_other) const noexcept;

    inline
    decltype(auto) getRadius() noexcept;

    inline
    decltype(auto) getRadius() const noexcept;

    inline
    decltype(auto) getCenter() noexcept;

    inline
    decltype(auto) getCenter() const noexcept;

    inline
    T getDistanceToPoint(const Vec3& a_x0) const noexcept override;

  protected:

    T m_radius;
    Vec3 m_center;

    inline
    decltype(auto) buildRitter(const std::vector<Vec3>& a_points) noexcept;
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
    decltype(auto) define(const std::vector<Vec3>& a_points) noexcept;

    inline
    decltype(auto) intersects(const AABB& a_other) const noexcept;

    inline
    T getDistanceToPoint(const Vec3& a_x0) const noexcept override;

    inline
    decltype(auto) getLowCorner() noexcept;

    inline
    decltype(auto) getLowCorner() const noexcept;

    inline
    decltype(auto) getHighCorner() noexcept;

    inline
    decltype(auto) getHighCorner() const noexcept;
    
  protected:

    Vec3 m_loCorner;
    Vec3 m_hiCorner;
  };
}

#include "dcel_BoundingVolumesI.H"

#endif
