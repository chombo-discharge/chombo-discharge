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

  class BoundingVolume {
  public:

    using Vec3 = Vec3T<double>;
    
    BoundingVolume();
    ~BoundingVolume();

    auto getPointLocation(const Vec3& a_p) const noexcept -> PointLocation;
    auto isPointInside(const Vec3& a_x0) const -> bool;
    auto isPointOutside(const Vec3& a_x0) const -> bool;
    
    virtual auto getDistanceToPoint(const Vec3& a_x0) const -> double = 0;
  };
  
  class BoundingSphere : public BoundingVolume {
  public:

    enum class Algorithm {
      Ritter,
    };

    using Vec3 = Vec3T<double>;

    BoundingSphere();
    BoundingSphere(const BoundingSphere& a_other);
    BoundingSphere(const std::vector<Vec3>& a_points, const Algorithm& a_alg);
    ~BoundingSphere();

    auto define(const std::vector<Vec3>& a_points, const Algorithm& a_algorithm) noexcept -> void;

    auto intersects(const BoundingSphere& a_other) const noexcept-> bool;
    
    auto getDistanceToPoint(const Vec3& a_x0) const noexcept -> double override;
    
    auto getRadius() noexcept -> double&;
    auto getRadius() const noexcept -> const double&;
    
    auto getCenter() noexcept -> Vec3&;
    auto getCenter() const noexcept -> const Vec3&;

  protected:

    double m_radius;
    Vec3 m_center;

    auto buildRitter(const std::vector<Vec3>& a_points) noexcept -> void;
  };

  class AABB : public BoundingVolume {
  public:

    using Vec3 = Vec3T<double>;

    AABB();
    AABB(const AABB& a_other);
    AABB(const std::vector<AABB>& a_others);
    AABB(const std::vector<Vec3>& a_points);
    ~AABB();

    auto define(const std::vector<Vec3>& a_points) noexcept -> void;

    auto intersects(const AABB& a_other) const noexcept -> bool;
    
    auto getDistanceToPoint(const Vec3& a_x0) const noexcept -> double override;

    auto getLowCorner() noexcept -> Vec3&;
    auto getLowCorner() const noexcept -> const Vec3&;

    auto getHighCorner() noexcept -> Vec3&;
    auto getHighCorner() const noexcept -> const Vec3&;

    
  protected:

    Vec3 m_loCorner;
    Vec3 m_hiCorner;
  };
}

#endif
