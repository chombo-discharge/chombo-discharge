/*!
  @file   dcel_BoundingVolumes.H
  @brief  Declaration of a various bounding volumes. 
  @author Robert Marskar
  @date   March 2021
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
  class BoundingVolumeT {
  public:

    using Vec3 = Vec3T<T>;
    
    BoundingVolumeT();
    ~BoundingVolumeT();

    inline
    decltype(auto) getPointLocation(const Vec3& a_p) const noexcept;

    inline
    decltype(auto) isPointInside(const Vec3& a_x0) const;

    inline
    decltype(auto) isPointOutside(const Vec3& a_x0) const;
    
    virtual T getSignedMinDistanceToPoint(const Vec3& a_x0) const = 0;

    virtual T getUnsignedMaxDistanceToPoint(const Vec3& a_x0) const = 0;
  };

  template <class T>
  class BoundingSphereT : public BoundingVolumeT<T> {
  public:

    using Vec3 = Vec3T<T>;

    BoundingSphereT();
    BoundingSphereT(const BoundingSphereT& a_other);
    BoundingSphereT(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_alg);
    ~BoundingSphereT();

    inline
    decltype(auto) define(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept;

    inline
    decltype(auto) intersects(const BoundingSphereT& a_other) const noexcept;

    inline
    decltype(auto) getRadius() noexcept;

    inline
    decltype(auto) getRadius() const noexcept;

    inline
    decltype(auto) getCenter() noexcept;

    inline
    decltype(auto) getCenter() const noexcept;

    inline
    T getSignedMinDistanceToPoint(const Vec3& a_x0) const noexcept override;

    inline
    T getUnsignedMaxDistanceToPoint(const Vec3& a_x0) const noexcept override;

  protected:

    T m_radius;
    Vec3 m_center;

    inline
    decltype(auto) buildRitter(const std::vector<Vec3>& a_points) noexcept;
  };

  template <class T>
  class AABBT : public BoundingVolumeT<T> {
  public:

    using Vec3 = Vec3T<T>;

    AABBT();
    AABBT(const AABBT& a_other);
    AABBT(const std::vector<AABBT>& a_others);
    AABBT(const std::vector<Vec3>& a_points);
    ~AABBT();

    inline
    decltype(auto) define(const std::vector<Vec3>& a_points) noexcept;

    inline
    decltype(auto) intersects(const AABBT& a_other) const noexcept;

    inline
    T getSignedMinDistanceToPoint(const Vec3& a_x0) const noexcept override;

    inline
    T getUnsignedMaxDistanceToPoint(const Vec3& a_x0) const noexcept override;

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
