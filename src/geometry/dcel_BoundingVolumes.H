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
    
    virtual T getDistance(const Vec3& a_x0) const = 0;
  };

  template <class T>
  class BoundingSphereT : public BoundingVolumeT<T> {
  public:

    using Vec3 = Vec3T<T>;

    BoundingSphereT();
    BoundingSphereT(const BoundingSphereT& a_other);
    ~BoundingSphereT();
    
    template <class P>
    BoundingSphereT(const std::vector<Vec3T<P> >& a_points, const BoundingVolumeAlgorithm& a_alg = BoundingVolumeAlgorithm::Ritter);

    template <class P>
    inline
    decltype(auto) define(const std::vector<Vec3T<P> >& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept;

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
    decltype(auto) getOverlappingVolume(const BoundingSphereT<T>& a_other) const noexcept;

    inline
    T getDistance(const Vec3& a_x0) const noexcept override;

  protected:

    T m_radius;
    Vec3 m_center;

    template <class P>
    inline
    decltype(auto) buildRitter(const std::vector<Vec3T<P> >& a_points) noexcept;
  };

  template <class T>
  class AABBT : public BoundingVolumeT<T> {
  public:

    using Vec3 = Vec3T<T>;

    AABBT();
    AABBT(const AABBT& a_other);
    AABBT(const std::vector<AABBT>& a_others);
    ~AABBT();
    
    template <class P>
    AABBT(const std::vector<Vec3T<P> >& a_points);

    template <class P>
    inline
    decltype(auto) define(const std::vector<Vec3T<P> >& a_points) noexcept;

    inline
    decltype(auto) intersects(const AABBT& a_other) const noexcept;

    inline
    decltype(auto) getLowCorner() noexcept;

    inline
    decltype(auto) getLowCorner() const noexcept;

    inline
    decltype(auto) getHighCorner() noexcept;

    inline
    decltype(auto) getHighCorner() const noexcept;

    inline
    decltype(auto) getOverlappingVolume(const AABBT<T>& a_other) const noexcept;

    inline
    T getDistance(const Vec3& a_x0) const noexcept override;
    
  protected:

    Vec3 m_loCorner;
    Vec3 m_hiCorner;
  };

  template <class T>
  decltype(auto) intersects(const BoundingSphereT<T>& a_u, const BoundingSphereT<T>& a_v) noexcept;

  template <class T>
  decltype(auto) getOverlappingVolume(const BoundingSphereT<T>& a_u, const BoundingSphereT<T>& a_v) noexcept;

  template <class T>
  decltype(auto) intersects(const AABBT<T>& a_u, const AABBT<T>& a_v) noexcept;

  template <class T>
  decltype(auto) getOverlappingVolume(const AABBT<T>& a_u, const AABBT<T>& a_v) noexcept;


}

#include "dcel_BoundingVolumesI.H"

#endif
