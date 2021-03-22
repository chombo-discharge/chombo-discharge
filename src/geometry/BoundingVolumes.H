/*!
  @file   BoundingVolumes.H
  @brief  Declaration of a various bounding volumes. 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BOUNDINGVOLUMES_H_
#define _BOUNDINGVOLUMES_H_

#include "Vec.H"

#include <vector>

enum class PointLocation{
  Inside,
  Outside,
  Boundary
};

enum class BoundingVolumeAlgorithm {
  Ritter,
};

template <class T>
class BoundingSphereT {
public:

  using Vec3 = Vec3T<T>;

  BoundingSphereT();
  BoundingSphereT(const BoundingSphereT& a_other);
  ~BoundingSphereT();
    
  template <class P>
  BoundingSphereT(const std::vector<Vec3T<P> >& a_points, const BoundingVolumeAlgorithm& a_alg = BoundingVolumeAlgorithm::Ritter);

  template <class P>
  inline
  void define(const std::vector<Vec3T<P> >& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept;

  inline
  bool intersects(const BoundingSphereT& a_other) const noexcept;

  inline
  T& getRadius() noexcept;

  inline
  const T& getRadius() const noexcept;

  inline
  Vec3& getCenter() noexcept;

  inline
  const Vec3& getCenter() const noexcept;

  inline
  T getOverlappingVolume(const BoundingSphereT<T>& a_other) const noexcept;

  inline
  T getDistance(const Vec3& a_x0) const noexcept;

  inline
  T getDistance2(const Vec3& a_x0) const noexcept;

  inline
  T getArea() const noexcept;

protected:

  T m_radius;
  Vec3 m_center;

  template <class P>
  inline
  void buildRitter(const std::vector<Vec3T<P> >& a_points) noexcept;
};

template <class T>
class AABBT {
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
  void define(const std::vector<Vec3T<P> >& a_points) noexcept;

  inline
  bool intersects(const AABBT& a_other) const noexcept;

  inline
  Vec3T<T>& getLowCorner() noexcept;

  inline
  const Vec3T<T>& getLowCorner() const noexcept;

  inline
  Vec3T<T>& getHighCorner() noexcept;

  inline
  const Vec3T<T>& getHighCorner() const noexcept;

  inline
  T getOverlappingVolume(const AABBT<T>& a_other) const noexcept;

  inline
  T getDistance(const Vec3& a_x0) const noexcept;

  inline
  T getDistance2(const Vec3& a_x0) const noexcept;

  inline
  T getArea() const noexcept;
    
protected:

  Vec3 m_loCorner;
  Vec3 m_hiCorner;
};

template <class T>
bool intersects(const BoundingSphereT<T>& a_u, const BoundingSphereT<T>& a_v) noexcept;

template <class T>
T getOverlappingVolume(const BoundingSphereT<T>& a_u, const BoundingSphereT<T>& a_v) noexcept;

template <class T>
bool intersects(const AABBT<T>& a_u, const AABBT<T>& a_v) noexcept;

template <class T>
T getOverlappingVolume(const AABBT<T>& a_u, const AABBT<T>& a_v) noexcept;

#include "BoundingVolumesI.H"

#endif
