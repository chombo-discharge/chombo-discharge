/*!
  @file   dcel_poly.H
  @brief  Declaration of a two-dimensional polygon class for embedding 3D dcel_face faces. 
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_H_
#define _DCEL_POLY_H_

#include "Vec.H"

#include <memory>
#include <vector>

namespace dcel {
    
  template <class T>
  class Polygon2D {
  public:

    using Vec2 = Vec2T<T>;
    using Vec3 = Vec3T<T>;

    Polygon2D() = delete; 
    Polygon2D(const Vec3& a_normal, const std::vector<Vec3>& a_points);
    ~Polygon2D() = default;

    inline
    decltype(auto) isPointInside(const Vec3& a_point, const InsideOutsideAlgorithm a_algorithm) const noexcept;

    inline
    decltype(auto) isPointInsidePolygonWindingNumber(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) isPointInsidePolygonSubtend(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) isPointInsidePolygonCrossingNumber(const Vec3& a_point) const noexcept;
    
  private:

    inline
    decltype(auto) projectPoint(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) define(const Vec3& a_normal, const std::vector<Vec3>& a_points);

    inline
    decltype(auto) computeWindingNumber(const Vec2& P) const noexcept;

    inline
    decltype(auto) computeCrossingNumber(const Vec2& P) const noexcept;

    inline
    T computeSubtendedAngle(const Vec2& P) const noexcept;


    
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    
    std::vector<Vec2> m_points;
  };
}

#include "dcel_polyI.H"

#endif
