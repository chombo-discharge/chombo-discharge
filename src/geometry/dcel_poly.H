/*!
  @file   dcel_poly2.H
  @brief  Declaration of a two-dimensional polygon class for embedding 3D dcel_face faces. 
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY2_H_
#define _DCEL_POLY2_H_

#include "dcel_vec.H"

#include <memory>
#include <vector>

namespace dcel {

  template <class T>
  class Polygon2D {
  public:
      
    Polygon2D() = delete; 
    Polygon2D(const Vec3<T>& a_normal, const std::vector<Vec3<T> >& a_points);
    ~Polygon2D() = default;

    inline
    bool isPointInside(const Vec3<T>& a_point, const InsideOutsideAlgorithm a_algorithm);
    
  private:

    inline
    Vec2<T> projectPoint(const Vec3<T>& a_point) const noexcept;

    inline
    void define(const Vec3<T>& a_normal, const std::vector<Vec3<T> >& a_points);

    inline
    int computeWindingNumber(const Vec2<T>& P) const noexcept;

    inline
    int computeCrossingNumber(const Vec2<T>& P) const noexcept;

    inline
    T computeSubtendedAngle(const Vec2<T>& P) const noexcept;

    inline
    bool isPointInsidePolygonWindingNumber(const Vec3<T>& a_point) const noexcept;

    inline
    bool isPointInsidePolygonSubtend(const Vec3<T>& a_point) const noexcept;

    inline
    bool isPointInsidePolygonCrossingNumber(const Vec3<T>& a_point) const noexcept;
    
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    
    std::vector<Vec2<T> > m_points;
  };
}

#include "dcel_polyI.H"

#endif
