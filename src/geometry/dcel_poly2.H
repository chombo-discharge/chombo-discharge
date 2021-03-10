/*!
  @file   dcel_poly2.H
  @brief  Declaration of a two-dimensional polygon class for embedding 3D dcel_face faces. 
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY2_H_
#define _DCEL_POLY2_H_

#include "RealVect.H"

#include <memory>
#include <vector>

namespace dcel {

  class face;

  class Polygon2D {
  public:
      
    Polygon2D() = delete; 
    Polygon2D(const face& a_face);
    ~Polygon2D() = default;

    bool isPointInside(const RealVect& a_point, const InsideOutsideAlgorithm a_algorithm);

    bool isPointInsidePolygonWindingNumber(const RealVect& a_point) const noexcept;
    bool isPointInsidePolygonSubtend(const RealVect& a_point) const noexcept;
    bool isPointInsidePolygonCrossingNumber(const RealVect& a_point) const noexcept;
    
  protected:

    struct Point2D {
      Point2D() = default;
      Point2D(const Real a_x, const Real a_y);
      ~Point2D() = default;
      Real x;
      Real y;

      Point2D& operator=(const Point2D& a_other) noexcept;
      Point2D& operator+=(const Point2D& a_other) noexcept;
      Point2D& operator-=(const Point2D& a_other) noexcept;
      Point2D operator+(const Point2D& a_other) const noexcept;
      Point2D operator-(const Point2D& a_other) const noexcept;
      Real dotProduct(const Point2D& a_other) noexcept;
      Real length() const noexcept;
      Real length2() const noexcept;
    };

    Point2D projectPoint(const RealVect& a_point) const noexcept;

    bool isPointOnEdge(const Point2D& a_point, const Point2D& a_endPoint1, const Point2D& a_endPoint2, const Real a_thresh) const noexcept;
    bool isPointOnBoundary(const Point2D& a_point, const Real a_thresh) const noexcept;

    void define(const face& a_face) noexcept;

    int computeWindingNumber(const Point2D& P) const noexcept;
    int computeCrossingNumber(const Point2D& P) const noexcept;
    
    Real computeSubtendedAngle(const Point2D& P) const noexcept;
    
  protected:
    
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    std::vector<Point2D> m_points;
  };

}

#endif
