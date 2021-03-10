/*!
  @file   dcel_poly2.H
  @brief  Declaration of a two-dimensional polygon class for embedding 3D dcel_face faces. 
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY2_H_
#define _DCEL_POLY2_H_

#include "RealVect.H"

#include "dcel_vec.H"

#include <memory>
#include <vector>

namespace dcel {

  class Polygon2D {
  public:
      
    Polygon2D() = delete; 
    Polygon2D(const Vec3<double>& a_normal, const std::vector<Vec3<double> >& a_points);
    ~Polygon2D() = default;

    bool isPointInside(const RealVect& a_point, const InsideOutsideAlgorithm a_algorithm);
    bool isPointInside(const Vec3<double>& a_point, const InsideOutsideAlgorithm a_algorithm);
    
  protected:

    Vec2<double> projectPoint(const Vec3<double>& a_point) const noexcept;

    bool isPointOnEdge(const Vec2<double>& a_point, const Vec2<double>& a_endPoint1, const Vec2<double>& a_endPoint2, const double a_thresh) const noexcept;
    bool isPointOnBoundary(const Vec2<double>& a_point, const double a_thresh) const noexcept;

    void define(const Vec3<double>& a_normal, const std::vector<Vec3<double> >& a_points);

    int computeWindingNumber(const Vec2<double>& P) const noexcept;
    int computeCrossingNumber(const Vec2<double>& P) const noexcept;
    
    double computeSubtendedAngle(const Vec2<double>& P) const noexcept;

    bool isPointInsidePolygonWindingNumber(const Vec3<double>& a_point) const noexcept;
    bool isPointInsidePolygonSubtend(const Vec3<double>& a_point) const noexcept;
    bool isPointInsidePolygonCrossingNumber(const Vec3<double>& a_point) const noexcept;
    
    
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    std::vector<Vec2<double> > m_points;


  };
}



#endif
