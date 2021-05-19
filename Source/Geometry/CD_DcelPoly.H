/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelPoly.H
  @brief  Declaration of a two-dimensional polygon class for embedding 3D polygon faces
  @author Robert Marskar
*/

#ifndef CD_DcelPoly_H
#define CD_DcelPoly_H

// Std includes
#include <memory>
#include <vector>

// Our includes
#include <CD_Vec.H>
#include <CD_NamespaceHeader.H>

namespace Dcel {
    
  template <class T>
  class Polygon2D {
  public:

    using Vec2 = Vec2T<T>;
    using Vec3 = Vec3T<T>;

    Polygon2D() = delete; 
    Polygon2D(const Vec3& a_normal, const std::vector<Vec3>& a_points);
    ~Polygon2D() = default;

    inline
    bool isPointInside(const Vec3& a_point, const InsideOutsideAlgorithm a_algorithm) const noexcept;

    inline
    bool isPointInsidePolygonWindingNumber(const Vec3& a_point) const noexcept;

    inline
    bool isPointInsidePolygonSubtend(const Vec3& a_point) const noexcept;

    inline
    bool isPointInsidePolygonCrossingNumber(const Vec3& a_point) const noexcept;
    
  private:

    inline
    Vec2 projectPoint(const Vec3& a_point) const noexcept;

    inline
    void define(const Vec3& a_normal, const std::vector<Vec3>& a_points);

    inline
    int computeWindingNumber(const Vec2& P) const noexcept;

    inline
    int computeCrossingNumber(const Vec2& P) const noexcept;

    inline
    T computeSubtendedAngle(const Vec2& P) const noexcept;


    
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    
    std::vector<Vec2> m_points;
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_DcelPolyImplem.H>

#endif
