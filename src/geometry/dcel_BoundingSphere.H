/*!
  @file   dcel_BoundingSphere.H
  @brief  Declaration of a BoundingSphere class. Sphere can be computed with various approximate or exact functions. 
  @author Robert Marskar
  @date   Ap. 2018
*/

#ifndef _DCEL_RITTER_
#define _DCEL_RITTER_

#include "RealVect.H"

#include <vector>

namespace dcel {

  enum class PointLocation{
    Inside,
    Outside,
    Boundary
  };
  
  /*!
    @brief Class for computing a bounding sphere for a set of points by using Ritter's method
  */
  class BoundingSphere {
  public:

    enum class Algorithm {
      Ritter,
    };

    BoundingSphere();
    BoundingSphere(const BoundingSphere& a_other);
    BoundingSphere(const std::vector<RealVect>& a_points, const Algorithm& a_alg);
    ~BoundingSphere();

    void define(const std::vector<RealVect>& a_points, const Algorithm& a_algorith);

    PointLocation getPointLocation(const RealVect& a_p) const;
    
    bool isPointInside(const RealVect& a_x0) const;
    bool isPointOutside(const RealVect& a_x0) const;
    
    Real getDistanceToPoint(const RealVect& a_x0) const;
    
    Real& getRadius();
    const Real& getRadius() const;
    
    RealVect& getCenter();
    const RealVect& getCenter() const;

  protected:

    Real m_radius;
    RealVect m_center;

    void buildRitter(const std::vector<RealVect>& a_points);
  };
}
#endif
