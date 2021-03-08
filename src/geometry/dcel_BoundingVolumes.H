/*!
  @file   dcel_BoundingSphere.H
  @brief  Declaration of a BoundingSphere class. Sphere can be computed with various approximate or exact functions. 
  @author Robert Marskar
  @date   Ap. 2018
*/

#ifndef _DCEL_BOUNDINGVOLUMES_H_
#define _DCEL_BOUNDINGVOLUMES_H_

#include "RealVect.H"

#include <vector>

namespace dcel {

  enum class PointLocation{
    Inside,
    Outside,
    Boundary
  };

  /*!
    @brief Generic bounding volume class
  */
  class BoundingVolume {
  public:
    
    BoundingVolume();
    ~BoundingVolume();

    PointLocation getPointLocation(const RealVect& a_p) const;
    
    bool isPointInside(const RealVect& a_x0) const;
    bool isPointOutside(const RealVect& a_x0) const;
    
    virtual Real getDistanceToPoint(const RealVect& a_x0) const = 0;
  };
  
  /*!
    @brief Bounding sphere for enclosing a set of points. 
  */
  class BoundingSphere : public BoundingVolume {
  public:

    enum class Algorithm {
      Ritter,
    };

    BoundingSphere();
    BoundingSphere(const BoundingSphere& a_other);
    BoundingSphere(const std::vector<RealVect>& a_points, const Algorithm& a_alg);
    ~BoundingSphere();

    void define(const std::vector<RealVect>& a_points, const Algorithm& a_algorith);

    bool intersects(const BoundingSphere& a_other);
    
    Real getDistanceToPoint(const RealVect& a_x0) const override;
    
    Real& getRadius();
    const Real& getRadius() const;
    
    RealVect& getCenter();
    const RealVect& getCenter() const;

  protected:

    Real m_radius;
    RealVect m_center;

    void buildRitter(const std::vector<RealVect>& a_points);
  };

  /*!
    @brief Axis-aligned bounding box for enclosing a set of points. 
  */
  class AABB : public BoundingVolume{
  public:

    AABB();
    AABB(const AABB& a_other);
    AABB(const std::vector<AABB>& a_others);
    AABB(const std::vector<RealVect>& a_points);
    ~AABB();

    void define(const std::vector<RealVect>& a_points);

    bool intersects(const AABB& a_other) const;
    
    Real getDistanceToPoint(const RealVect& a_x0) const override;

    RealVect& getLowCorner();
    const RealVect& getLowCorner() const;

    RealVect& getHighCorner();
    const RealVect& getHighCorner() const;

  protected:

    RealVect m_loCorner;
    RealVect m_hiCorner;
  };
}

#endif
