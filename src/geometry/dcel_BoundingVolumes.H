/*!
  @file   dcel_BoundingSphere.H
  @brief  Declaration of a BoundingSphere class. Sphere can be computed with various approximate or exact functions. 
  @author Robert Marskar
  @date   Ap. 2018
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

  class BoundingVolume {
  public:
    
    BoundingVolume();
    ~BoundingVolume();

    PointLocation getPointLocation(const Vec3T<double>& a_p) const;
    
    bool isPointInside(const Vec3T<double>& a_x0) const;
    bool isPointOutside(const Vec3T<double>& a_x0) const;
    
    virtual double getDistanceToPoint(const Vec3T<double>& a_x0) const = 0;
  };
  
  class BoundingSphere : public BoundingVolume {
  public:

    enum class Algorithm {
      Ritter,
    };

    BoundingSphere();
    BoundingSphere(const BoundingSphere& a_other);
    BoundingSphere(const std::vector<Vec3T<double>>& a_points, const Algorithm& a_alg);
    ~BoundingSphere();

    void define(const std::vector<Vec3T<double>>& a_points, const Algorithm& a_algorith);

    bool intersects(const BoundingSphere& a_other);
    
    double getDistanceToPoint(const Vec3T<double>& a_x0) const override;
    
    double& getRadius();
    const double& getRadius() const;
    
    Vec3T<double>& getCenter();
    const Vec3T<double>& getCenter() const;

  protected:

    double m_radius;
    Vec3T<double> m_center;

    void buildRitter(const std::vector<Vec3T<double>>& a_points);
  };

  class AABB : public BoundingVolume {
  public:

    AABB();
    AABB(const AABB& a_other);
    AABB(const std::vector<AABB>& a_others);
    AABB(const std::vector<Vec3T<double>>& a_points);
    ~AABB();

    void define(const std::vector<Vec3T<double>>& a_points);

    bool intersects(const AABB& a_other) const;
    
    double getDistanceToPoint(const Vec3T<double>& a_x0) const override;

    Vec3T<double>& getLowCorner();
    const Vec3T<double>& getLowCorner() const;

    Vec3T<double>& getHighCorner();
    const Vec3T<double>& getHighCorner() const;

  protected:

    Vec3T<double> m_loCorner;
    Vec3T<double> m_hiCorner;
  };
}

#endif
