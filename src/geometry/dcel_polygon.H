/*!
  @file   dcel_poly.H
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_
#define _DCEL_POLY_

#include <ritter_sphere.H>
#include <dcel_BoundingSphere.H>

#include "RealVect.H"

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;

  class polygon {
  public:
  
    polygon();
    polygon(const polygon& a_otherPolygon);
    ~polygon();

    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge);

    void setEdge(const std::shared_ptr<edge>& a_edge);
    void setNormal(const RealVect& a_normal);

    void normalizeNormalVector();
    void computeArea();
    void computeCentroid();
    void computeNormal(const bool a_outwardNormal);
    void computeBoundingBox();
    void computeBoundingSphere();

    std::shared_ptr<edge>& getEdge();
    const std::shared_ptr<edge>& getEdge() const;

    std::vector<RealVect> getAllVertexCoordinates();
    std::vector<std::shared_ptr<vertex> > getVertices();
    std::vector<std::shared_ptr<edge> > getEdges();

    RealVect& getCentroid();
    const RealVect& getCentroid() const;

    RealVect& getNormal();
    const RealVect& getNormal() const;


    Real& getArea();
    const Real& getArea() const;

    RealVect& getBoundingBoxLo();
    const RealVect& getBoundingBoxLo() const;

    RealVect& getBoundingBoxHi();
    const RealVect& getBoundingBoxHi() const;
    
    // Real signed_distance(const RealVect a_x0) const; Should be this signature. 
    Real signedDistance(const RealVect a_x0);

#if 1 // Just so kd will shut up
    RealVect get_coord() const {
      return m_centroid;
    }
#endif
  
  protected:
  
    std::shared_ptr<edge> m_edge;

    Real m_area;
    
    RealVect m_normal;
    RealVect m_centroid;
    RealVect m_lo;
    RealVect m_hi;

    BoundingSphere m_boundingSphere;
    ritter_sphere m_ritter;
  };

}

#endif
