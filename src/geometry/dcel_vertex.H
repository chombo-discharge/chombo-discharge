/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "RealVect.H"

#include <vector>
#include <memory>

namespace dcel {

  class edge;
  class polygon;

  class vertex {
  public:

    vertex();
    vertex(const RealVect& a_pos);
    vertex(const RealVect& a_pos, const RealVect& a_normal);
    vertex(const vertex& a_otherVertex);
    ~vertex();

    void define(const RealVect& m_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal = RealVect::Zero) noexcept;
    void setPosition(const RealVect& a_pos) noexcept;
    void setNormal(const RealVect& a_normal) noexcept;
    void setEdge(const std::shared_ptr<edge>& a_edge) noexcept;
    void addPolygonToCache(const std::shared_ptr<polygon>& a_poly) noexcept;
    void normalizeNormalVector() noexcept;
    void clearPolygonCache() noexcept;

    RealVect& getPosition() noexcept;
    const RealVect& getPosition() const noexcept;
    
    std::shared_ptr<edge>& getEdge() noexcept;
    const std::shared_ptr<edge>& getEdge() const noexcept;

    RealVect& getNormal() noexcept;
    const RealVect& getNormal() const noexcept;

    std::vector<std::shared_ptr<polygon> > getPolygons() noexcept;

    std::vector<std::shared_ptr<polygon> >& getPolycache() noexcept;
    const std::vector<std::shared_ptr<polygon> >& getPolycache() const noexcept;
  
    Real signedDistance(const RealVect& a_x0) const noexcept;
    Real unsignedDistance2(const RealVect& a_x0) const noexcept; // Returns squared distance to point

  protected:

    RealVect m_pos;
    RealVect m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<polygon> > m_polycache;
  };
}

#endif
