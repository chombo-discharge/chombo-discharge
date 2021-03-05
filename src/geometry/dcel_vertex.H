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
    vertex(const vertex& a_otherVertex);
    ~vertex();

    void define(const RealVect& m_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal = RealVect::Zero);
    void setPosition(const RealVect& a_pos);
    void setNormal(const RealVect& a_normal);
    void setEdge(const std::shared_ptr<edge>& a_edge);
    void addPolygon(const std::shared_ptr<polygon>& a_poly);

    std::shared_ptr<edge>& getEdge();
    const std::shared_ptr<edge>& getEdge() const;

    std::vector<std::shared_ptr<polygon> > getPolygons();

    std::vector<std::shared_ptr<polygon> >& getPolycache();
    const std::vector<std::shared_ptr<polygon> >& getPolycache() const;

    RealVect& getPosition();
    const RealVect& getPosition() const;

    RealVect& getNormal();
    const RealVect& getNormal() const;
  
    Real signedDistance(const RealVect a_x0) const;

  protected:

    RealVect m_pos;
    RealVect m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<polygon> > m_polycache;
  };
}

#endif
