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
    ~vertex();

    void define(const RealVect& m_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal = RealVect::Zero);
    void set_pos(const RealVect& a_pos);
    void set_normal(const RealVect& a_normal);
    void set_edge(const std::shared_ptr<edge>& a_edge);
    void add_polygon(const std::shared_ptr<polygon>& a_poly);

    std::shared_ptr<edge>& get_edge();
    const std::shared_ptr<edge>& get_edge() const;

    std::vector<std::shared_ptr<polygon> > get_polygons();

    std::vector<std::shared_ptr<polygon> >& get_polycache();
    const std::vector<std::shared_ptr<polygon> >& get_polycache() const;

    RealVect& position();
    const RealVect& position() const;

    RealVect& normal();
    const RealVect& normal() const;
  
    Real signed_distance(const RealVect a_x0) const;

  protected:

    RealVect m_pos;
    RealVect m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<polygon> > m_polycache;
  };
}

#endif
