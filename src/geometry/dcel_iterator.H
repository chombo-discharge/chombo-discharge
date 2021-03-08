/*!
  @file   dcel_iterator.H
  @brief  Declaration of an iterators for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ITERATOR_
#define _DCEL_ITERATOR_

#include "RealVect.H"

#include <memory>

namespace dcel {
  
  class vertex;
  class edge;
  class polygon;

  class edge_iterator {
  public:

    edge_iterator() = delete;
  
    edge_iterator(polygon& a_poly);
    edge_iterator(const polygon& a_poly);

    edge_iterator(vertex& a_vert);
    edge_iterator(const vertex& a_vert);

    std::shared_ptr<edge>& operator() ();

    void reset();

    void operator ++();

    bool ok();
  
  protected:

    enum class IterationMode {
      Vertex,
      Polygon
    };

    bool m_polyMode;
    bool m_full_loop;

    IterationMode m_iterMode;
  
    std::shared_ptr<edge> m_startEdge;   // Starting edge
    std::shared_ptr<edge> m_curEdge; // Current edge
  };
}

#endif
