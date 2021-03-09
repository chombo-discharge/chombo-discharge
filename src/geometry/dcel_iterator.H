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
  class face;

  class edge_iterator {
  public:

    edge_iterator() = delete;
  
    edge_iterator(face& a_face);
    edge_iterator(const face& a_face);

    edge_iterator(vertex& a_vert);
    edge_iterator(const vertex& a_vert);

    std::shared_ptr<edge>& operator() ();

    void reset();

    void operator ++();

    bool ok();
  
  protected:

    enum class IterationMode {
      Vertex,
      Face
    };

    bool m_fullLoop;

    IterationMode m_iterMode;
  
    std::shared_ptr<edge> m_startEdge;   // Starting edge
    std::shared_ptr<edge> m_curEdge; // Current edge
  };
}

#endif
