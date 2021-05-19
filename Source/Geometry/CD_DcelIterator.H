/*!
  @file   dcel_iterator.H
  @brief  Declaration of an iterators for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ITERATOR_
#define _DCEL_ITERATOR_

#include <memory>

#include "CD_NamespaceHeader.H"
namespace Dcel {
  
  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;

  template <class T>
  class edgeIteratorT {
  public:

    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using facePtr   = std::shared_ptr<face>;

    edgeIteratorT() = delete;
  
    edgeIteratorT(face& a_face);
    edgeIteratorT(const face& a_face);

    edgeIteratorT(vertex& a_vert);
    edgeIteratorT(const vertex& a_vert);

    inline
    edgePtr& operator() () noexcept;

    inline
    const edgePtr& operator() () const noexcept;

    inline
    void reset() noexcept;

    inline
    void operator++() noexcept;

    inline
    bool ok() const noexcept;
  
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
#include "CD_NamespaceFooter.H"

#include <CD_DcelIteratorImplem.H>

#endif
