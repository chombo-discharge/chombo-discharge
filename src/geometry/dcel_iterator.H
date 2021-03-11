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
  
  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;

  template <class T>
  class edgeIteratorT {
  public:

    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    edgeIteratorT() = delete;
  
    edgeIteratorT(face& a_face);
    edgeIteratorT(const face& a_face);

    edgeIteratorT(vertex& a_vert);
    edgeIteratorT(const vertex& a_vert);

    inline
    auto operator() () noexcept -> std::shared_ptr<edge>&;

    inline
    auto operator() () const noexcept -> const std::shared_ptr<edge>&;

    inline
    auto reset() noexcept -> void;

    inline
    auto operator ++() noexcept -> void;

    inline
    auto ok() const noexcept -> bool;
  
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

#include "dcel_iteratorI.H"

#endif
