/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelIterator.H
  @brief  Declaration of an iterators for DCEL surface tesselations
  @author Robert Marskar
*/

#ifndef CD_DcelIterator_H
#define CD_DcelIterator_H

// Std includes
#include <memory>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Dcel {
  
  template <class T> class VertexT;
  template <class T> class EdgeT;
  template <class T> class FaceT;

  template <class T>
  class edgeIteratorT {
  public:

    using Vertex = VertexT<T>;
    using Edge   = EdgeT<T>;
    using Face   = FaceT<T>;

    using VertexPtr = std::shared_ptr<Vertex>;
    using EdgePtr   = std::shared_ptr<Edge>;
    using FacePtr   = std::shared_ptr<Face>;

    edgeIteratorT() = delete;
  
    edgeIteratorT(Face& a_face);
    edgeIteratorT(const Face& a_face);

    edgeIteratorT(Vertex& a_vert);
    edgeIteratorT(const Vertex& a_vert);

    inline
    EdgePtr& operator() () noexcept;

    inline
    const EdgePtr& operator() () const noexcept;

    inline
    void reset() noexcept;

    inline
    void operator++() noexcept;

    inline
    bool ok() const noexcept;
  
  protected:

    enum class IterationMode {
      Vertices,
      Faces
    };

    bool m_fullLoop;

    IterationMode m_iterMode;
  
    std::shared_ptr<Edge> m_startEdge;   // Starting edge
    std::shared_ptr<Edge> m_curEdge; // Current edge
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_DcelIteratorImplem.H>

#endif
