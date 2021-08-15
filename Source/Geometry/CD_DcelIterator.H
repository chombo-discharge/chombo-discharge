/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DcelIterator.H
  @brief  Declaration of an iterators for DCEL surface Tesselations
  @author Robert Marskar
*/

#ifndef CD_DcelIterator_H
#define CD_DcelIterator_H

// Std includes
#include <memory>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Dcel {

  // Forward declare classes.
  template <class T> class VertexT;
  template <class T> class EdgeT;
  template <class T> class FaceT;

  /*!
    @brief Class which can iterate through edges and vertices around a DCEL polygon face.
    @details This class can be used so that it either visits all the half-edges in a face, or all the outgoing half-edges from a vertex. 
  */
  template <class T>
  class edgeIteratorT {
  public:

    /*!
      @brief Alias to cut down on typing
    */
    using Vertex = VertexT<T>;

    /*!
      @brief Alias to cut down on typing
    */    
    using Edge = EdgeT<T>;

    /*!
      @brief Alias to cut down on typing
    */    
    using Face = FaceT<T>;

    /*!
      @brief Alias to cut down on typing
    */    
    using VertexPtr = std::shared_ptr<Vertex>;
    
    /*!
      @brief Alias to cut down on typing
    */        
    using EdgePtr = std::shared_ptr<Edge>;

    /*!
      @brief Alias to cut down on typing
    */        
    using FacePtr = std::shared_ptr<Face>;

    /*!
      @brief Default construction is not allowed. Use one of the full constructors
    */
    edgeIteratorT() = delete;

    /*!
      @brief Constructor, taking a face as argument. The iterator begins at the half-edge pointer contained in the face
      @param[in] a_face DCEL polygon face
      @note This constructor will will iterate through the half-edges in the polygon face. 
    */
    edgeIteratorT(Face& a_face);

    /*!
      @brief Constructor, taking a face as argument. The iterator begins at the half-edge pointer contained in the face
      @param[in] a_face DCEL polygon face
      @note This constructor will will iterate through the half-edges in the polygon face. 
    */    
    edgeIteratorT(const Face& a_face);

    /*!
      @brief Constructor, taking a vertex as argument. The iterator begins at the outgoing half-edge from the vertex
      @param[in] a_vertex DCEL vertex
      @note This constructor will will iterate through the outgoing half-edges from a vertex. 
    */        
    edgeIteratorT(Vertex& a_vertex);

    /*!
      @brief Constructor, taking a vertex as argument. The iterator begins at the outgoing half-edge from the vertex
      @param[in] a_vertex DCEL vertex
      @note This constructor will will iterate through the outgoing half-edges from a vertex. 
    */            
    edgeIteratorT(const Vertex& a_vertex);

    /*!
      @brief Operator returning a pointer to the current half-edge
    */
    inline
    EdgePtr& operator() () noexcept;

    /*!
      @brief Operator returning a pointer to the current half-edge
    */
    inline
    const EdgePtr& operator() () const noexcept;

    /*!
      @brief Reset function for the iterator. This resets the iterator so that it begins from the starting half-edge
    */
    inline
    void reset() noexcept;

    /*!
      @brief Incrementation operator, bringing the iterator to the next half-edge
    */ 
    inline
    void operator++() noexcept;

    /*!
      @brief Function which checks if the iteration can be continued. 
      @return Returns true unless the current half-edge is a nullptr (i.e., a broken polygon face) OR a full loop has been made around the polygon face (i.e. all half-edges
      have been visited)
    */
    inline
    bool ok() const noexcept;
  
  protected:

    /*!
      @brief Iteration mode, used to distinguish between the two constructors (face- or vertex-based iteration)
    */
    enum class IterationMode {
      Vertices,
      Faces
    };

    /*!
      @brief If true, a full loop has been made around the polygon face
    */
    bool m_fullLoop;

    /*!
      @brief Iteration mode. Set in constructor
    */
    IterationMode m_iterMode;

    /*!
      @brief Starting half-edge
    */
    std::shared_ptr<Edge> m_startEdge;

    /*!
      @brief Current half-edge
    */
    std::shared_ptr<Edge> m_curEdge; 
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_DcelIteratorImplem.H>

#endif
