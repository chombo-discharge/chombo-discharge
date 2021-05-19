/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelMesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
*/

#ifndef CD_DcelMesh_H
#define CD_DcelMesh_H

// Std includes
#include <vector>
#include <memory>
#include <functional>
#include <map>

// Our includes
#include <CD_DcelAlgorithms.H>
#include "CD_NamespaceHeader.H"

namespace Dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class FaceT;

  template <class T>
  class MeshT {
  public:

    using Vec3   = Vec3T<T>;
    
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using Face   = FaceT<T>;

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using FacePtr   = std::shared_ptr<Face>;
    
    using Mesh   = MeshT<T>;

    MeshT();
    MeshT(const Mesh& a_otherMesh) = delete;
    MeshT(std::vector<FacePtr >&   a_faces,
	  std::vector<edgePtr >&   a_edges,
	  std::vector<vertexPtr >& a_vertices);

    ~MeshT();

    inline
    void define(std::vector<FacePtr >&   a_faces,
		std::vector<edgePtr >&   a_edges,
		std::vector<vertexPtr >& a_vertices) noexcept;
    
    inline
    void sanityCheck() const noexcept;

    inline
    void setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept;

    inline
    void setInsideOutsideAlgorithm(const Dcel::InsideOutsideAlgorithm a_algorithm) noexcept;

    inline
    void reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept;

    inline
    std::vector<vertexPtr>& getVertices() noexcept;

    inline
    const std::vector<vertexPtr>& getVertices() const noexcept;

    inline
    std::vector<edgePtr>& getEdges() noexcept;

    inline
    const std::vector<edgePtr>& getEdges() const noexcept;

    inline
    std::vector<FacePtr>& getFaces() noexcept;

    inline
    const std::vector<FacePtr>& getFaces() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T signedDistance(const Vec3& a_x0, SearchAlgorithm a_algorithm) const noexcept;

  protected:

    SearchAlgorithm m_algorithm;

    std::vector<vertexPtr >  m_vertices;
    std::vector<edgePtr >    m_edges;
    std::vector<FacePtr >    m_faces;

    inline
    std::vector<Vec3T<T> > getAllVertexCoordinates() const noexcept;

    inline
    void reconcileFaces() noexcept;

    inline
    void reconcileEdges() noexcept;

    inline
    void reconcileVertices(VertexNormalWeight a_weight) noexcept;

    inline
    T DirectSignedDistance(const Vec3& a_point) const noexcept;

    inline
    T DirectSignedDistance2(const Vec3& a_point) const noexcept;

    inline
    void incrementWarning(std::map<std::string, int>& a_warnings, const std::string& a_warn) const noexcept;

    inline
    void printWarnings(const std::map<std::string, int>& a_warnings) const noexcept;
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_DcelMeshImplem.H>

#endif
