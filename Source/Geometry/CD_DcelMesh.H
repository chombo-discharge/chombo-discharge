/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelMesh.H
  @brief  Declaration of a mesh class for handling surface Tesselations
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

  template <class T> class VertexT;
  template <class T> class EdgeT;
  template <class T> class FaceT;

  template <class T>
  class MeshT {
  public:

    using Vec3   = Vec3T<T>;
    
    using Vertex = VertexT<T>;
    using Edge   = EdgeT<T>;
    using Face   = FaceT<T>;

    using VertexPtr = std::shared_ptr<Vertex>;
    using EdgePtr   = std::shared_ptr<Edge>;
    using FacePtr   = std::shared_ptr<Face>;
    
    using Mesh   = MeshT<T>;

    MeshT();
    MeshT(const Mesh& a_otherMesh) = delete;
    MeshT(std::vector<FacePtr >&   a_faces,
	  std::vector<EdgePtr >&   a_edges,
	  std::vector<VertexPtr >& a_vertices);

    ~MeshT();

    inline
    void define(std::vector<FacePtr >&   a_faces,
		std::vector<EdgePtr >&   a_edges,
		std::vector<VertexPtr >& a_vertices) noexcept;
    
    inline
    void sanityCheck() const noexcept;

    inline
    void setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept;

    inline
    void setInsideOutsideAlgorithm(const Dcel::InsideOutsideAlgorithm a_algorithm) noexcept;

    inline
    void reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept;

    inline
    std::vector<VertexPtr>& getVertices() noexcept;

    inline
    const std::vector<VertexPtr>& getVertices() const noexcept;

    inline
    std::vector<EdgePtr>& getEdges() noexcept;

    inline
    const std::vector<EdgePtr>& getEdges() const noexcept;

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

    std::vector<VertexPtr >  m_vertices;
    std::vector<EdgePtr >    m_edges;
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
