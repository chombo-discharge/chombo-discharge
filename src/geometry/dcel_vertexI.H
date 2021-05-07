/*!
  @file   dcel_vertexI.H
  @brief  Implementaion of dcel_vertex.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERTEXI_
#define _DCEL_VERTEXI_

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

namespace ChomboDischarge {
  namespace dcel {

    template <class T>
    inline
    vertexT<T>::vertexT(){
      m_pos    = Vec3::zero();
      m_normal = Vec3::zero();

      m_faces.resize(0);
    }

    template <class T>
    inline
    vertexT<T>::vertexT(const Vec3& a_pos) : vertexT() {
      m_pos    = a_pos;
    }

    template <class T>
    inline
    vertexT<T>::vertexT(const Vec3& a_pos, const Vec3& a_normal) : vertexT() {
      m_pos    = a_pos;
      m_normal = a_normal;
    }

    template <class T>
    inline
    vertexT<T>::vertexT(const vertexT<T>& a_otherVertex){
      m_pos          = a_otherVertex.m_pos;
      m_normal       = a_otherVertex.m_m_normal;
      m_outgoingEdge = a_otherVertex.m_outgoingEdge;
    }

    template <class T>
    inline
    vertexT<T>::~vertexT(){

    }

    template <class T>
    inline
    void vertexT<T>::define(const Vec3& a_pos, const edgePtr& a_edge, const Vec3& a_normal) noexcept {
      m_pos    = a_pos;
      m_outgoingEdge   = a_edge;
      m_normal = a_normal;
    }

    template <class T>
    inline
    void vertexT<T>::setPosition(const Vec3& a_pos) noexcept {
      m_pos = a_pos;
    }

    template <class T>
    inline
    void vertexT<T>::setEdge(const edgePtr& a_edge) noexcept {
      m_outgoingEdge = a_edge;
    }

    template <class T>
    inline
    void vertexT<T>::setNormal(const Vec3& a_normal) noexcept {
      m_normal = a_normal;
    }

    template <class T>
    inline
    void vertexT<T>::addFace(const facePtr& a_face) noexcept {
      m_faces.emplace_back(a_face);
    }

    template <class T>
    inline
    void vertexT<T>::normalizeNormalVector() noexcept {
      m_normal = m_normal/m_normal.length();
    }

    template <class T>
    inline
    void vertexT<T>::computeVertexNormalAverage() noexcept {
      this->computeVertexNormalAverage(m_faces);
    }

    template <class T>
    inline
    void vertexT<T>::computeVertexNormalAverage(const std::vector<facePtr>& a_faces) noexcept {
      m_normal = Vec3::zero();
    
      for (const auto& f : a_faces){
	m_normal += f->getNormal();
      }

      this->normalizeNormalVector();
    }

    template <class T>
    inline
    void vertexT<T>::computeVertexNormalAngleWeighted() noexcept {
      computeVertexNormalAngleWeighted(m_faces);
    }

    template <class T>
    inline
    void vertexT<T>::computeVertexNormalAngleWeighted(const std::vector<facePtr>& a_faces) noexcept {
      m_normal = Vec3::zero();

      // This is more complicated than it looks, and it happens because we want the two half edges that has the current vertex
      // as a mutual vertex (i.e. the "incoming" and "outgoing" edges into this vertex). Normally we'd just iterate through edges,
      // but if it happens that an input face is flipped, this will result in infinite iteration. Instead, we have stored the pointers
      // to each face connected to this vertex. We look through each face to find the endpoints of the edges the have the current vertex
      // as the common vertex, and then compute the subtended angle between those. Sigh...

      const vertexPtr& originVertex = m_outgoingEdge->getVertex(); // AKA 'this'
    
      for (const auto& f : a_faces){
      
	std::vector<vertexPtr> inoutVertices(0);
	for (edgeIterator edgeIt(f->getHalfEdge()); edgeIt.ok(); ++edgeIt){
	  const auto& e = edgeIt();

	  const auto& v1 = e->getVertex();
	  const auto& v2 = e->getOtherVertex();

	  if(v1 == originVertex || v2 == originVertex){
	    if(v1 == originVertex){
	      inoutVertices.emplace_back(v2);
	    }
	    else if(v2 == originVertex){
	      inoutVertices.emplace_back(v1);
	    }
	    else{
	      std::cerr << "In file 'dcel_vertexI.H' function vertexT<T>::computeVertexNormalAngleWeighted() - logic bust.\n";
	    }
	  }
	}

	if(inoutVertices.size() != 2) {
	  std::cerr << "In file 'dcel_vertexI.H' function vertexT<T>::computeVertexNormalAngleWeighted() - logic bust 2.\n";
	}

	const Vec3& x0 = originVertex    ->getPosition();
	const Vec3& x1 = inoutVertices[0]->getPosition();
	const Vec3& x2 = inoutVertices[1]->getPosition();

	if(x0 == x1 || x0 == x2 || x1 == x2){
	  std::cerr << "In file 'dcel_vertexI.H' function vertexT<T>::computeVertexNormalAngleWeighted() - logic bust 3.\n";
	}

	Vec3 v1 = x1-x0;
	Vec3 v2 = x2-x0;

	v1 = v1/v1.length();
	v2 = v2/v2.length();

	const Vec3& norm = f->getNormal();

	const T alpha = acos(v1.dot(v2));

	m_normal += alpha*norm;
      }

      this->normalizeNormalVector();
    }
  
    template <class T>
    inline
    Vec3T<T>& vertexT<T>::getPosition() noexcept {
      return (m_pos);
    }

    template <class T>
    inline
    const Vec3T<T>& vertexT<T>::getPosition() const noexcept {
      return (m_pos);
    }

    template <class T>
    inline
    Vec3T<T>& vertexT<T>::getNormal() noexcept {
      return (m_normal);
    }

    template <class T>
    inline
    const Vec3T<T>& vertexT<T>::getNormal() const noexcept {
      return (m_normal);
    }

    template <class T>
    inline
    std::shared_ptr<edgeT<T> >& vertexT<T>::getOutgoingEdge() noexcept {
      return (m_outgoingEdge);
    }

    template <class T>
    inline
    const std::shared_ptr<edgeT<T> >& vertexT<T>::getOutgoingEdge() const noexcept {
      return (m_outgoingEdge);
    }

    template <class T>
    inline
    std::vector<std::shared_ptr<faceT<T> > >& vertexT<T>::getFaces() noexcept {
      return (m_faces);
    }

    template <class T>
    inline
    const std::vector<std::shared_ptr<faceT<T> > >& vertexT<T>::getFaces() const noexcept {
      return (m_faces);
    }

    template <class T>
    inline
    T vertexT<T>::signedDistance(const Vec3& a_x0) const noexcept {
      const auto delta = a_x0 - m_pos;
      const T dist     = delta.length();
      const T dot      = m_normal.dot(delta);
      const int sign   = (dot > 0.) ? 1 : -1;
  
      return dist*sign;
    }

    template <class T>
    inline
    T vertexT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
      const auto d = a_x0 - m_pos;

      return d.dot(d);
    }
  }
}

#endif
