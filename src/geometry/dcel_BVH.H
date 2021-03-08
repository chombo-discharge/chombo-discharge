/*!
  @file   dcel_BVH.H
  @brief  Declaration of boundary volume hierarchies for use with dcel meshes.
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

namespace dcel {
  namespace BVH {
        
    template <class BV>
    class Node {
    public:

      Node();
      Node(std::shared_ptr<Node<BV> >& a_parent);
      ~Node();

      void setParent(const std::shared_ptr<Node<BV> >& a_parent) noexcept;
      void setLeft(const std::shared_ptr<Node<BV> >& a_left) noexcept;
      void setRight(const std::shared_ptr<Node<BV> >& a_right) noexcept;

      PointLocation getPointLocation(const RealVect& a_p) const noexcept;

      bool isRoot() const noexcept;
      bool isLeaf() const noexcept;
    
      bool isPointInside(const RealVect& a_x0) const noexcept;
      bool isPointOutside(const RealVect& a_x0) const noexcept;
      
      Real getDistanceToPoint(const RealVect& a_x0) const noexcept;

      BV& getBoundingVolume() noexcept;
      const BV& getBoundingVolume() const noexcept;

    protected:

      BV m_bv;

      std::shared_ptr<Node<BV> > m_parent;
      std::shared_ptr<Node<BV> > m_left;
      std::shared_ptr<Node<BV> > m_right;
    };

    template <class BV>
    class Tree {
      Tree();
      ~Tree();

    protected:

      std::shared_ptr<Node<BV> > m_root;
    };
  }
}

#include "dcel_BVHI.H"

#endif
