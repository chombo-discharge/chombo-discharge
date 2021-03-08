/*!
  @file   dcel_BVH.H
  @brief  Declaration of boundary volume hierarchies for use with dcel meshes.
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

#include <memory>
#include <vector>

namespace dcel {
  namespace BVH {
        
    template <class T, class BV>
    class Node {
    public:

      Node();
      Node(std::shared_ptr<Node<T, BV> >& a_parent);
      ~Node();

      void setParent(const std::shared_ptr<Node<T, BV> >& a_parent) noexcept;
      void setLeft(const std::shared_ptr<Node<T, BV> >& a_left) noexcept;
      void setRight(const std::shared_ptr<Node<T, BV> >& a_right) noexcept;

      PointLocation getPointLocation(const RealVect& a_p) const noexcept;

      bool isRoot() const noexcept;
      bool isLeaf() const noexcept;
    
      bool isPointInside(const RealVect& a_x0) const noexcept;
      bool isPointOutside(const RealVect& a_x0) const noexcept;
      
      Real getDistanceToPoint(const RealVect& a_x0) const noexcept;

      BV& getBoundingVolume() noexcept;
      const BV& getBoundingVolume() const noexcept;

      std::vector<std::shared_ptr<T> >& getObjects() noexcept;
      const std::vector<std::shared_ptr<T> >& getObjects() const noexcept;

    protected:

      BV m_bv;

      std::vector<std::shared_ptr<T> > m_objects;

      std::shared_ptr<Node<T, BV> > m_parent;
      std::shared_ptr<Node<T, BV> > m_left;
      std::shared_ptr<Node<T, BV> > m_right;
    };

    template <class T, class BV>
    class Tree {
      Tree();
      ~Tree();

    protected:

      std::shared_ptr<Node<T, BV> > m_root;
    };
  }
}

#include "dcel_BVHI.H"

#endif
