/*!
  @file   dcel_BVH.H
  @brief  Declaration of boundary volume hierarchies for use with dcel meshes.
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

#include "dcel_vec.H"

#include <memory>
#include <vector>
#include <functional>

namespace dcel {
  namespace BVH {
        
    template <class T, class BV>
    class Node {
    public:

      Node();
      Node(std::shared_ptr<Node<T, BV> >& a_parent);
      ~Node();

      inline
      void setParent(const std::shared_ptr<Node<T, BV> >& a_parent) noexcept;

      inline
      void setLeft(const std::shared_ptr<Node<T, BV> >& a_left) noexcept;

      inline
      void setRight(const std::shared_ptr<Node<T, BV> >& a_right) noexcept;

      inline
      PointLocation getPointLocation(const Vec3T<double>& a_p) const noexcept;

      inline
      bool isRoot() const noexcept;

      inline
      bool isLeaf() const noexcept;

      inline
      bool isPointInside(const Vec3T<double>& a_x0) const noexcept;

      inline
      bool isPointOutside(const Vec3T<double>& a_x0) const noexcept;

      inline
      double getDistanceToPoint(const Vec3T<double>& a_x0) const noexcept;

      inline
      BV& getBoundingVolume() noexcept;

      inline
      const BV& getBoundingVolume() const noexcept;

      inline
      std::vector<std::shared_ptr<T> >& getObjects() noexcept;

      inline
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
    public:
      Tree();
      ~Tree();

      
    protected:

      std::shared_ptr<Node<T, BV> > m_root;
    };
  }
}

#include "dcel_BVHI.H"

#endif
