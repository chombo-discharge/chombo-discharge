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
    class NodeT {
    public:

      using Vec3 = Vec3T<T>;
      using Node = NodeT<T, BV>;

      NodeT();
      NodeT(std::shared_ptr<Node>& a_parent);
      NodeT(std::vector<std::shared_ptr<T> >& a_data);
      ~NodeT();

      inline
      decltype(auto) setParent(const std::shared_ptr<Node>& a_parent) noexcept;

      inline
      decltype(auto) setLeft(const std::shared_ptr<Node>& a_left) noexcept;

      inline
      decltype(auto) setRight(const std::shared_ptr<Node>& a_right) noexcept;

      inline
      decltype(auto) getParent() noexcept;

      inline
      decltype(auto) getParent() const noexcept;

      inline
      decltype(auto) isRoot() const noexcept;

      inline
      decltype(auto) isLeaf() const noexcept;

      inline
      decltype(auto) isPointInside(const Vec3& a_x0) const noexcept;

      inline
      decltype(auto) isPointOutside(const Vec3& a_x0) const noexcept;

      inline
      decltype(auto) getDistanceToPoint(const Vec3& a_x0) const noexcept;

      inline
      decltype(auto) getBoundingVolume() noexcept;

      inline
      decltype(auto) getBoundingVolume() const noexcept;

      inline
      decltype(auto) getObjects() noexcept;

      inline
      decltype(auto) getObjects() const noexcept;

    protected:

      BV m_bv;

      std::vector<std::shared_ptr<T> > m_objects;

      std::shared_ptr<Node> m_parent;
      std::shared_ptr<Node> m_left;
      std::shared_ptr<Node> m_right;
    };

    template <class T, class BV>
    class Tree {
    public:

      using Node = NodeT<T, BV>;
      
      Tree();
      Tree(std::vector<std::shared_ptr<T> >& a_objects);
      ~Tree();

      
    protected:

      std::shared_ptr<Node> m_root;
    };
  }
}

#include "dcel_BVHI.H"

#endif
