/*!
  @file   bvh_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BVH_IF_
#define _BVH_IF_

#include "BVH.H"
#include "dcel_face.H"

#include <memory>
#include <chrono>

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"

using namespace dcel;
using namespace BVH;

template <class T, class BV>
class bvh_if : public BaseIF {
public:

  using Vec3  = Vec3T<T>;
  using face  = faceT<T>;
  using Node  = NodeT<T, face, BV>;
  

  //  bvh_if() = delete;
  bvh_if();
  bvh_if(const std::shared_ptr<Node>& a_root, const bool a_flipInside);
  bvh_if(const bvh_if& a_object);
  ~bvh_if();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

protected:
  std::shared_ptr<Node> m_root;

  mutable long m_numCalled;

  mutable std::chrono::duration<double> m_timespan;
  
  bool m_flipInside;
};
#include "CD_NamespaceFooter.H"

#include "bvh_ifI.H"

#endif
