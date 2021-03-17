/*!
  @file   bvh_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BVH_IF_
#define _BVH_IF_

#include "dcel_BVH.H"
#include "dcel_mesh.H"

#include <memory>

#include <BaseIF.H>

using namespace dcel;
using namespace dcel::BVH;

template <class T, class BV>
class bvh_if : public BaseIF {
public:

  using Vec3  = Vec3T<T>;
  using mesh  = meshT<T>;
  using face  = faceT<T>;
  using Node  = NodeT<T, face, BV>;

  // Lambdas for interfacing into dcel_BVH
  static StopFunctionT<T, face, BV> defaultStopFunction;
  static PartitionFunctionT<face>   defaultPartitionFunction;
  static BVConstructorT<face, BV>   defaultBVConstructor;



  bvh_if() = delete;

  bvh_if(const std::shared_ptr<Node>& a_root, const bool a_flipInside);
  bvh_if(const bvh_if& a_object);
  ~bvh_if();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

protected:
  std::shared_ptr<Node> m_root;
  
  bool m_flipInside;
};

#include "bvh_ifI.H"

#endif
