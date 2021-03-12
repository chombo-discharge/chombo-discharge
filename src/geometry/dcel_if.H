/*!
  @file   dcel_if.H
  @brief  Declaration of an implicit-function class that gets its value function from a DCEL surface tesselation
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_IF_
#define _DCEL_IF_

#include "dcel_BVH.H"
#include "dcel_mesh.H"

#include <memory>

#include <BaseIF.H>

using namespace dcel;


/*!
  @brief Class that gets it's value function from a surface tesselation 
*/
template <class T, class BV>
class dcel_if : public BaseIF {
public:

  using Vec3 = dcel::Vec3T<T>;
  using mesh = dcel::meshT<T>;
  using face = dcel::faceT<T>;
  using AABB = dcel::AABBT<T>;
  using Node = dcel::BVH::NodeT<face, BV>;

  static SortingFunctor<face> defaultSortingFunctor;
  static CostFunctor<face> defaultCostFunctor;

  dcel_if() = delete;

  dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside);
  dcel_if(const dcel_if& a_object);
  ~dcel_if();

  Real value(const RealVect& a_point) const override;
  
  BaseIF* newImplicitFunction() const override;

  inline
  decltype(auto) buildBVH(SortingFunctor<face> sortFunc = defaultSortingFunctor, CostFunctor<face> costFunc = defaultCostFunctor) noexcept;

protected:

  std::shared_ptr<mesh> m_mesh;
  std::shared_ptr<Node> m_root;
  
  bool m_flipInside;
};

#include "dcel_ifI.H"

#endif
